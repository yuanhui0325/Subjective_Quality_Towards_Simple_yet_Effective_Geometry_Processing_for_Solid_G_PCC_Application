/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2019, ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the ISO/IEC nor the names of its contributors
 *   may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <stdint.h>
#include <vector>

#include "entropy.h"
#include "hls.h"
#include "PCCTMC3Common.h"
#include "motionWip.h"

namespace pcc {

//============================================================================

std::vector<int> sortedPointCloud(
  const int attribCount,
  const PCCPointSet3& pointCloud,
  std::vector<int64_t>& mortonCode,
  std::vector<attr_t>& attributes,
  bool isEncoder = false);

//---------------------------------------------------------------------------

void sortedPointCloud(
  const int attribCount,
  const PCCPointSet3& pointCloud,
  const std::vector<int>& indexOrd,
  std::vector<attr_t>& attributes);

//============================================================================

struct MSOctree;
struct EncoderParams;

struct AttributeInterPredParams: InterPredParams {
  bool enableAttrInterPred;
#if 0
  int getPointCount() const { return referencePointCloud.getPointCount(); }
  void clear() { referencePointCloud.clear(); }
#endif
  bool hasLocalMotion() const { return compensatedPointCloud.getPointCount() > 0; }

  void findMotion(
    const EncoderParams* params,
    const EncodeMotionSearchParams& msParams,
    const ParameterSetMotion& mvPS,
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    MotionEntropyEncoder& motionEncoder,
    PCCPointSet3& pointCloud
  );

  int64_t encodeMotionAndBuildCompensated(
    const AttributeParameterSet& aps,
    MotionEntropyEncoder& motionEncoder,
    uint64_t& ioSizeMotionBits
  );

  int64_t buildCompensatedSlabBlock(
    const AttributeParameterSet& aps,
    const ParameterSetMotion& mvPS
  );

  void prepareEncodeMotion(
    const ParameterSetMotion& mvPS,
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    PCCPointSet3& pointCloud,
    const point_t& origin
  );

  void prepareDecodeMotion(
    const ParameterSetMotion& mvPS,
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    PCCPointSet3& pointCloud,
    const point_t& origin
  );

  int64_t decodeMotionAndBuildCompensated(
    const AttributeParameterSet& aps,
    MotionEntropyDecoder& motionDecoder
  );

  void extractMotionForSlabBlock(point_t start, point_t size)
  {
    slabBlockStart = start;
    slabBlockSize = size;
    dualMotion = MVField(
        mvField,
        start,
        start + size);
  }

  void copyMotion()
  {
    dualMotion = mvField;
  }

  void setSlabBlockStart(const point_t& start) { slabBlockStart = start; }
  void setSlabBlockSize(const point_t& size) { slabBlockSize = size; }

  const point_t& getSlabBlockStart() const { return slabBlockStart; }
  const point_t& getSlabBlockSize() const { return slabBlockSize; }

protected:
  point_t slabBlockStart = -1;
  point_t slabBlockSize = -1;
  MVField  dualMotion;
  MSOctree mSOctreeCurr;
  std::vector<std::pair<int/*dualMotion's puNodeIdx*/, int/*mSOctreeCurr's nodeIdx*/> > motionPUTrees;
};

//============================================================================

class AttributeContexts {
public:
  void reset();

protected:
  AdaptiveBitModel ctxZeroBlock[72][2];
  AdaptiveBitModel ctxZeroCoeffs[144][2];
  AdaptiveBitModel ctxCoeffGtN[2][7][6];
  AdaptiveBitModel ctxCoeffRemPrefix[2][12];
  AdaptiveBitModel ctxCoeffRemSuffix[2][12];
  AdaptiveBitModel ctxSkip[4];
};

//----------------------------------------------------------------------------

inline void
AttributeContexts::reset()
{
  this->~AttributeContexts();
  new (this) AttributeContexts;
}

//============================================================================

static int constexpr kNumCtxMode = 7 * 2 * 2 * 4;

class PredModeContexts {
public:
  void reset();

  std::array<AdaptiveBitModel, kNumCtxMode> modeIsIntra;
};

//----------------------------------------------------------------------------

inline void
PredModeContexts::reset()
{
  this->~PredModeContexts();
  new (this) PredModeContexts;
}

//============================================================================

namespace attr {

//============================================================================

enum PredMode:int8_t
{
  Null = 0,
  Intra = 1,
  Inter = 2
};

//============================================================================

template <typename T>
inline bool isNull(T mode) { return mode == PredMode::Null; }
template <typename T>
inline bool isIntra(T mode) { return mode == PredMode::Intra; }
template <typename T>
inline bool isInter(T mode) { return mode == PredMode::Inter; }

//============================================================================

static inline
int getPredCtxMode(
  int childrenCount,
  int cousinPred,
  int cousinPredW,
  int extreme)
{
  int ctx = (childrenCount > 4) + 2 * cousinPred;
  ctx = 2 * ctx + (cousinPredW > 5);
  ctx = 4 * ctx + extreme;

  return ctx;
}

//============================================================================

} // namespace attr

//============================================================================

//============================================================================

namespace RAHT {

struct BlockBoundaries;
struct BlockRefBoundaries;

//============================================================================
// For Raster Scan Ordered RAHT
struct RSO_RAHT {
  //static constexpr int numBitsPerDim = 21;
  static constexpr int numBitsPerDim = 16; // For debugging
  //static constexpr int64_t maxPosXYZ = (1ULL << 3 * numBitsPerDim) - 1;
  static constexpr int64_t maxPosXYZ = std::numeric_limits<int64_t>::max();
  static constexpr int64_t maxPosYZ = (1LL << 2 * numBitsPerDim) - 1;
  static constexpr int64_t maxPosZ = (1LL << numBitsPerDim) - 1;

  static constexpr int64_t maskZ = maxPosZ;
  static constexpr int64_t maskY = maxPosYZ - maxPosZ;
  // x keeps sign bit for debugging
  static constexpr int64_t maskX = /*maxPosXYZ*/ -1 - maxPosYZ;

  static constexpr int64_t maskXY = maskX + maskY;

  // do not apply mask,
  // - values are expected to fit
  // - it propagates sign such that negative offsets will decrease higher order
  //   coordinate by one, ensuring ordering properties
  static constexpr int64_t packX(int x) { return int64_t(x) << 2 * numBitsPerDim; }
  static constexpr int64_t packY(int y) { return int64_t(y) << numBitsPerDim; }
  static constexpr int64_t packZ(int z) { return int64_t(z); }

  static int X(int64_t pos) { return (pos & maskX) >> 2 * numBitsPerDim; }
  static int Y(int64_t pos) { return (pos & maskY) >> numBitsPerDim; }
  static int Z(int64_t pos) { return (pos & maskZ); }

  static Vec3<int> pos3D(int64_t pos)
    { return Vec3<int>{X(pos), Y(pos), Z(pos)}; }

  // addition is used for sign propagation
  static int64_t pack(const Vec3<int>& pos)
    { return packX(pos[0]) + packY(pos[1]) + packZ(pos[2]); }

  static constexpr int64_t pack(int x, int y, int z)
    { return packX(x) + packY(y) + packZ(z); }

  static int64_t pos(const Vec3<int>& pos)
    { return packX(pos[0]) + packY(pos[1]) + packZ(pos[2]); }
};

//============================================================================
// Compute mask of boundaries belongings
struct BlockBoundariesCheckerRSO {

  BlockBoundariesCheckerRSO(
    point_t blockStart,
    point_t blockSizeMinusOne,
    int level
  ) : level(level)
    , lowBoundLevel(RSO_RAHT::pos(blockStart >> level))
    , highBoundLevel(RSO_RAHT::pos(blockStart + blockSizeMinusOne >> level))
    , highBoundChildLevel(RSO_RAHT::pos(blockStart + blockSizeMinusOne >> level - 1))
  {
  }

  int
  computeBoundariesPosMask(int64_t pos) const
  {
    // current position (discard extra precision)
    const int64_t posLevel = RSO_RAHT::pos(RSO_RAHT::pos3D(pos) >> level);
    const auto lowCheck = posLevel ^ lowBoundLevel;
    const auto highCheck = posLevel ^ highBoundLevel;
    int res = !(lowCheck & RSO_RAHT::maskZ)
      | !(lowCheck & RSO_RAHT::maskY) << 1
      | !(lowCheck & RSO_RAHT::maskX) << 2
      | !(highCheck & RSO_RAHT::maskZ) << 3
      | !(highCheck & RSO_RAHT::maskY) << 4
      | !(highCheck & RSO_RAHT::maskX) << 5;
    return res;
  }

  const int level;
  const uint64_t lowBoundLevel;
  const uint64_t highBoundLevel;
  const uint64_t highBoundChildLevel;
};


static constexpr int64_t mask_pos_boundary_rso[3] = {
  RSO_RAHT::maskX | RSO_RAHT::maskY,
  RSO_RAHT::maskX | RSO_RAHT::maskZ,
  RSO_RAHT::maskY | RSO_RAHT::maskZ
};

struct RahtBoundaryNode {
  RahtBoundaryNode(int64_t pos)
    : pos(pos)//, firstChildIdx(0xFFFFFFFF)
    , childsWintra(128)
    {}

  int64_t pos; // packed position of a node (see RSO_RAHT)
  int32_t reconstructedAttr[3]; //reconstructed  attribute
  uint32_t occupancy :8; // 8-bit occupancy
  uint32_t CntNonZero :4; // nb of non zero coeffs, as in node
  uint32_t mode :2; // [0..2]
  uint32_t child_mode :2; // [0..2]
  uint32_t childsWintra :8; // [0..255]
};

struct BlockRefBoundaries {
  // Nodes belonging to the {ZTop, YRight, XFront} faces, per successive layer
  const std::vector<RahtBoundaryNode>* nodes[3] = {};
  // layerEnd[axis][i] is the index after the last node belonging to layer i
  const std::vector<int>* layerEnd[3] = {};

  decltype(nodes[0]->cbegin())
    layer_cbegin(int faceIdx, int layer) const
    {
      return nodes[faceIdx]
        ? nodes[faceIdx]->cbegin() + ((layer + 1 < layerEnd[faceIdx]->size()) ? (*layerEnd[faceIdx])[layer + 1] : 0)
        : decltype(nodes[0]->cbegin())();
    }

  decltype(nodes[0]->cbegin())
    layer_cend(int faceIdx, int layer) const
    {
      return nodes[faceIdx]
        ? nodes[faceIdx]->cbegin() + (*layerEnd[faceIdx])[layer]
        : decltype(nodes[0]->cbegin())();
    }
};

struct BlockBoundaries {
  // Nodes belonging to the {ZTop, YRight, XFront} faces, per successive layer
  std::vector<RahtBoundaryNode> nodes[3];
  // layerEnd[axis][i] is the index after the last node belonging to layer i
  std::vector<int> layerEnd[3];

  // temporary buffer for keeping nodes in lexicographic order
  std::vector<RahtBoundaryNode> nodesNextRow[3];

  void init(point_t slabBlockSize, int numAttrs) {
    // worst case
    int maxBlockSize = std::max(slabBlockSize[2], std::max(slabBlockSize[1], slabBlockSize[0]));
    int maxSize = 0;
    int numLayers = ilog2(uint32_t(maxBlockSize - 1)) + 2;
    for (int d = 0; d < numLayers; ++d) {
      maxSize += maxBlockSize * maxBlockSize; //(1 << d) * (1 << d);
      maxBlockSize = maxBlockSize + (maxBlockSize & 1) >> 1;
    }

    for (int k = 0; k < 3; ++k) {
      nodes[k].clear();
      layerEnd[k].clear();
      nodesNextRow[k].clear();

      nodes[k].reserve(maxSize);
      layerEnd[k].resize(numLayers, 0);
      nodesNextRow[k].reserve(slabBlockSize[k]);
    }
  }

  void flushNextRow(int faceIdx)
  {
    if (nodesNextRow[faceIdx].size()) {
      nodes[faceIdx].insert(
        nodes[faceIdx].end(), nodesNextRow[faceIdx].begin(), nodesNextRow[faceIdx].end());
      nodesNextRow[faceIdx].clear();
    }
  }
};

//============================================================================
// Initialize boundaries to be filled and retrieve boundaries to be used as
// reference

inline
BlockBoundaries& prepareBoundariesInfo(
  point_t slabBlockStart,
  point_t slabBlockSize,
  int numAttr,
  ringbuf<std::pair<point_t/*slab block index*/,RAHT::BlockBoundaries>>& allBoundaries,
  BlockRefBoundaries& boundariesToBeUsed)
{
  boundariesToBeUsed = BlockRefBoundaries(); // clear

  // HACK: allocate if needed for a new element, to not invalidate the pointers
  // TODO: reserve a sufficient buffer size in an upper level function
  allBoundaries.emplace_back();
  allBoundaries.pop_back();

  auto cmp = [&](decltype(*allBoundaries.begin()) a, point_t b) {
    return (a.first[0] < b[0]
      || a.first[0] == b[0] && (
        a.first[1] < b[1]
        || a.first[1] == b[1] && a.first[2] < b[2]));
  };

  // TODO: search may be simplified if critical (keep iterators)
  // as they are ordered in raster scan order

  // x-1
  auto xMinus1 = slabBlockStart - point_t(slabBlockSize[0], 0, 0);
  auto itFoundX = std::lower_bound(
    allBoundaries.begin(), allBoundaries.end(), xMinus1, cmp);
  if (itFoundX != allBoundaries.end() && itFoundX->first == xMinus1) {
    auto& foundBounds = itFoundX->second;
    boundariesToBeUsed.nodes[2] = &foundBounds.nodes[2];
    boundariesToBeUsed.layerEnd[2] = &foundBounds.layerEnd[2];
  }

  // clear expired elements
  while(allBoundaries.begin() != allBoundaries.end()
      && cmp(allBoundaries.front(), xMinus1))
    allBoundaries.pop_front();

  // y-1
  auto yMinus1 = slabBlockStart - point_t(0, slabBlockSize[1], 0);
  auto itFoundY = std::lower_bound(
    allBoundaries.begin(), allBoundaries.end(), yMinus1, cmp);
  if (itFoundY != allBoundaries.end() && itFoundY->first == yMinus1) {
    auto& foundBounds = itFoundY->second;
    boundariesToBeUsed.nodes[1] = &foundBounds.nodes[1];
    boundariesToBeUsed.layerEnd[1] = &foundBounds.layerEnd[1];
  }

  // z-1
  auto zMinus1 = slabBlockStart - point_t(0, 0, slabBlockSize[2]);
  auto itFoundZ = std::lower_bound(
    allBoundaries.begin(), allBoundaries.end(), zMinus1, cmp);
  if (itFoundZ != allBoundaries.end() && itFoundZ->first == zMinus1) {
    auto& foundBounds = itFoundZ->second;
    boundariesToBeUsed.nodes[0] = &foundBounds.nodes[0];
    boundariesToBeUsed.layerEnd[0] = &foundBounds.layerEnd[0];
  }

  // add new element
  // +1 because capacity maybe uint max if empty
  // buf must have a sufficient size otherwise the pointers in boundariesToBeUsed
  // becomes invalidated
  assert(allBoundaries.capacity() + 1 > allBoundaries.size() + 1);
  allBoundaries.emplace_back(slabBlockStart, BlockBoundaries());
  auto& boundariesToBeSet = allBoundaries.back().second;
  // allocate and initialize
  boundariesToBeSet.init(slabBlockSize, numAttr);

  return boundariesToBeSet;
}

//============================================================================

struct UrahtNodeRSOEncoder {
  int64_t pos;

  int64_t sumAttr[3];
  int64_t sumAttrInter[3];

  int32_t weight;

  int32_t reconstructedAttr[3];
  int64_t reconstructedRahtDC[3];

  std::array<int16_t, 2> qp;

  bool decoded = false;
  attr::PredMode mode = attr::PredMode::Null;
  attr::PredMode child_mode = attr::PredMode::Null;
  uint8_t Wintra = 128;
  uint8_t childsWintra = 128;

  uint8_t occupancy = 0;

  uint8_t CntNonZero = 0;
  uint8_t infoParent = 0;
};

struct UrahtNodeDecoder {
  int64_t pos;

  // 32 bits is not enough for 16 bits internal precision or we shall reduce precision of SumAttr*
  int64_t sumAttrInter[3] = { 0, 0, 0 };

  int32_t weight;

  int32_t reconstructedAttr[3];
  int64_t reconstructedRahtDC[3];

  std::array<int16_t, 2> qp;

  bool decoded = false;
  attr::PredMode mode = attr::PredMode::Null;
  attr::PredMode child_mode = attr::PredMode::Null;
  uint8_t Wintra = 128;
  uint8_t childsWintra = 128;

  uint8_t occupancy = 0;

  uint8_t CntNonZero = 0;
  uint8_t infoParent = 0;
};

constexpr int kSizeUrahtNodeDecoder = sizeof(UrahtNodeDecoder);

struct UrahtNodeDecoderHaar {
  int64_t pos;

  int32_t reconstructedAttr[3];
  int16_t sumAttrInter[3] = { 0,0,0 };

  bool decoded = false;
  attr::PredMode mode = attr::PredMode::Null;
  attr::PredMode child_mode = attr::PredMode::Null;
  uint8_t Wintra = 128;
  uint8_t childsWintra = 128;

  uint8_t occupancy = 0;

  uint8_t CntNonZero = 0;
  uint8_t infoParent = 0;

  // hack for template functions compatibility when non existing members are
  // referred conditionally (if constexpr could be used instead, if C++17
  // was allowed)
  static int weight;
  static std::array<int16_t, 2> qp;
  static int64_t reconstructedRahtDC[3];
};

constexpr int kSizeUrahtNodeDecoderHaar = sizeof(UrahtNodeDecoderHaar);

//============================================================================

// Fixed Point precision for RAHT
static constexpr int kFPFracBits = 15;
static constexpr int64_t kFPOne = 1 << kFPFracBits;
static constexpr int64_t kFPOneHalf = 1 << kFPFracBits - 1;
static constexpr int64_t kFPDecMask = ((1 << kFPFracBits) - 1);
static constexpr int64_t kFPIntMask = ~kFPDecMask;

static_assert(kFPFracBits >= 0 && kFPFracBits < 16);

//============================================================================
// Encapsulation of a RAHT transform stage.

class RahtKernel {
public:
  inline
  RahtKernel(int64_t weightLeft, int64_t weightRight, bool scaledWeights = false)
  {
    if (scaledWeights) {
      _a = weightLeft;
      _b = weightRight;
    }
    else {
      int64_t w = weightLeft + weightRight;
      int64_t isqrtW = fastIrsqrt(w);
      _a = fastIsqrt(weightLeft) * isqrtW >> 40 + kFISqrtFracBits - kFPFracBits;
      _b = fastIsqrt(weightRight) * isqrtW >> 40 + kFISqrtFracBits - kFPFracBits;
    }
  }

  void fwdTransform(int64_t& lf, int64_t& hf)
  {
    auto tmp = lf * _b;
    lf = fpReduce<kFPFracBits>(lf * _a + hf * _b);
    hf = fpReduce<kFPFracBits>(hf * _a - tmp);
  }

  void invTransform(int64_t& left, int64_t& right)
  {
    int64_t tmp = right * _b;
    right = fpReduce<kFPFracBits>(right * _a + left * _b);
    left = fpReduce<kFPFracBits>(left * _a - tmp);
  }

  int64_t getW0() { return _a; }
  int64_t getW1() { return _b; }

private:
  int64_t _a, _b;
};

//============================================================================
// Encapsulation of an Integer Haar transform stage.

class HaarKernel {
public:
  inline
  HaarKernel(int weightLeft, int weightRight, bool scaledWeights=false) { }

  void fwdTransform(int64_t& lf, int64_t& hf)
  {
    hf -= lf;
    lf += (hf >> 1 + kFPFracBits) << kFPFracBits;
  }

  void invTransform(int64_t& left, int64_t& right)
  {
    left -= (right >> 1 + kFPFracBits) << kFPFracBits;
    right += left;
  }

  int64_t getW0() { return 1; }
  int64_t getW1() { return 1; }

};

//============================================================================
// In-place transform a set of sparse 2x2x2 blocks each using the same weights

template<bool haarFlag>
struct FwdTransformBlock222 {
  template<class Kernel = RahtKernel>
  static inline void
  apply(
    const int numBufs, int64_t* buf, const int64_t weights[])
  {
    static const int a[4 + 4 + 4] = {0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3};
    static const int b[4 + 4 + 4] = {1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7};
    const int64_t* w = &weights[32];
    for (int i = 0; i < 12; i++) {
      int64_t w0 = *w++;
      int64_t w1 = *w++;

      if (w0 || w1) {
        int i0 = a[i];
        int i1 = b[i];

        if (!w0) {
          for (int k = 0; k < numBufs; k++)
            std::swap(buf[8 * k + i0], buf[8 * k + i1]);
        }
        else if (w1) { // only one occupied, propagate to next level
          // actual transform
          Kernel kernel(w0, w1, true);
          for (int k = 0; k < numBufs; k++) {
            kernel.fwdTransform(buf[8 * k + i0], buf[8 * k + i1]);
          }
        }
      }
    }
  }

  template<int numBufs, class Kernel = RahtKernel>
  static void
  apply(int64_t* buf, const int64_t weights[])
  {
    apply<Kernel>(numBufs, buf, weights);
  }
};

template<>
struct FwdTransformBlock222<true> {
  template<class Kernel = HaarKernel>
  static inline void
  apply(const int numBufs, int64_t* buf, const bool weights[])
  {
    static const int a[4 + 4 + 4] = { 0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3 };
    static const int b[4 + 4 + 4] = { 1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7 };
    const bool* w = &weights[32];
    for (int i = 0; i < 12; i++) {
      bool w0 = *w++;
      bool w1 = *w++;

      if (w0 || w1) {
        int i0 = a[i];
        int i1 = b[i];

        if (!w0) {
          for (int k = 0; k < numBufs; k++)
            std::swap(buf[8 * k + i0], buf[8 * k + i1]);
        }
        else if (w1) { // only one occupied, propagate to next level
          for (int k = 0; k < numBufs; k++) {
            buf[8 * k + i1] -= buf[8 * k + i0];
            buf[8 * k + i0] += (buf[8 * k + i1] >> 1 + kFPFracBits) << kFPFracBits;
          }
        }
      }
    }
  }

  template<int numBufs, class Kernel = HaarKernel>
  static void
  apply(int64_t* buf, const bool weights[])
  {
    apply<Kernel>(numBufs, buf, weights);
  }
};



template<bool haarFlag>
struct InvTransformBlock222 {
  template<int numBufs, bool computekernel = false, class Kernel = RahtKernel>
  static void
  apply(int64_t* buf, const int64_t weights[])
  {
    static const int a[4 + 4 + 4] = { 0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3 };
    static const int b[4 + 4 + 4] = { 1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7 };
    const int64_t* w = &weights[33 + 22];
    for (int i = 11; i >= 0; i--) {

      int64_t w1 = *(w--);
      int64_t w0 = *(w--);

      if (w0 || w1) {
        int i0 = a[i];
        int i1 = b[i];

        if (!w0) {
          for (int k = 0; k < numBufs; k++)
            buf[8 * k + i1] = buf[8 * k + i0];
        }
        else if (w1) { // only one occupied, propagate to next level
          // actual transform
          Kernel kernel(w0, w1, !computekernel);
          for (int k = 0; k < numBufs; k++) {
            kernel.invTransform(buf[8 * k + i0], buf[8 * k + i1]);
          }
        }
      }
    }
  }
};

template<>
struct InvTransformBlock222<true> {
  template<int numBufs, bool computekernel = false, class Kernel = HaarKernel>
  static void
  apply(int64_t* buf, const bool weights[])
  {
    static const int a[4 + 4 + 4] = { 0, 2, 4, 6, 0, 4, 1, 5, 0, 1, 2, 3 };
    static const int b[4 + 4 + 4] = { 1, 3, 5, 7, 2, 6, 3, 7, 4, 5, 6, 7 };
    const bool* w = &weights[33 + 22];
    for (int i = 11; i >= 0; i--) {

      bool w1 = *(w--);
      bool w0 = *(w--);

      if (w0 || w1) {
        int i0 = a[i];
        int i1 = b[i];

        if (!w0) {
          for (int k = 0; k < numBufs; k++)
            buf[8 * k + i1] = buf[8 * k + i0];
        }
        else if (w1) { // only one occupied, propagate to next level
          for (int k = 0; k < numBufs; k++) {
            buf[8 * k + i0] -= (buf[8 * k + i1] >> 1 + kFPFracBits) << kFPFracBits;
            buf[8 * k + i1] += buf[8 * k + i0];
          }
        }
      }
    }
  }
};

//============================================================================

template <typename NodeIterator>
struct OneLevelRAHTNodesTraversal {

  //typedef typename std::vector<RahtNode>::const_iterator NodeIterator;

  const int level;
  const int parentLevel;

  const NodeIterator itEnd;

  const int64_t maskXlevel = RSO_RAHT::maskX & ~RSO_RAHT::packX((1 << level) - 1);
  const int64_t maskXparentLevel = RSO_RAHT::maskX & ~RSO_RAHT::packX((1 << parentLevel) - 1);
  const int64_t maskYparentLevel = RSO_RAHT::maskY & ~RSO_RAHT::packY((1 << parentLevel) - 1);
  const int64_t maskZparentLevel = RSO_RAHT::maskZ & ~RSO_RAHT::packZ((1 << parentLevel) - 1);

  const int64_t maskXYparentLevel = maskXparentLevel | maskYparentLevel;
  const int64_t maskXYZparentLevel = maskXYparentLevel | maskZparentLevel;
  // const int64_t maskXYZlevel = maskXYZparentLevel | RSO_RAHT::pos(1 << level);

  const std::array<int64_t, 8> nodePosOffset = {
    0,
    RSO_RAHT::packZ(1 << level),
    RSO_RAHT::packY(1 << level),
    RSO_RAHT::packZ(1 << level) + RSO_RAHT::packY(1 << level),
    RSO_RAHT::packX(1 << level),
    RSO_RAHT::packX(1 << level) + RSO_RAHT::packZ(1 << level),
    RSO_RAHT::packX(1 << level) + RSO_RAHT::packY(1 << level),
    RSO_RAHT::packX(1 << level) + RSO_RAHT::packY(1 << level) + RSO_RAHT::packZ(1 << level),
  };

  std::array<NodeIterator, 8> nodeIt;
  std::array<int64_t, 8> nodeItPos;

  // finished tubes and slices, as bit flags
  // should be wrong now: tubes or slices are already finished or will finish after the current node
  int endedTubes;
  int endedSlices;

  int64_t nodePos;
  uint8_t nodeOccupancy;

  bool finished() const { return nodeIt[0] >= itEnd; }

  inline int64_t getNodeItPos(NodeIterator nodeIt) const
  { return nodeIt != itEnd ? nodeIt->pos : RSO_RAHT::maxPosXYZ; }

  inline void updateNodeItPos(int itIdx)
  { nodeItPos[itIdx] = getNodeItPos(nodeIt[itIdx]); }

  // return X from already set nodeItPos[0]
  inline int64_t computeNodePosX() const
  { return nodeItPos[0] & maskXparentLevel; }

  // return X,Y from already set nodeItPos[0] and nodeItPos[4]
  inline int64_t computeNodePosXY() const
  {
    // => hacked nodeItPos to be set to RSO_RAHT::maxPosXYZ
    return std::min(nodeItPos[0] & maskXYparentLevel, nodeItPos[4] & maskXYparentLevel);
  }

  // return X,Y,Y from already set nodeItPos[0], nodeItPos[2], nodeItPos[4] and nodeItPos[6]
  inline int64_t computeNodePosXYZ() const
  {
    // => hacked nodeItPos to be set to RSO_RAHT::maxPosXYZ
    return std::min(
      std::min(nodeItPos[0] & maskXYZparentLevel, nodeItPos[2] & maskXYZparentLevel),
      std::min(nodeItPos[4] & maskXYZparentLevel, nodeItPos[6] & maskXYZparentLevel));
  }

  //--------------------------------------------------------------------------

  OneLevelRAHTNodesTraversal(int nodeLevel, const NodeIterator& begin, const NodeIterator& end)
    : level(nodeLevel)
    , parentLevel(nodeLevel + 1)
    , itEnd(end)
  {
    assert(nodeLevel >= 0);

    // It is assumed:
    //  - 'rasterNodes' are ordered in raster scan (i.e lexicographic) order,
    //  - there are no duplicated nodes at level 'nodeLevel',
    //  - for each node of 'rasterNodes', (node.pos & ~maskXYZlevel) == 0
    //    (i.e. lower than nodeLevel bits are equal to zero)

    // each iterator will point to the node corresponding to current child node
    // position, or the next one if it does not exist

    // -- init --

    nodeIt[0] = begin;

    nodeItPos[0] = nodeIt[0] != itEnd ? nodeIt[0]->pos : RSO_RAHT::maxPosXYZ;

    int64_t nodePosX = computeNodePosX();

    searchSecondSliceIt(nodePosX);
    initEndedSlices(nodePosX);

    int64_t nodePosXY = computeNodePosXY();

    searchSecondTubeEachSliceIt(nodePosXY);
    initEndedTubes(nodePosXY);

    int64_t nodePosXYZ = computeNodePosXYZ();

    searchSecondNodeEachTubeIt(nodePosXYZ);

    nodePos = nodePosXYZ;
  }

  // determine occupancy and position of first node
  void determineOccupancyAndPosition() {

    int64_t nodePosXYZ = computeNodePosXYZ();

    nodeOccupancy = 0;
    // second node each tube
    for (int nodeIdx = 0; nodeIdx < 8; ++nodeIdx) {
      int64_t childNodePos = nodePosXYZ + nodePosOffset[nodeIdx];
      nodeOccupancy |= (nodeItPos[nodeIdx] == childNodePos) << nodeIdx;
    }

    nodePos = nodePosXYZ;
  }

  // occupancy and position of first node is provided (e.g. from parent node)
  void setOccupancyAndPosition(uint8_t occupancy, int64_t position) {
    nodePos = position;
    nodeOccupancy = occupancy;
  }

  // advance next knowing the occupancy
  void next() {

    for (int tubeIdx = 0; tubeIdx < 4; ++tubeIdx) {
      // update first node iterator of each tube
      auto occTube = 3 & nodeOccupancy >> 2 * tubeIdx;
      if (occTube) {
        nodeIt[2 * tubeIdx] = nodeIt[2 * tubeIdx + 1] + (occTube >= 2);
        updateNodeItPos(2 * tubeIdx);
        auto endedTube = (nodeItPos[2 * tubeIdx] & RSO_RAHT::maskXY) > nodePos + nodePosOffset[2 * tubeIdx];
        endedTubes += endedTube << tubeIdx;
        if (endedTube) {
          // hack: save position as second tube's node pos
          nodeItPos[2 * tubeIdx + 1] = nodeItPos[2 * tubeIdx];
          nodeItPos[2 * tubeIdx] = RSO_RAHT::maxPosXYZ;
        }
      }
    }

    // check if we have changed of parent tube
    bool tubeChange = endedTubes == 0xF;

    if (tubeChange) {
      // parent tube is changing for next one
      for (int sliceIdx = 0; sliceIdx < 2; ++sliceIdx) {
        if (!(endedSlices & 1 << sliceIdx)) {
          nodeIt[4 * sliceIdx] = nodeIt[4 * sliceIdx + 2];
          // hack retrieve position saved as second tube's node pos
          nodeItPos[4 * sliceIdx] = nodeItPos[4 * sliceIdx + 2 + 1];
          auto endedSlice = (nodeItPos[4 * sliceIdx] & RSO_RAHT::maskX) > nodePos + nodePosOffset[4 * sliceIdx];
          endedSlices += endedSlice << sliceIdx;
          if (endedSlice) {
            nodeItPos[4 * sliceIdx] = RSO_RAHT::maxPosXYZ;
          }
        }
      }
    }

    // check if we have changed of parent slice
    bool sliceChange = endedSlices == 0x3;

    if (sliceChange) {
      // parent slice is changing for next one
      nodeIt[0] = nodeIt[4];
      updateNodeItPos(0);

      int64_t nodePosX = computeNodePosX();

      searchSecondSliceIt(nodePosX);
      initEndedSlices(nodePosX);
    }

    if (tubeChange) {
      int64_t nodePosXY = computeNodePosXY();

      searchSecondTubeEachSliceIt(nodePosXY);
      initEndedTubes(nodePosXY);
    }

    int64_t nodePosXYZ = computeNodePosXYZ();

    searchSecondNodeEachTubeIt(nodePosXYZ);

    nodePos = nodePosXYZ;
  }

private:

  void searchSecondSliceIt(int64_t nodePosX)
  {
    // second slice
    int64_t secondSliceMinPos = nodePosX + nodePosOffset[4];
    nodeIt[4] = std::lower_bound(nodeIt[0], itEnd, secondSliceMinPos,
      [](decltype(*nodeIt[0]) v, decltype(secondSliceMinPos) end) {
        return v.pos < end;});

    updateNodeItPos(4);
  }

  void initEndedSlices(int64_t nodePosX)
  {
    endedSlices = 0;
    for (int sliceIdx = 0; sliceIdx < 2; ++sliceIdx) {
      auto endedSlice = (nodeItPos[4 * sliceIdx] & RSO_RAHT::maskX) > nodePosX + nodePosOffset[4 * sliceIdx];
      endedSlices += endedSlice << sliceIdx;
      if (endedSlice) {
        for (int sliceTubeIdx = 0; sliceTubeIdx < 2; ++sliceTubeIdx) {
          nodeItPos[4* sliceIdx + 2 * sliceTubeIdx] = RSO_RAHT::maxPosXYZ;
        }
      }
    }
  }

  void searchSecondTubeEachSliceIt(int64_t nodePosXY)
  {
    // second tube each slice
    int64_t secondTubeMinPos = nodePosXY + nodePosOffset[2];
    nodeIt[2] = nodeIt[0];
    while (nodeIt[2] != itEnd && nodeIt[2]->pos < secondTubeMinPos) ++nodeIt[2];

    updateNodeItPos(2);

    secondTubeMinPos = nodePosXY + nodePosOffset[6];
    nodeIt[6] = nodeIt[4];
    while (nodeIt[6] != itEnd && nodeIt[6]->pos < secondTubeMinPos) ++nodeIt[6];

    updateNodeItPos(6);
  }

  void initEndedTubes(int64_t nodePosXY)
  {
    endedTubes = 0;
    for (int tubeIdx = 0; tubeIdx < 4; ++tubeIdx) {
      auto endedTube = (nodeItPos[2 * tubeIdx] & RSO_RAHT::maskXY) > nodePosXY + nodePosOffset[2 * tubeIdx];
      endedTubes += endedTube << tubeIdx;
      if (endedTube) {
        // hack: save position as second tube's node pos
        nodeItPos[2 * tubeIdx + 1] = nodeItPos[2 * tubeIdx];
        nodeItPos[2 * tubeIdx] = RSO_RAHT::maxPosXYZ;
      }
    }
  }

  void searchSecondNodeTubeIt(int64_t nodePosXYZ, int tubeIdx)
  {
    int64_t secondNodeMinPos = nodePosXYZ + nodePosOffset[2 * tubeIdx + 1];
    nodeIt[2 * tubeIdx + 1] = nodeIt[2 * tubeIdx] + (nodeItPos[2 * tubeIdx] < secondNodeMinPos);
    updateNodeItPos(2 * tubeIdx + 1);
  }

  void searchSecondNodeEachTubeIt(int64_t nodePosXYZ)
  {
    for (int tubeIdx = 0; tubeIdx < 4; ++tubeIdx) {
      searchSecondNodeTubeIt(nodePosXYZ, tubeIdx);
    }
  }

};

//============================================================================

} /* namespace RAHT */

//============================================================================

}  // namespace pcc
