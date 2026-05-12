/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
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
#include "TMC3.h"
#include <vector>

#include "PCCPointSet.h"
#include "entropy.h"
#include "hls.h"

#include <queue>
#include <tuple>

#include <queue>
#include <tuple>

#include <queue>
#include <tuple>

#include "ringbuf.h"
#include<unordered_map>

namespace pcc {

//============================================================================
static const unsigned int motionParamPrec = 16;
static const unsigned int motionParamScale = 1 << motionParamPrec;
static const unsigned int motionParamOffset = 1 << (motionParamPrec - 1);
struct PCCOctree3Node;
struct MSOctree;
struct ParameterSetMotion;
struct EncodeMotionSearchParams;
struct EncodeSkipParams;
struct EncoderParams;
//============================================================================

struct MotionEntropy {
  AdaptiveBitModel splitPu;
  StaticBitModel mvSign;
  AdaptiveBitModel mvAllZero;
  AdaptiveBitModel mvIsZero;
  AdaptiveBitModel mvIsOne;
  AdaptiveBitModel mvIsTwo;
  AdaptiveBitModel mvIsThree;
  AdaptiveBitModel _ctxLocalMV;

  void reset();
};

//----------------------------------------------------------------------------
class MotionEntropyEncoder : protected MotionEntropy {
public:
  MotionEntropyEncoder(
    const MotionEntropy& motionEntropy, EntropyEncoder* arithmeticEncoder)
    : MotionEntropy(motionEntropy)
    , _arithmeticEncoder(arithmeticEncoder)
  {}

  MotionEntropyEncoder(EntropyEncoder* arithmeticEncoder)
    : _arithmeticEncoder(arithmeticEncoder)
  {}

  const MotionEntropy& getCtx() const { return *this; }

  MotionEntropy& getCtx() { return *this; }

  EntropyEncoder* getArithmeticEncoder() const { return _arithmeticEncoder; }

  // local
  void encodeSplitPu(int symbol);
  void encodeVector(const point_t& mv);

private:
  EntropyEncoder* _arithmeticEncoder;
};

//----------------------------------------------------------------------------
class MotionEntropyDecoder : protected MotionEntropy {
public:
  MotionEntropyDecoder(
    const MotionEntropy& motionEntropy, EntropyDecoder* arithmeticDecoder)
    : MotionEntropy(motionEntropy)
    , _arithmeticDecoder(arithmeticDecoder)
  {}

  MotionEntropyDecoder(EntropyDecoder* arithmeticDecoder)
    : _arithmeticDecoder(arithmeticDecoder)
  {}

  const MotionEntropy& getCtx() const { return *this; }
  EntropyDecoder* getArithmeticDecoder() const { return _arithmeticDecoder; }

  //local
  bool decodeSplitPu();
  void decodeVector(point_t* mv);

private:
  EntropyDecoder* _arithmeticDecoder;
};

//==========================================================================

//----------------------------------------- LOCAL MOTION -------------------
struct MVField {
  static constexpr uint64_t kMaskX = 0x0FFFFF0000000000ULL;
  static constexpr uint64_t kMaskY = 0x000000FFFFF00000ULL;
  static constexpr uint64_t kMaskZ = 0x00000000000FFFFFULL;
  static constexpr int kOffsetX = 40;
  static constexpr int kOffsetY = 20;
  static constexpr int kOffsetZ = 0;

  static constexpr uint32_t kNotSetMVIdx = 0xFFFFFF;

  struct PUNode {
    PUNode()
      : _mvIdx(kNotSetMVIdx)
      , _childsMask(0)
    {}
    // should be 16 bytes aligned
    // TODO: use masks getter/setters to ensure bit order
    struct {
      union {
        uint32_t _firstChildIdx:24; // if _childsMask != 0
        uint32_t _mvIdx:24; // if _childsMask == 0, kNotSetMVIdx means not set
      };
      uint32_t _childsMask:8; // 0 means leaf node or not set
    };
    struct {
      uint32_t _reserved:27;
      uint32_t _puSizeLog2:5;
    };
    uint64_t _packedPos0; // z | y << 20 | x << 40

    point_t pos0() const
    { return point_t {
        int32_t((_packedPos0 & kMaskX) >> kOffsetX),
        int32_t((_packedPos0 & kMaskY) >> kOffsetY),
        int32_t((_packedPos0 & kMaskZ) >> kOffsetZ)};
    }

    void set_pos0(point_t _pos)
    {
      _packedPos0 =
        (int64_t(_pos[0]) << kOffsetX)
        + (int64_t(_pos[1]) << kOffsetY)
        + (int64_t(_pos[2]) << kOffsetZ);
    }
  };
  //uint32_t puSizeLog2;
  uint32_t numRoots = 0;
  // (numRoots PUs in first layer)
  std::vector<PUNode> puNodes;
  // mvPool stores the motion vectors for child Nodes
  std::vector<point_t> mvPool;

  MVField() = default;
  MVField(const MVField&) = default;
  MVField(MVField&&) = default;
  MVField(const MVField& from, point_t begin, point_t end)
  {
    puNodes.reserve(from.puNodes.size());
    mvPool.reserve(from.mvPool.size());

    auto pos1 = end - 1;
    auto pos0 = begin;
    // For now taking entire node
    // TODO: generate new by only taking the intersection ? => not so easy
    for (int i = 0; i < from.numRoots; ++i) {
      const auto& nodeFrom = from.puNodes[i];
      auto nodePos0 = nodeFrom.pos0();
      auto nodePos1 = nodePos0 + ((1 << nodeFrom._puSizeLog2) - 1);
      int intersect = 0;
      for (int k=0; k < 3; ++k) {
        intersect |= std::min(nodePos1[k], pos1[k]) - std::max(nodePos0[k], pos0[k]);
      }
      if (intersect >= 0)
        puNodes.push_back(nodeFrom);
    }
    numRoots = puNodes.size();

    int N = puNodes.size();
    // TODO: depth first would be better
    for (int i = 0; i < N; ++i) {
      auto& node = puNodes[i];
      if (node._childsMask) {
        uint32_t childIdxFrom = node._firstChildIdx;
        node._firstChildIdx = N;
        for (int c = 0; c < 8; ++c) {
          if (node._childsMask & (1 << c)) {
            const auto& childNodeFrom =  from.puNodes[childIdxFrom++];
            puNodes.push_back(childNodeFrom);
            ++N;
          }
        }
      } else {
        uint32_t mvIdxFrom = node._mvIdx;
        node._mvIdx = mvPool.size();
        mvPool.push_back(from.mvPool[mvIdxFrom]);
      }
    }
  }
  MVField & operator =(const MVField&) = default;
  MVField & operator =(MVField&&) = default;

  void clear()
  {
    numRoots = 0;
    puNodes.clear();
    mvPool.clear();
  }

  point_t estimateMotion(point_t pos) const
  {
    if (!numRoots)
      return 0;

    // take motion from approximate nearest neighbor node to point pos

    // identify closest tree root
    int i_min;
    int d_min = std::numeric_limits<int>::max();
    for (int i = 0; i < numRoots; ++i) {
      const auto& node = puNodes[i];
      auto nodePos0 = node.pos0();
      auto nodePos1 = nodePos0 + (1 << node._puSizeLog2) - 1;
      bool overlap = true;
      for (int k=0; k < 3; ++k) {
        overlap &= pos[k] <= nodePos1[k] && pos[k] >= nodePos0[k];
      }
      if (overlap) {
        i_min = i;
        d_min = 0;
        break;
      } else {
        const auto dPos0 = nodePos0 - pos;
        const auto dPos1 = nodePos1 + dPos0;
        int local_d_min
          = (dPos0[0] > 0 ? dPos0[0] : 0)
          + (dPos0[1] > 0 ? dPos0[1] : 0)
          + (dPos0[2] > 0 ? dPos0[2] : 0)
          - (dPos1[0] < 0 ? dPos1[0] : 0)
          - (dPos1[1] < 0 ? dPos1[1] : 0)
          - (dPos1[2] < 0 ? dPos1[2] : 0);

        if (local_d_min < d_min) {
          d_min = local_d_min;
          i_min = i;
        }
      }
    }

    // identify closest leaf node
    while (true) {
      auto& node = puNodes[i_min];
      if (node._childsMask) {
        auto childIdx = node._firstChildIdx;
        d_min = std::numeric_limits<int>::max();
        for (int c = 0; c < 8; ++c) {
          if (node._childsMask & (1 << c)) {
            const auto& childNode =  puNodes[childIdx++];
            auto nodePos0 = childNode.pos0();
            auto nodePos1 = nodePos0 + (1 << childNode._puSizeLog2) - 1;

            const auto dPos0 = nodePos0 - pos;
            const auto dPos1 = nodePos1 + dPos0;
            int local_d_min
              = (dPos0[0] > 0 ? dPos0[0] : 0)
              + (dPos0[1] > 0 ? dPos0[1] : 0)
              + (dPos0[2] > 0 ? dPos0[2] : 0)
              - (dPos1[0] < 0 ? dPos1[0] : 0)
              - (dPos1[1] < 0 ? dPos1[1] : 0)
              - (dPos1[2] < 0 ? dPos1[2] : 0);

            if (local_d_min < d_min) {
              d_min = local_d_min;
              i_min = childIdx - 1;
            }
          }
        }
      } else {
        auto mvIdx = node._mvIdx;
        return mvPool[mvIdx];
      }
    }
  }
};

bool apply_splitPU_MC_color(
  const MSOctree& mSOctree,
  const AttributeParameterSet& aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  PCCPointSet3* compensatedPointCloud,
  int64_t& sumDGeom);

template <bool mcap>
int64_t encode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  const AttributeParameterSet* aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  MotionEntropyEncoder& motionEncoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2 = false,
  int S = -1,
  int S2 = -1);

// motion decoder

template <bool mcap>
int64_t decode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  const AttributeParameterSet* aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  MotionEntropyDecoder& motionDecoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2 = false,
  int S = -1,
  int S2 = -1);

//============================================================================

bool activateSkipMode(
  const PCCPointSet3& pointCloud,
  PCCOctree3Node& node0,
  PCCPointSet3& compensatedPointCloud,
  AdaptiveBitModel ctxCopyMode,
  std::unordered_map<uint32_t, std::array<int, 2>>& DepthNodeNum,
  const EncoderParams& param
  );

//============================================================================

void computeSkipNodesStatsEncoder(
  const MSOctree& octree,
  const int skipNodeSize,
  std::unordered_map<uint32_t, std::array<int, 2>>& DepthNodeNum
);

//============================================================================

class MotionEntropyEstimate;
struct MSOctree {
  MSOctree& operator=(const MSOctree&) = default;
  MSOctree& operator=(MSOctree&&) = default;
  MSOctree() = default;
  MSOctree(const MSOctree&) = default;
  MSOctree(MSOctree&&) = default;
  MSOctree(
    PCCPointSet3* predPointCloud,
    point_t offsetOrigin,
    uint32_t leafSizeLog2 = 0,
    uint32_t minDepth = 0
    );

  struct MSONode {
    uint32_t start;
    uint32_t end;
    std::array<uint32_t, 8> child = {}; // 0 means none
    int32_t sizeMinus1;
    point_t pos0;
    uint32_t parent;
    uint32_t reserved; // to align to 64 bytes (otherwise could be avoided)

    uint32_t numPoints() const { return end - start; }
  };

  point_t offsetOrigin; // offset applied to origin while construction the octree
  uint32_t maxDepth; // depth of full octree to get unitary sized nodes
  uint32_t depth; // depth of the motion search octree
  PCCPointSet3* pointCloud;
  std::vector<MSONode> nodes;

  void allocRingBuffers() {
    a = ringbuf<int>(nodes.size());
    b = ringbuf<int>(nodes.size());
  }

  enum MaxNumNN {
    kSingleNN = 1,
    kMaxNumNNForMCAP = 4,
  };

  // 12 bits precision fixed point division approximation
  static const int lutDivMCAP_fp12[kMaxNumNNForMCAP + 1];

  int
  nearestNeighbour_updateDMax(point_t pos, int32_t& d_max, bool approximate = false) const;

  inline
  int
  iNearestNeighbour_updateDMax(const point_t& pos, int32_t& d_max) const;

  inline
  int
  iApproximateNearestNeighbour_updateDMax(const point_t& pos, int32_t& d_max) const;

  template <MSOctree::MaxNumNN NN> inline
  int
  iApproxNearestNeighbourAttr(const point_t& pos, int nearestIdx[]) const;

  double
  find_motion(
    const EncodeMotionSearchParams& param,
    const ParameterSetMotion& mvPS,
    MotionEntropyEstimate& motionEntropy,
    const PCCPointSet3& Block0,
    const point_t& xyz0,
    int local_size,
    MVField& mvField,
    uint32_t puNodeIdx, // node Idx in mvField
    const point_t V00 = 0
  ) const;

  uint64_t
  find_dist(
    const int step,
    const PCCPointSet3& Block0
  ) const;

  void
  apply_motion(
    const point_t currNodePos0,
    const point_t currNodePos1,
    point_t Mvd,
    PCCOctree3Node* node0,
    PCCPointSet3* compensatedPointCloud,
    uint32_t depthMax = UINT32_MAX,
    bool flagNonPow2 = false,
    int S = -1,
    int S2 = -1
  ) const;

  template <MSOctree::MaxNumNN NN> inline
  int64_t
  apply_recolor_motion(
    point_t Mvd,
    PCCOctree3Node* node0,
    PCCPointSet3& pointCloud
  ) const;

  // find the smallest tree node containing the overall intersection between
  // the tree and a given cuboid, or noting if there is no intersection
  int
  nodeIdxIfIntersects(point_t pos0, uint32_t nodeSizeLog2) const {
    const point_t begin = pos0 - offsetOrigin;
    const point_t end = begin + (1 << nodeSizeLog2);
    if (  (end[0] > 0) & (begin[0] < 1 << maxDepth)
        & (end[1] > 0) & (begin[1] < 1 << maxDepth)
        & (end[2] > 0) & (begin[2] < 1 << maxDepth)
    ) {
      const point_t intersect_begin = {
        std::max(0, begin[0]),
        std::max(0, begin[1]),
        std::max(0, begin[2])
      };
      const point_t intersect_end = {
        std::min(1 << maxDepth, end[0]),
        std::min(1 << maxDepth, end[1]),
        std::min(1 << maxDepth, end[2])
      };
      int32_t nodeIdx = 0;
      const uint32_t targetNodeSizeLog2 = ilog2(uint32_t(std::max({
        intersect_end[0] - intersect_begin[0],
        intersect_end[1] - intersect_begin[1],
        intersect_end[2] - intersect_begin[2] }) - 1)) + 1;
        //  ^--- if true cuboid would occupy a wider node in the tree
      const int depthMax = std::min(depth, maxDepth - targetNodeSizeLog2);
      int pointChildMask = 1 << maxDepth - 1;
      for (int currDepth = 0; currDepth < depthMax; ++currDepth) {
        const int childIdx
          = (!!((intersect_begin[2]) & pointChildMask))
          | (!!((intersect_begin[1]) & pointChildMask) << 1)
          | (!!((intersect_begin[0]) & pointChildMask) << 2);

        const auto & node = nodes[nodeIdx];
        if (!node.child[childIdx])
          return -1;

        nodeIdx = node.child[childIdx];
        pointChildMask >>= 1;
      }
      return nodeIdx;
    }
    return -1;
  }

  mutable ringbuf<int> a; // for search
  mutable ringbuf<int> b; // for search
};

//----------------------------------------------------------------------------

bool
motionSearchForNode(
  const MSOctree& mSOctreeOrig,
  const MSOctree& mSOctree,
  const EncodeMotionSearchParams& msParams,
  const ParameterSetMotion& mvPS,
  int pointIdxNodeStart,
  int pointIdxNodeEnd,
  Vec3<int32_t> nodePos0,
  int nodeSize,
  MotionEntropyEncoder& motionEncoder,
  MVField& mvField,
  const MVField* refMvField,
  uint32_t puNodeIdx, // node Idx in mvField
  bool flagNonPow2 = false,
  int S = -1,
  int S2 = -1
);

//============================================================================

struct InterPredParams {
  PCCPointSet3 referencePointCloud;
  MSOctree mSOctreeRef;
  PCCPointSet3 compensatedPointCloud;
  // Motion
  MVField mvField;
  // TMP hack
  mutable std::vector<attr_t> attributes_mc;
};

//============================================================================

}  // namespace pcc
