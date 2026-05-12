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
#include "TMC3.h"

#include "motionWip.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <set>
#include <map>
#include <vector>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "geometry_octree.h"
#include "PCCTMC3Encoder.h"

namespace pcc {

//============================================================================

struct MotionEntropyEstimate {
  MotionEntropyEstimate(const MotionEntropy& contexts);

  double estimateVector(const point_t& mv);

  double estimateSplit(bool splitFlag)
  {
    const double inv2powkFPP = 1.525878906250000e-05;  // 2^-kFPP ; kFPP =16
    double R = localEncoder.getRateEstimate(
        [&] () { localMotionEncoder.encodeSplitPu(splitFlag); },
        localMotionEntropy.splitPu
      ) * inv2powkFPP;
    return R;
  }

  EntropyEncoder localEncoder;
  MotionEntropyEncoder localMotionEncoder;
  MotionEntropy& localMotionEntropy;
};

//============================================================================

void
MotionEntropy::reset()
{
  this->~MotionEntropy();
  new (this) MotionEntropy;
}

//----------------------------------------------------------------------------
MotionEntropyEstimate::MotionEntropyEstimate(const MotionEntropy& contexts)
: localEncoder()
, localMotionEncoder(contexts, &localEncoder)
, localMotionEntropy(localMotionEncoder.getCtx())
{
  localEncoder.setBuffer(4096, nullptr);
  localEncoder.start();
}

//----------------------------------------------------------------------------
inline void
MotionEntropyEncoder::encodeSplitPu(int symbol)
{
  _arithmeticEncoder->encode(symbol, splitPu);
}

//----------------------------------------------------------------------------

inline bool
MotionEntropyDecoder::decodeSplitPu()
{
  return _arithmeticDecoder->decode(splitPu);
}

//----------------------------------------------------------------------------

inline void
MotionEntropyEncoder::encodeVector(
  const point_t& mv)
{
  if (mv[0] == 0 && mv[1] == 0 && mv[2] == 0) {
    _arithmeticEncoder->encode(1, mvAllZero);
    return;
  }
  else
    _arithmeticEncoder->encode(0, mvAllZero);

  for (int comp = 0; comp < 3; comp++) {
    int v = mv[comp];
    if (v == 0 && !(comp==2 && mv[0] == 0 && mv[1] == 0)) {
      _arithmeticEncoder->encode(1, mvIsZero);
    }
    else {
      _arithmeticEncoder->encode(0, mvIsZero);

      _arithmeticEncoder->encode(v < 0, mvSign);
      if (v < 0)
        v = -v;
      v--;
      _arithmeticEncoder->encode(v == 0, mvIsOne);

      if (!v) {
        continue;
      }
      v--;

      _arithmeticEncoder->encode(v == 0, mvIsTwo);
      if (!v) {
        continue;
      }
      v--;

      _arithmeticEncoder->encode(v == 0, mvIsThree);
      if (!v) {
        continue;
      }
      v--;

      // expGolomb on |v|-1 with truncation
      _arithmeticEncoder->encodeExpGolomb(uint32_t(v), 1, _ctxLocalMV);
    }
  }
}

//----------------------------------------------------------------------------

inline void
MotionEntropyDecoder::decodeVector(point_t* mv)
{
  if (_arithmeticDecoder->decode(mvAllZero)) {
    (*mv)[0] = 0;
    (*mv)[1] = 0;
    (*mv)[2] = 0;
    return;
  }

  for (int comp = 0; comp < 3; comp++) {
    if (_arithmeticDecoder->decode(mvIsZero) && !(comp == 2 && mv[0] == 0 && mv[1] ==0)) {
      (*mv)[comp] = 0;
      continue;
    }
    bool sign = _arithmeticDecoder->decode(mvSign);

    if (_arithmeticDecoder->decode(mvIsOne)) {
      (*mv)[comp] = sign ? -1 : 1;
      continue;
    }
    if (_arithmeticDecoder->decode(mvIsTwo)) {
      (*mv)[comp] = sign ? -2 : 2;
      continue;
    }
    if (_arithmeticDecoder->decode(mvIsThree)) {
      (*mv)[comp] = sign ? -3 : 3;
      continue;
    }

    int v = 4 + _arithmeticDecoder->decodeExpGolomb(1, _ctxLocalMV);
    if (sign)
      v = -v;
    (*mv)[comp] = v;
  }
}

//----------------------------------------------------------------------------

double
MotionEntropyEstimate::estimateVector(
  const point_t& mv)
{
  const double inv2powkFPP = 1.525878906250000e-05;  // 2^-kFPP; kFPP =16
  double R = localEncoder.getRateEstimate(
      [&] () { localMotionEncoder.encodeVector(mv); },
      localMotionEntropy
    ) * inv2powkFPP;//
  return R;
}



//----------------------------------------- LOCAL MOTION -------------------

//============================================================================

bool
apply_splitPU_MC_color(
  const MSOctree& mSOctree,
  const AttributeParameterSet& aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  PCCPointSet3* compensatedPointCloud,
  int64_t& sumDGeom)
{
  auto& puNode = mvField.puNodes[puNodeIdx];
  // --------------  non-split / terminal case  ----------------
  if (!puNode._childsMask) {
    // use MV
    point_t MV = mvField.mvPool[puNode._mvIdx];

    // motion compensated attribute projection
    if (aps.lossless_flag)
      sumDGeom += mSOctree.apply_recolor_motion<MSOctree::kSingleNN>(MV, node0, *compensatedPointCloud);
    else
      sumDGeom += mSOctree.apply_recolor_motion<MSOctree::kMaxNumNNForMCAP>(MV, node0, *compensatedPointCloud);

    return true;
  }

  // --------------- split case ----------------------
  return false;
}

//----------------------------------------------------------------------------

template <bool mcap>
int64_t
encode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  const AttributeParameterSet* aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  MotionEntropyEncoder& motionEncoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2)
{
  int64_t sumDGeom = 0;
  auto& puNode = mvField.puNodes[puNodeIdx];
  // --------------  non-split / terminal case  ----------------
  if (nodeSize <= param.motion_min_pu_size || !puNode._childsMask) {
    if (nodeSize > param.motion_min_pu_size) {
      motionEncoder.encodeSplitPu(0);
    }

    // encode MV
    point_t MV = mvField.mvPool[puNode._mvIdx];
    motionEncoder.encodeVector(MV);

    if (!mcap) {
      mSOctree.apply_motion(
        node0->pos * nodeSize,
        (node0->pos + 1) * nodeSize - 1,
        MV, node0, compensatedPointCloud, mSOctree.depth, flagNonPow2, S, S2);
    } else {
      if (aps->lossless_flag)
        sumDGeom = mSOctree.apply_recolor_motion<MSOctree::kSingleNN>(MV, node0, *compensatedPointCloud);
      else
        sumDGeom = mSOctree.apply_recolor_motion<MSOctree::kMaxNumNNForMCAP>(MV, node0, *compensatedPointCloud);
    }
    node0->isCompensated = true;
    return sumDGeom;
  }

  // --------------- split case ----------------------
  motionEncoder.encodeSplitPu(1);
  return sumDGeom;
}

// instanciate for geometry
template int64_t encode_splitPU_MV_MC<false>(
  const MSOctree& mSOctree,
  const AttributeParameterSet* aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  MotionEntropyEncoder& motionEncoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2);

// instanciate for mcap
template int64_t encode_splitPU_MV_MC<true>(
  const MSOctree& mSOctree,
  const AttributeParameterSet* aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  MotionEntropyEncoder& motionEncoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2);

//----------------------------------------------------------------------------

template <bool mcap>
int64_t
decode_splitPU_MV_MC(
  const MSOctree& mSOctree,
  const AttributeParameterSet* aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  MotionEntropyDecoder& motionDecoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2)
{
  int64_t sumDGeom = 0;
  // Note:
  // geometry mvField is stored only for attributes and local attributes processing
  // when attributes do not use dual motion
  auto& puNode = mvField.puNodes[puNodeIdx];
  puNode.set_pos0(node0->pos * nodeSize);
  puNode._puSizeLog2 = ilog2(uint32_t(nodeSize - 1)) + 1; // TODO: check if we need the true size at some point or clean

  // decode split flag
  bool split = false;
  if (nodeSize > param.motion_min_pu_size)
    split = motionDecoder.decodeSplitPu();

  if (!split) {  // not split
                 // decode MV
    point_t MV = 0;
    motionDecoder.decodeVector(&MV);
    puNode._childsMask = 0;
    puNode._mvIdx = mvField.mvPool.size();
    mvField.mvPool.emplace_back(MV);

    if (!mcap) {
      mSOctree.apply_motion(
        node0->pos * nodeSize,
        (node0->pos + 1) * nodeSize - 1,
        MV, node0, compensatedPointCloud, mSOctree.depth, flagNonPow2, S, S2);
    } else {
      if (aps->lossless_flag)
        sumDGeom = mSOctree.apply_recolor_motion<MSOctree::kSingleNN>(MV, node0, *compensatedPointCloud);
      else
        sumDGeom = mSOctree.apply_recolor_motion<MSOctree::kMaxNumNNForMCAP>(MV, node0, *compensatedPointCloud);
    }

    node0->isCompensated = true;
    return sumDGeom;
  }

  // split; nothing to do
  return sumDGeom;
}

//============================================================================

struct MSOctreeStackElt {
  MSOctreeStackElt(int32_t nodeIdx, int32_t childIdx)
    : nodeIdx(nodeIdx), childIdx(childIdx) {}
  int32_t nodeIdx;
  int32_t childIdx;
};

//============================================================================

void computeSkipNodesStatsEncoder(
  const MSOctree& octree,
  const int skipNodeSize,
  std::unordered_map<uint32_t, std::array<int, 2>>& DepthNodeNum
) {
  // depth first traversal of the tree
  std::vector<MSOctreeStackElt> stack;
  stack.reserve(32);

  stack.push_back(MSOctreeStackElt(0,-1));

  decltype(DepthNodeNum.begin()) it;
  bool skipSizeReached = false;

  while (!stack.empty()) {
    auto& curr = stack.back();

    const int currDepth = stack.size() - 1;
    const int nodeSize = 1 << octree.maxDepth - currDepth;

    auto& node0 = octree.nodes[curr.nodeIdx];

    if (curr.childIdx == -1) {
      if (nodeSize == skipNodeSize) {
        skipSizeReached = true;
        it = DepthNodeNum.emplace(
          std::make_pair(curr.nodeIdx, std::array<int, 2> {1, 0})).first;
      } else if (skipSizeReached && nodeSize < skipNodeSize) {
        ++(it->second[0]);
        if (currDepth == octree.depth - 1)
          for (int i = 0; i < 8; ++i)
            it->second[1] += node0.child[i] != 0;
      }

      if (currDepth == octree.depth - 1) {
        stack.pop_back();
        continue;
      }

      curr.childIdx = 0;
    }

    while (curr.childIdx < 8 && !node0.child[curr.childIdx])
      curr.childIdx++;

    if (curr.childIdx < 8) {
      stack.push_back(MSOctreeStackElt(node0.child[curr.childIdx], -1));
      curr.childIdx++;
    } else {
      stack.pop_back();
    }
  }
}

//============================================================================

bool activateSkipMode(
  const PCCPointSet3& pointCloud,
  PCCOctree3Node& node0,
  PCCPointSet3& compensatedPointCloud,
  AdaptiveBitModel ctxSkipMode,
  std::unordered_map<uint32_t, std::array<int, 2>>& DepthNodeNum,
  const EncoderParams& param)
{
  // compute skip distortion
  PCCPointSet3 BlockC, BlockP;
  BlockC.appendPartition(pointCloud, node0.start, node0.end);
  BlockP.appendPartition(compensatedPointCloud, node0.predStart, node0.predEnd);
  MSOctree BlockC_mSOctree = MSOctree(&BlockC, 0, 2);
  MSOctree BlockP_mSOctree = MSOctree(&BlockP, 0, 2);

  int64_t D1 = BlockP_mSOctree.find_dist(param.skip.subsampleStep, BlockC);
  int64_t D2 = BlockC_mSOctree.find_dist(param.skip.subsampleStep, BlockP);
  double D = 2. * std::max(D1, D2);

  // no skip if distortion significantly bigger than Trisoup
  double QstepTriSoup = param.skip.QstepTriSoup;
  double TriSoupSize = param.gbh.trisoup_node_size;
  const double meanSecondDerivativeEstimate = 0.015 * param.skip.strength;
  double TriSoupPointDistEstimate =  0.25*QstepTriSoup + TriSoupSize * TriSoupSize * 7.8125e-03 * meanSecondDerivativeEstimate;
  double DTrisoupEstimate1 = (BlockC.getPointCount() + BlockP.getPointCount()) * TriSoupPointDistEstimate;

  if (D > DTrisoupEstimate1)
    return false;

  int64_t H[2];
  ctxSkipMode.getEntropy<16>(H);
  double R0 = DepthNodeNum[node0.mSOctreeNodeIdx][0] * param.skip.occRate + DepthNodeNum[node0.mSOctreeNodeIdx][1] * param.skip.triRate;
  bool isSkip = D < param.skip.lambda * (R0 + fpToDouble<16>(H[0] - H[1]));
  return isSkip;
}

//============================================================================

// instanciate for geometry
template int64_t decode_splitPU_MV_MC<false>(
  const MSOctree& mSOctree,
  const AttributeParameterSet* aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  MotionEntropyDecoder& motionDecoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2);

// instanciate for mcap
template int64_t decode_splitPU_MV_MC<true>(
  const MSOctree& mSOctree,
  const AttributeParameterSet* aps,
  PCCOctree3Node* node0,
  MVField& mvField,
  uint32_t puNodeIdx,
  const ParameterSetMotion& param,
  int nodeSize,
  MotionEntropyDecoder& motionDecoder,
  PCCPointSet3* compensatedPointCloud,
  bool flagNonPow2,
  int S,
  int S2);

//============================================================================

MSOctree::MSOctree(
    PCCPointSet3* predPointCloud,
    point_t offsetOrigin,
    uint32_t leafSizeLog2,
    uint32_t minDepth
  )
  : pointCloud(predPointCloud)
  , offsetOrigin(offsetOrigin)
{
  PCCPointSet3& pointCloud(*this->pointCloud);
  if (!pointCloud.size())
    return;

  nodes.reserve(pointCloud.size());

  auto bbox = pointCloud.computeBoundingBox();

  assert(bbox.min[0] >= 0 && bbox.min[1] >= 0 && bbox.min[2] >= 0);

  const auto maxPos = std::max(std::max(bbox.max[0],bbox.max[1]), bbox.max[2]);

  maxDepth = leafSizeLog2;
  depth = 0;
  while (1 << maxDepth <= maxPos) {
    ++depth;
    ++maxDepth;
  }
  int numExtraDepths = 0;
  while (depth < minDepth) {
    ++numExtraDepths;
    ++depth;
    ++maxDepth;
  }

  // push the first node(s)
  for (int parentDepth = -1; parentDepth < numExtraDepths; ++parentDepth) {
    nodes.emplace_back();
    MSONode& node00 = nodes.back();
    node00.start = uint32_t(0);
    node00.end = uint32_t(pointCloud.getPointCount());
    node00.pos0 = point_t{0} + offsetOrigin;
    node00.parent = std::max(0, parentDepth);
    node00.sizeMinus1 = int32_t(1 << maxDepth - (parentDepth + 1)) - 1;
    if (parentDepth >= 0)
      nodes[parentDepth].child[0] = parentDepth + 1;
  }

  // constructing depth first might be better for memory
  std::vector<MSOctreeStackElt> stack;
  stack.reserve(32);

  if (depth) {
    stack.push_back(MSOctreeStackElt(nodes.size() - 1, -1));
  }

  while (!stack.empty()) {
    auto& curr = stack.back();

    const int currDepth = stack.size() + numExtraDepths - 1;
    const int childSizeLog2 = maxDepth - currDepth - 1;
    const int pointSortMask = 1 << childSizeLog2;
    const int childSizeMinus1 = pointSortMask - 1;

    if (curr.childIdx == -1) {
      MSONode& node0 = nodes[curr.nodeIdx]; // adress may be invalidated after emplace_back....
      std::array<int32_t, 8> childCounts = {};

      countingSort(
        PCCPointSet3::iterator(&pointCloud, node0.start),
        PCCPointSet3::iterator(&pointCloud, node0.end),
        childCounts, [=](const PCCPointSet3::Proxy& proxy) {
          const auto & point = *proxy;
          return !!(int(point[2]) & pointSortMask)
            | (!!(int(point[1]) & pointSortMask) << 1)
            | (!!(int(point[0]) & pointSortMask) << 2);
        });

      int32_t childNodeIdx = nodes.size();
      for (int i = 0; i < 8; ++i) {
        if (childCounts[i]) {
          node0.child[i] = childNodeIdx++;
        }
      }
      // create child nodes
      uint32_t childStart = node0.start;
      auto pos0 = node0.pos0;
      for (int i = 0; i < 8; ++i) {
        if (childCounts[i]) {
          nodes.emplace_back(); // node0 may be invalidated after that
          MSONode& child0 = nodes.back();
          uint32_t childEnd = childStart + childCounts[i];
          child0.start = childStart;
          child0.end = childEnd;
          childStart = childEnd;
          child0.pos0 = pos0 + (Vec3<int32_t>{i >> 2, (i >> 1) & 1, i & 1} << childSizeLog2);
          child0.sizeMinus1 = childSizeMinus1;
          child0.parent = curr.nodeIdx;
        }
      }

      if (currDepth == depth-1) {
        stack.pop_back();
        continue;
      }
      curr.childIdx = 0;
    }

    MSONode& node0 = nodes[curr.nodeIdx];

    while (curr.childIdx < 8 && !node0.child[curr.childIdx])
      curr.childIdx++;

    if (curr.childIdx < 8) {
      stack.push_back(MSOctreeStackElt(node0.child[curr.childIdx], -1));
      curr.childIdx++;
    }
    else {
      stack.pop_back();
    }
  }

  // now everything is built, we can offset the points
  for (int i=0; i < pointCloud.size(); ++i) {
    pointCloud[i] += offsetOrigin;
  }

  allocRingBuffers();
}

//----------------------------------------------------------------------------

int
MSOctree::nearestNeighbour_updateDMax(point_t pos, int32_t& d_max, bool approximate) const {
  if (approximate)
    return iApproximateNearestNeighbour_updateDMax(pos, d_max);
  else
    return iNearestNeighbour_updateDMax(pos, d_max);
}

//----------------------------------------------------------------------------

struct NNStackElt {
  NNStackElt() = default;
  NNStackElt(const NNStackElt&) = default;
  NNStackElt(int32_t nodeIdx, int16_t firstChildIdx)
    : firstChildIdx(firstChildIdx), childIdx(0), nodeIdx(nodeIdx) {}
  int32_t nodeIdx;
  int16_t childIdx;
  int16_t firstChildIdx;
};

inline
int
MSOctree::iNearestNeighbour_updateDMax(const point_t& pos, int32_t& d_max) const {
  std::array<NNStackElt,32> stack;
  int stack_last = -1;

  const int depthMax = depth;
  int32_t nodeIdx = 0;
  const MSONode* node = &nodes[nodeIdx];
  int currDepth = 0;
  const point_t posOf = pos - offsetOrigin;
  if (  (posOf[0] >= 0) & (posOf[0] < 1 << maxDepth)
      & (posOf[1] >= 0) & (posOf[1] < 1 << maxDepth)
      & (posOf[2] >= 0) & (posOf[2] < 1 << maxDepth)
  ) {
    int pointChildMask = 1 << maxDepth - 1;
    for (; currDepth < depthMax; ++currDepth) {
      const int childIdx
        = (!!((posOf[2]) & pointChildMask))
        | (!!((posOf[1]) & pointChildMask) << 1)
        | (!!((posOf[0]) & pointChildMask) << 2);

      if (!node->child[childIdx])
        break;

      stack[++stack_last] = NNStackElt(nodeIdx, childIdx);
      nodeIdx = node->child[childIdx];
      node = &nodes[nodeIdx];
      pointChildMask >>= 1;
    }
  }
  for (; currDepth < depthMax; ++currDepth) {
    int d_min = INT32_MAX;
    int i_min;
    for(int i = 0; i < 8; ++i)
      if (node->child[i]) {
        const MSONode& child = nodes[node->child[i]];

        const auto dPos0 = child.pos0 - pos;
        const auto dPos1 = child.sizeMinus1 + dPos0;
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
    const int childNodeIdx = node->child[i_min];
    stack[++stack_last] = NNStackElt(nodeIdx, i_min);
    nodeIdx = childNodeIdx;
    node = &nodes[childNodeIdx];
  }

  int32_t local_d_max = INT32_MAX;

  int nearestIdx = node->start;

  for (int i = node->start; i < node->end; ++i) {
    auto dPoint = pos - (*pointCloud)[i];
    int32_t d
      = std::abs(dPoint[0])
      + std::abs(dPoint[1])
      + std::abs(dPoint[2]);
    if (d < local_d_max) {
      local_d_max = d;
      nearestIdx = i;
    }
  }

  if (local_d_max < d_max)
    d_max = local_d_max;

  if (!local_d_max)
    return nearestIdx;

  while (stack_last >= 0) {
    auto& sn = stack[stack_last];

    const MSONode& node = nodes[sn.nodeIdx];

    while (sn.childIdx < 8
        && (!node.child[sn.childIdx] || sn.childIdx == sn.firstChildIdx))
      ++sn.childIdx;

    if (sn.childIdx < 8) {
      const auto childNodeIdx = node.child[sn.childIdx];
      const MSONode& child = nodes[childNodeIdx];

      const auto dPos0 = child.pos0 - pos;
      const auto dPos1 = child.sizeMinus1 + dPos0;
      const auto d_min
        = (dPos0[0] > 0 ? dPos0[0] : 0)
        + (dPos0[1] > 0 ? dPos0[1] : 0)
        + (dPos0[2] > 0 ? dPos0[2] : 0)
        - (dPos1[0] < 0 ? dPos1[0] : 0)
        - (dPos1[1] < 0 ? dPos1[1] : 0)
        - (dPos1[2] < 0 ? dPos1[2] : 0);

      if (d_min < local_d_max) {
        if (stack_last < depthMax-1) {
          stack[++stack_last] = NNStackElt(childNodeIdx, -1);
        }
        else {
          for (int i = child.start; i < child.end; ++i) {
            auto dPoint = pos - (*pointCloud)[i];
            int32_t d
              = std::abs(dPoint[0])
              + std::abs(dPoint[1])
              + std::abs(dPoint[2]);
            if (d < local_d_max) {
              local_d_max = d;
              nearestIdx = i;
            }
          }
        }
      }
      ++sn.childIdx;
    } else {
      --stack_last;
    }
  }

  if (local_d_max < d_max)
    d_max = local_d_max;

  return nearestIdx;
}

//----------------------------------------------------------------------------

inline
int
MSOctree::iApproximateNearestNeighbour_updateDMax(const point_t& pos, int32_t& d_max) const {
  const int depthMax = depth;
  int32_t nodeIdx = 0;
  const MSONode* node = &nodes[nodeIdx];
  int currDepth = 0;
  const point_t posOf = pos - offsetOrigin;
  if (  (posOf[0] >= 0) & (posOf[0] < 1 << maxDepth)
      & (posOf[1] >= 0) & (posOf[1] < 1 << maxDepth)
      & (posOf[2] >= 0) & (posOf[2] < 1 << maxDepth)
  ) {
    int pointChildMask = 1 << maxDepth - 1;
    for (; currDepth < depthMax; ++currDepth) {
      const int childIdx
        = (!!((posOf[2]) & pointChildMask))
        | (!!((posOf[1]) & pointChildMask) << 1)
        | (!!((posOf[0]) & pointChildMask) << 2);

      if (!node->child[childIdx])
        break;

      nodeIdx = node->child[childIdx];
      node = &nodes[nodeIdx];
      pointChildMask >>= 1;
    }
  }
  for (; currDepth < depthMax; ++currDepth) {
    int d_min = INT32_MAX;
    int i_min;
    for(int i = 0; i < 8; ++i)
      if (node->child[i]) {
        const MSONode& child = nodes[node->child[i]];

        const auto dPos0 = child.pos0 - pos;
        const auto dPos1 = child.sizeMinus1 + dPos0;
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
    const int childNodeIdx = node->child[i_min];
    nodeIdx = childNodeIdx;
    node = &nodes[childNodeIdx];
  }

  int32_t local_d_max = INT32_MAX;

  int nearestIdx = node->start;

  for (int i = node->start; i < node->end; ++i) {
    auto dPoint = pos - (*pointCloud)[i];
    int32_t d
      = std::abs(dPoint[0])
      + std::abs(dPoint[1])
      + std::abs(dPoint[2]);
    if (d < local_d_max) {
      local_d_max = d;
      nearestIdx = i;
    }
  }

  if (local_d_max < d_max)
    d_max = local_d_max;

  return nearestIdx;
}

//----------------------------------------------------------------------------

template <MSOctree::MaxNumNN NN> inline
int
MSOctree::iApproxNearestNeighbourAttr(const point_t& pos, int nearestIdx[]) const {
  const int depthMax = depth;
  int32_t nodeIdx = 0;
  const MSONode* node = &nodes[nodeIdx];
  int currDepth = 0;
  const point_t posOf = pos - offsetOrigin;
  if (  (posOf[0] >= 0) & (posOf[0] < 1 << maxDepth)
      & (posOf[1] >= 0) & (posOf[1] < 1 << maxDepth)
      & (posOf[2] >= 0) & (posOf[2] < 1 << maxDepth)
  ) {
    int pointChildMask = 1 << maxDepth - 1;
    for (; currDepth < depthMax; ++currDepth) {
      const int childIdx
        = (!!((posOf[2]) & pointChildMask))
        | (!!((posOf[1]) & pointChildMask) << 1)
        | (!!((posOf[0]) & pointChildMask) << 2);

      if (!node->child[childIdx])
        break;

      nodeIdx = node->child[childIdx];
      node = &nodes[nodeIdx];
      pointChildMask >>= 1;
    }
  }
  for (; currDepth < depthMax; ++currDepth) {
    int d_min = INT32_MAX;
    int i_min;
    for(int i = 0; i < 8; ++i)
      if (node->child[i]) {
        const MSONode& child = nodes[node->child[i]];

        const auto dPos0 = child.pos0 - pos;
        const auto dPos1 = child.sizeMinus1 + dPos0;
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
    const int childNodeIdx = node->child[i_min];
    nodeIdx = childNodeIdx;
    node = &nodes[childNodeIdx];
  }

  int32_t d_max = INT32_MAX;

  int numN = 0;
  for (int i = node->start; i < node->end; ++i) {
    auto dPoint = pos - (*pointCloud)[i];
    int32_t d
      = std::abs(dPoint[0])
      + std::abs(dPoint[1])
      + std::abs(dPoint[2]);
    if (d < d_max) {
      d_max = d;
      numN = 1;
      nearestIdx[0] = i;
    } else if (NN > 1 && d == d_max) {
      numN -= numN == NN;
      nearestIdx[numN] = i;
      ++numN;
    }
  }

  return numN;
}

//----------------------------------------------------------------------------

uint64_t
MSOctree::find_dist(
  const int step,
  const PCCPointSet3& Block0) const
{
  uint64_t Dist = 0;
  for (int Nb = 0, idx = 0; Nb < Block0.size(); Nb += step, ++idx) {
    int32_t min_d = INT32_MAX;
    auto p = Block0[Nb];
    nearestNeighbour_updateDMax(p, min_d, false);
    Dist += min_d;
  }
  return Dist * step;
}

//----------------------------------------------------------------------------

double
MSOctree::find_motion(
  const EncodeMotionSearchParams& param,
  const ParameterSetMotion& mvPS,
  MotionEntropyEstimate& motionEntropy,
  const PCCPointSet3& Block0,
  const point_t& xyz0,
  int local_size,
  MVField& mvField,
  uint32_t puNodeIdx,
  const point_t V00) const // node Idx in mvField
{
  auto RDO = motionEntropy.localEncoder.makeRDO(
    param.lambda, motionEntropy.localMotionEntropy);

  RDO.start();
  // ---------------------------- test no split --------------------
  RDO.startAlternative();

  //const int max_distance = 3 * wSize +10;//3 * wSize;
  point_t bestV_NoSplit = V00;
  int jumpBlock = 1 + (Block0.size() >> param.decimate);  // (kind of) random sampling of the original block to code

  // TODO: buffer vector difference or dmax to estimate dmax after motion is applied
  const bool hasColor = Block0.hasColors();
  const int NtestedPoints = (Block0.size() + jumpBlock - 1) / jumpBlock;
  std::vector<int32_t> min_d0(NtestedPoints);
  std::vector<int32_t> min_dK(NtestedPoints);
  std::vector<int32_t> min_dTmp(NtestedPoints);
  std::vector<int32_t> min_start;

  std::vector<std::pair<pcc::point_t,pcc::Vec3<int>>> Buffer(NtestedPoints);

  // initial vector
  point_t startMV = V00;
  point_t V0 = startMV;
  int Dist = 0;
  int dColor_forMinD = 0;
  for (int Nb = 0, idx = 0; Nb < Block0.size(); Nb += jumpBlock, ++idx) {
    int32_t min_d = INT32_MAX;
    auto p = Block0[Nb];
    Buffer[idx].first = p;
    int nearestPointIdx = nearestNeighbour_updateDMax(p + V0, min_d, param.approximate_nn);

    if (hasColor) {
      auto cW = pointCloud->getColor(nearestPointIdx);
      auto cB = Block0.getColor(Nb);
      Buffer[idx].second[0] = cB[0];
      Buffer[idx].second[1] = cB[1];
      Buffer[idx].second[2] = cB[2];

      dColor_forMinD += std::abs(cW[0] - cB[0]);
      dColor_forMinD += std::abs(cW[1] - cB[1]);
      dColor_forMinD += std::abs(cW[2] - cB[2]);
    }

    Dist += min_d;
    min_d0[idx] = min_d;
    startMV += (*pointCloud)[nearestPointIdx] - p;
  }

  auto Vtemp = startMV.abs();
  int maxV = std::max(std::max(Vtemp[0], Vtemp[1]), Vtemp[2]);
  if (maxV)
    for (int k = 0; k < 3; k++)
      startMV[k] = 2 * Vtemp[k] < maxV ? 0 : (startMV[k] > 0 ? 1 : -1);

  const double scaleToBlock = double(Block0.size()) / NtestedPoints;

  double d = (param.d_geom_factor * Dist + param.d_color_factor * dColor_forMinD) * scaleToBlock;
  double cost = d + param.lambda * motionEntropy.estimateVector(V00);
  bestV_NoSplit = V0;

  // set loop search parameters
  double best_d = d;
  double best_cost = cost;
  std::set<int64_t> list_tested = { ((V00[0] + 32768LL) << 32) + ((V00[1] + 32768LL) << 16) + V00[2] + 32768LL };

  point_t VPrev = V00;
  V0 = startMV;

  while(true) {
    int64_t V1D = ((V0[0] + 32768LL) << 32) + ((V0[1] + 32768LL) << 16) + V0[2] + 32768LL;
    if (!list_tested.emplace(V1D).second)
      break;

    point_t meanV = 0;

    Dist = 0;
    dColor_forMinD = 0;
    auto offset = V0 - VPrev;
    int absOffset = std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]);
    for (int idx = 0; idx < NtestedPoints; ++idx) {
      auto p = Buffer[idx].first;
      int32_t min_d = min_d0[idx] + absOffset;
      int nearestPointIdx = nearestNeighbour_updateDMax(p + V0, min_d, param.approximate_nn);

      if (hasColor) {
        auto cW = pointCloud->getColor(nearestPointIdx);
        auto cB = Buffer[idx].second;
        dColor_forMinD += std::abs(cW[0] - cB[0]);
        dColor_forMinD += std::abs(cW[1] - cB[1]);
        dColor_forMinD += std::abs(cW[2] - cB[2]);
      }

      Dist += min_d;
      min_dTmp[idx] = min_d;
      meanV += (*pointCloud)[nearestPointIdx] - p;
    } // loop on points of block

    auto Vtemp = meanV.abs();
    int maxV = std::max(std::max(Vtemp[0], Vtemp[1]), Vtemp[2]);
    if (maxV)
      for (int kk = 0; kk < 3; kk++)
        meanV[kk] = 2 * Vtemp[kk] < maxV ? 0 : (meanV[kk] > 0 ? 1 : -1);

    d = (param.d_geom_factor * Dist + param.d_color_factor * dColor_forMinD) * scaleToBlock;
    cost = d + param.lambda * motionEntropy.estimateVector(V0);

    if (cost < best_cost) {
      best_cost = cost;
      best_d = d;
      bestV_NoSplit = V0;
      std::swap(min_d0,min_dTmp);

      // start next step adding meanV
      VPrev = V0;
      V0 += meanV;
    }
    else
      break;
  }

  const int searchPattern[3 * 18] = { 1,0,0, -1,0,0,  0,0,1, 0,0,-1,  0,1,0, 0,-1,0,  1,1,0, -1,-1,0,  0,1,1, 0,-1,-1,  1,0,1, -1,0,-1,  1,-1,0,  -1,1,0,  0,-1,1,  0,1,-1,  1,0,-1, -1,0,1, };
  //const int searchPattern[3 * 6] = {   1,0,0,   0,1,0,  0,-1,0,    -1,0,0,    0,0,1,   0,0,-1};
  //const int searchPattern[3 * 26] = { 1,0,0,   0,1,0,  0,-1,0,   -1,0,0,   0,0,1,   0,0,-1,  1,1,0,  1,-1,0,  -1,1,0,   -1,-1,0,  1,0,1,  0,1,1,  0,-1,1,  -1,0,1,  1,0,-1,  0,1,-1,  0,-1,-1,  -1,0,-1,  1,1,1,  1,1,-1,  1,-1,1,  1,-1,-1,  -1,1,1,  -1,1,-1,  -1,-1,1,  -1,-1,-1 };

  // loop MV search
  int Amotion = param.Amotion0;
  int lastT = 0;
  while (Amotion >= 1) {

    // loop on searchPattern
    const int* pSearch = searchPattern -3;
    bool flagBetter = false;
    for (int t = 0; t < 19; t++, pSearch += 3) {

      point_t V;
      if (t)
        V = point_t(pSearch[0] * Amotion, pSearch[1] * Amotion, pSearch[2] * Amotion);
      else
        V = point_t(searchPattern[3 * lastT] * Amotion, searchPattern[3 * lastT + 1] * Amotion, searchPattern[3 * lastT + 2] * Amotion);
      V0 = bestV_NoSplit + V;

      int64_t V1D = ((V0[0] + 32768LL) << 32) + ((V0[1] + 32768LL) << 16) + V0[2] + 32768LL;
      if(!list_tested.emplace(V1D).second)
        continue;

      Dist = 0;
      dColor_forMinD = 0;
      auto offset = V;
      int absOffset = std::abs(offset[0]) + std::abs(offset[1]) + std::abs(offset[2]);
      for (int idx = 0; idx < NtestedPoints; ++idx) {
        auto p = Buffer[idx].first;
        int32_t min_d = min_d0[idx] + absOffset;
        int nearestPointIdx = nearestNeighbour_updateDMax(p + V0, min_d, param.approximate_nn);

        if (hasColor) {
          auto cW = pointCloud->getColor(nearestPointIdx);
          auto cB = Buffer[idx].second;
          dColor_forMinD += std::abs(cW[0] - cB[0]);
          dColor_forMinD += std::abs(cW[1] - cB[1]);
          dColor_forMinD += std::abs(cW[2] - cB[2]);
        }

        min_dK[idx] = min_d;
        Dist += min_d;
      } // loop on points of block

      d = (param.d_geom_factor * Dist + param.d_color_factor * dColor_forMinD) * scaleToBlock;
      cost = d + param.lambda * motionEntropy.estimateVector(V0);

      // keep 2 best MV
      if (cost < best_cost) {
        best_cost = cost;
        best_d = d;
        bestV_NoSplit = V0;
        flagBetter = true;
        std::swap(min_d0,min_dK);

        if (t)
          lastT = t - 1;
        break;
      }
    }  // end loop on searchPattern

    // log reduction of search range
    if (!flagBetter) {
      Amotion >>= 1;
      if (!bestV_NoSplit[0] && !bestV_NoSplit[1] && !bestV_NoSplit[2]) {
        bool flag = Amotion > 1;
        Amotion >>= 1;
        if (flag)
          Amotion = std::max(Amotion, 1);
      }
    }

  }  // end loop MV search

  // cost split flag
  if (local_size > mvPS.motion_min_pu_size) {
    // cost no split flag
    motionEntropy.localMotionEncoder.encodeSplitPu(false);
  }
  motionEntropy.localMotionEncoder.encodeVector(bestV_NoSplit);

  RDO.finishAlternative(best_d).first;

  auto node = &mvField.puNodes[puNodeIdx];

  // ---------------------------- test split --------------------
  auto numPUsBeforeSplit = mvField.puNodes.size();
  auto numMVsBeforeSplit = mvField.mvPool.size();

  double d_Split = 0;
  bool split_taken = false;

  if (local_size > mvPS.motion_min_pu_size && Block0.size() >= 8) {
    RDO.startAlternative();

    // condition on number of points for search acceleration
    int local_size1 = local_size >> 1;

    std::array<point_t, 8> list_xyz = {
      point_t(0, 0, 0) + xyz0,
      point_t(0, 0, local_size1) + xyz0,
      point_t(0, local_size1, 0) + xyz0,
      point_t(0, local_size1, local_size1) + xyz0,
      point_t(local_size1, 0, 0) + xyz0,
      point_t(local_size1, 0, local_size1) + xyz0,
      point_t(local_size1, local_size1, 0) + xyz0,
      point_t(local_size1, local_size1, local_size1) + xyz0
    };

    // loop on 8 child PU

    PCCPointSet3 Block1;
    Block1.reserve(Block0.size());

    // NOTE: points shall be already ordered
    //       and all belonging to the current node
    std::array<int32_t, 8> childCounts = {};
    // child idx
    int childIdx = 0;
    // child PU coordinates
    point_t xyz1 = list_xyz[childIdx];
    // block for child PU
    point_t xyz1High = xyz1 + local_size1;
    // TODO: might be simplified / accelerated
    // by dichotomic search or suited comparizon according to index
    // => childIdx outside and index in block0
    for (const auto& b : Block0) {
      while (childIdx < 7 && (
        b[2] < xyz1[2] || b[2] >= xyz1High[2]
        || b[1] < xyz1[1] || b[1] >= xyz1High[1]
        || b[0] < xyz1[0] || b[0] >= xyz1High[0]
      )) {
        ++childIdx;
        xyz1 = list_xyz[childIdx];
        xyz1High = xyz1 + local_size1;
      }
      ++childCounts[childIdx];
    }

    node->_firstChildIdx = numPUsBeforeSplit;
    node->_childsMask = 0;
    for (int childIdx = 0; childIdx < 8; childIdx++) {
      if(childCounts[childIdx]) {
        node->_childsMask += 1 << childIdx;
        // add node for current pu
        mvField.puNodes.emplace_back();
        auto& childNode = mvField.puNodes.back();
        // address may have changed with emplace_back()
        node = &mvField.puNodes[puNodeIdx];
        xyz1 = list_xyz[childIdx];
        childNode.set_pos0(xyz1);
        childNode._puSizeLog2 = node->_puSizeLog2 - 1;
      }
    }
    assert(node->_childsMask);

    // cost split flag
    motionEntropy.localMotionEncoder.encodeSplitPu(true);
    int childStart = 0;
    uint32_t childNodeIdx = node->_firstChildIdx;
    for (int childIdx = 0; childIdx < 8; childIdx++) {
      // TODO: check the true cost and do not apply for attributes
      //cost_Split += 1.0 * param.lambda; // the cost due to not coding the occupancy with inter pred (only for geo !!)
      if(!childCounts[childIdx]) {  // empty PU
        continue;
      }
      Block1.resize(0);
      Block1.appendPartition(
        Block0, childStart, childStart + childCounts[childIdx]);
      childStart += childCounts[childIdx];

      xyz1 = list_xyz[childIdx];

      d_Split += find_motion(
        param, mvPS, motionEntropy, Block1, xyz1, local_size1, mvField,
        childNodeIdx, bestV_NoSplit);

      ++childNodeIdx;
    }
    auto rdo_res = RDO.finishAlternative(d_Split);
    split_taken = rdo_res.second;
  }
  RDO.finish();

  // ---------------------------- choose split vs no split --------------------
  if (!split_taken) {  // no split
    mvField.puNodes.resize(numPUsBeforeSplit);
    mvField.mvPool.resize(numMVsBeforeSplit + 1);
    // address may have changed
    node = &mvField.puNodes[puNodeIdx];
    //node._firstChildIdx = 0;
    node->_mvIdx = numMVsBeforeSplit;
    node->_childsMask = 0; // not split
    mvField.mvPool.back() = bestV_NoSplit;
    return best_d;
  }
  else {
    return d_Split;
  }
}

//----------------------------------------------------------------------------

bool
motionSearchForNode(
  const MSOctree& mSOctreeOrig,
  const MSOctree& mSOctree,
  const EncodeMotionSearchParams& param,
  const ParameterSetMotion& mvPS,
  int pointIdxNodeStart,
  int pointIdxNodeEnd,
  Vec3<int32_t> nodePos0,
  int nodeSize,
  MotionEntropyEncoder& motionEncoder,
  MVField& mvField,
  const MVField* refMvField,
  uint32_t puNodeIdx, // node Idx in mvField
  bool flagNonPow2,
  int S,
  int S2)
{
  PCCPointSet3 Block0;
  Block0.appendPartition(*mSOctreeOrig.pointCloud, pointIdxNodeStart, pointIdxNodeEnd);

  // entropy estimates
  MotionEntropyEstimate mcEstimate(motionEncoder.getCtx());

  // motion search

  // scale/undice point position according to trisoup node size
  if (flagNonPow2) {
    int maskS = (1 << S2) - 1;
    for (int i = 0; i < Block0.size(); ++i) {
      Block0[i][0] = ((Block0[i][0] >> S2) * S) + (Block0[i][0] & maskS);
      Block0[i][1] = ((Block0[i][1] >> S2) * S) + (Block0[i][1] & maskS);
      Block0[i][2] = ((Block0[i][2] >> S2) * S) + (Block0[i][2] & maskS);
    }
  }

  point_t mv00 = 0;

  if (refMvField) {
    mv00 = refMvField->estimateMotion(nodePos0 + (nodeSize >> 1));
  }

  auto& node = mvField.puNodes[puNodeIdx];
  node.set_pos0(nodePos0);
  node._puSizeLog2 = ilog2(uint32_t(nodeSize - 1)) + 1; // TODO: check if we need the true size at some point or clean

  // MV search
  mSOctree.find_motion(
    param, mvPS, mcEstimate, Block0, nodePos0, nodeSize, mvField, puNodeIdx, mv00);

  return true;
}

void
MSOctree::apply_motion(
  const point_t currNodePos0,
  const point_t currNodePos1,
  const point_t MVd,
  PCCOctree3Node* node0,
  PCCPointSet3* compensatedPointCloud,
  uint32_t depthMax,
  bool flagNonPow2,
  int S,
  int S2) const
{
  auto &fifo = a;
  assert(fifo.empty());
  fifo.clear();
  depthMax = std::min(depthMax, depth);
  const auto node0Pos0 = currNodePos0 + MVd;
  const auto node0Pos1 = currNodePos1 + MVd;
  const auto minNodeSizeMinus1 = (1 << maxDepth - depthMax) - 1;

  auto &local = b;
  assert(local.empty());
  local.clear();
  int addedPointCount = 0;

  fifo.push(0);
  while (!fifo.empty()) {
    const MSONode& node = nodes[fifo.front()];
    //
    const auto nodeSizeMinus1 = node.sizeMinus1;
    if (minNodeSizeMinus1 == nodeSizeMinus1)
      break;

    const auto nodePos0 = node.pos0;
    const auto nodePos1 = nodePos0 + nodeSizeMinus1;
    const auto dPos0 = nodePos0 - node0Pos0;
    const auto dPos1 = node0Pos1 - nodePos1;

    if ( (dPos0[0] | dPos0[1] | dPos0[2]
        | dPos1[0] | dPos1[1] | dPos1[2]) >= 0
    ) {
      local.push(fifo.front());
      addedPointCount += node.end - node.start;
    } else {
      int intersect = 0;
      for (int k=0; k < 3; ++k) {
        intersect |= std::min(nodePos1[k], node0Pos1[k]) - std::max(nodePos0[k], node0Pos0[k]);
      }
      if (intersect >= 0)
        for(int i = 0; i < 8; ++i)
          if (node.child[i])
            fifo.push(node.child[i]);
    }
    fifo.pop();
  }

  std::vector<int> indices;
  indices.reserve(pointCloud->size());
  indices.resize(addedPointCount);
  int i = 0;
  while (!local.empty()) {
    const MSONode& node = nodes[local.front()];
    for (int k=node.start; k<node.end; ++k)
      indices[i++] = k;
    local.pop();
  }

  while (!fifo.empty()) {
    const MSONode& node = nodes[fifo.front()];
    for (int i=node.start; i < node.end; ++i) {
      auto const & pt = (*pointCloud)[i];
      const auto dPos0 = pt - node0Pos0;
      const auto dPos1 = node0Pos1 - pt;
      if ( (dPos0[0] | dPos0[1] | dPos0[2]
          | dPos1[0] | dPos1[1] | dPos1[2]) >= 0)
        indices.push_back(i);
    }
    fifo.pop();
  }

  node0->predStart = compensatedPointCloud->size();
  compensatedPointCloud->appendPartition(*pointCloud, indices, false);
  node0->predEnd = compensatedPointCloud->size();

  for (int i = node0->predStart; i < node0->predEnd; ++i) {
    auto& predPoint = (*compensatedPointCloud)[i];
    predPoint[0] -= MVd[0];
    predPoint[1] -= MVd[1];
    predPoint[2] -= MVd[2];
  }
  // align points with octree nodes
  static const int LUTdicingBits[33] = {0, 1, 1, 23, 2, 22, 22, 23, 3, 27, 23, 23, 23, 27, 24, 27, 4,
    28, 22, 25, 24, 28, 24, 26, 24, 26, 28, 25, 25, 26, 28, 29, 5};
  static const int64_t LUTdicingFactor[33] = {0, 2, 1, 2796203, 1, 838861, 699051, 1198373, 1, 14913081, 838861, 762601,
    699051, 10324441, 1198373, 8947849, 1, 15790321, 233017, 1766023, 838861, 12641, 762601, 2917777, 699051,
    2684355, 10324441, 1242757, 1198373, 2314099, 8947849, 17318417, 1 };

  if (flagNonPow2) {
    // TODO:
    //  should we apply at octree node level to avoid division and multiplications ?
    const int factorS = (1 << S2) - S;
    for (int i = node0->predStart; i < node0->predEnd; ++i) {
      auto& predPoint = (*compensatedPointCloud)[i];
      int temp0 = int64_t(predPoint[0]) * LUTdicingFactor[S] >> LUTdicingBits[S];
      int temp1 = int64_t(predPoint[1]) * LUTdicingFactor[S] >> LUTdicingBits[S];
      int temp2 = int64_t(predPoint[2]) * LUTdicingFactor[S] >> LUTdicingBits[S];
      predPoint[0] += temp0 * factorS;
      predPoint[1] += temp1 * factorS;
      predPoint[2] += temp2 * factorS;
    }
  }
}

const int MSOctree::lutDivMCAP_fp12[MSOctree::kMaxNumNNForMCAP + 1] =
  { 0, 4096, 2048, 1365, 1024 };


template <MSOctree::MaxNumNN NN> inline
int64_t
MSOctree::apply_recolor_motion(
  point_t Mvd,
  PCCOctree3Node* node0,
  PCCPointSet3& pointCloud) const
{
  int64_t sumDGeom = 0;
  for (int i = node0->start; i < node0->end; ++i) {
    int nearestPointIdx[NN];
    const auto p = pointCloud[i] + Mvd;
    int nN = iApproxNearestNeighbourAttr<NN>(p, nearestPointIdx);

    // Note: considering max precision for attribute is 16, 32 bits is enough
    // if bigger attributes were used, 64 bits should be used instead
    static_assert(std::numeric_limits<attr_t>::max() < 2 << 16);
    Vec3<int32_t> A = this->pointCloud->getColor(nearestPointIdx[0]);
    if (NN > 1 && nN > 1) {
      for (int j = 1; j < nN; ++j) {
        A += this->pointCloud->getColor(nearestPointIdx[j]);
      }
      A += nN >> 1;
      A *= lutDivMCAP_fp12[nN];
      A >>= 12;
    }
    pointCloud.setColor(i, A);

    auto dPoint = p - (*this->pointCloud)[nearestPointIdx[0]];
    const auto d =
      std::abs(dPoint[0])
      + std::abs(dPoint[1])
      + std::abs(dPoint[2]);
    sumDGeom += d;
  }
  return sumDGeom;
}

//----------------------------------------------------------------------------
}  // namespace pcc
