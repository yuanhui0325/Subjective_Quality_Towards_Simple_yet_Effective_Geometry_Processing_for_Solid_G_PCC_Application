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

#include "geometry.h"

#include "OctreeNeighMap.h"
#include "geometry_octree.h"
#include "geometry_trisoup_encoder.h"
#include "io_hls.h"
#include "tables.h"
#include "TMC3.h"
#include "PCCTMC3Encoder.h"
#include "motionWip.h"

namespace pcc {

//============================================================================

class GeometryOctreeEncoder {
public:
  GeometryOctreeContexts& ctx;

  GeometryOctreeEncoder(
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    GeometryOctreeContexts& ctxtMem,
    EntropyEncoder* arithmeticEncoder);

  GeometryOctreeEncoder(const GeometryOctreeEncoder&) = default;
  GeometryOctreeEncoder(GeometryOctreeEncoder&&) = default;

  // dynamic OBUF
  void clearMap() { ctx.clearMap(); };
  void resetMap(bool forTrisoup) { ctx.resetMap(forTrisoup); }

  void encodePositionLeafNumPoints(int count);

  void encodeOccupancyFullNeihbourgs(
    const RasterScanContext::occupancy& occ,
    int occupancy,
    int predOccupancy[8],
    bool isInter,
    bool trisoup);

  // local QU
  void createQU(localQU& qu, Vec3<int32_t> pos, int size, int baseQP);
  void encodeQU(localQU& qu, int baseQP);

  const GeometryOctreeContexts& getCtx() const { return ctx; }

public:
  EntropyEncoder* _arithmeticEncoder;
};

//============================================================================

GeometryOctreeEncoder::GeometryOctreeEncoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  GeometryOctreeContexts& ctxtMem,
  EntropyEncoder* arithmeticEncoder)
  : ctx(ctxtMem)
  , _arithmeticEncoder(arithmeticEncoder)
{
}

//============================================================================
// Encode the number of points in a leaf node of the octree.

void
GeometryOctreeEncoder::encodePositionLeafNumPoints(int count)
{
  int dupPointCnt = count - 1;
  _arithmeticEncoder->encode(dupPointCnt > 0, ctx._ctxDupPointCntGt0);
  if (dupPointCnt <= 0)
    return;

  _arithmeticEncoder->encodeExpGolomb(dupPointCnt - 1, 0, ctx._ctxDupPointCntEgl);
  return;
}

//-------------------------------------------------------------------------
// encode node occupancy bits
//

static void (*pointer2FunctionContext[8])(
  OctreeNeighours&, int, int&, int&, bool&) = {
  makeGeometryAdvancedNeighPattern0,
  makeGeometryAdvancedNeighPattern1,
  makeGeometryAdvancedNeighPattern2,
  makeGeometryAdvancedNeighPattern3,
  makeGeometryAdvancedNeighPattern4,
  makeGeometryAdvancedNeighPattern5,
  makeGeometryAdvancedNeighPattern6,
  makeGeometryAdvancedNeighPattern7};

//-------------------------------------------------------------------------
// encode node occupancy bits
//
void
GeometryOctreeEncoder::encodeOccupancyFullNeihbourgs(
  const RasterScanContext::occupancy& occ,
  int occupancy,
  int predOcc[8],
  bool isInter,
  bool trisoup)
{
  int coded0 = 0;
  OctreeNeighours octreeNeighours;
  prepareGeometryAdvancedNeighPattern(occ, octreeNeighours);

  // loop on occupancy bits from occupancy map
  for (int i = 0; i < 8; i++) {
    if (coded0 >= 7)  // bit is 1
      continue;

    // OBUF contexts
    int ctx1, ctx2;
    bool Sparse;
    (*pointer2FunctionContext[i])(octreeNeighours, occupancy, ctx1, ctx2, Sparse);

    int ctxTable = 0;
    if (!isInter) // INTRA
      ctxTable = Sparse;
    else { //  INTER
      ctxTable = 2 + predOcc[i];
    }

    // encode
    int bit = (occupancy >> i) & 1;
    if (Sparse) {
      if (trisoup)
        ctx2 >>= 4;

      auto obufIdx = ctx._MapOccupancySparse[isInter][i].getEvolve(
        bit, ctx2, ctx1, &ctx._OBUFleafNumber, ctx._BufferOBUFleaves);
      _arithmeticEncoder->encode(
        bit, obufIdx >> 3, ctx._CtxMapDynamicOBUF[ctxTable][obufIdx],
        ctx._CtxMapDynamicOBUF[ctxTable].obufSingleBound);

    }
    else {
      if (trisoup)
        ctx2 >>= 2;

      auto obufIdx = ctx._MapOccupancy[isInter][i].getEvolve(
        bit, ctx2, ctx1, &ctx._OBUFleafNumber, ctx._BufferOBUFleaves);
      _arithmeticEncoder->encode(
        bit, obufIdx >> 3, ctx._CtxMapDynamicOBUF[ctxTable][obufIdx],
        ctx._CtxMapDynamicOBUF[ctxTable].obufSingleBound);
    }

    // update partial occupancy of current node
    coded0 += !bit;
  }
}


//-------------------------------------------------------------------------
void
GeometryOctreeEncoder::createQU(
  localQU& qu, Vec3<int32_t> pos, int size, int baseQP)
{
  // simple encoder decision (for testing purpose):
  // any node with starting position having a 'y' lower than 800
  // is encoded with a lower quality
  // TODO: add encoding parameters for more versatile quality definition
  if (pos[1]  < 800) {   // [1] is vertical pos
    qu.isBaseParameters = false;
    qu.localQP = baseQP + 6;
  }
  else {
    qu.isBaseParameters = true;
    qu.localQP = baseQP;
  }
}

//-------------------------------------------------------------------------
void
GeometryOctreeEncoder::encodeQU(localQU& qu, int baseQP)
{
  _arithmeticEncoder->encode(qu.isBaseParameters,  ctx._ctxQUflag);
  if (qu.isBaseParameters)
    return;

  int diffQP = qu.localQP - baseQP;
  // sign
  _arithmeticEncoder->encode(diffQP > 0, ctx._ctxQUSign);

  // magnitude
  _arithmeticEncoder->encodeExpGolomb(
    std::abs(diffQP) - 1, 1, ctx._ctxQUQPpref, ctx._ctxQUQPsuf);
}

//-------------------------------------------------------------------------

template <bool forTrisoup>
void
encodeGeometryOctree(
  const EncoderParams& encParams,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& rootEncoder,
  TrisoupEncoder* trisoup
  )
{
  PCCPointSet3& compensatedPointCloud = interPredParams.compensatedPointCloud;

  uint64_t fp16IOSizeMotionBits = 0;

  PCCPointSet3& refPointCloud = interPredParams.referencePointCloud;
  refPointCloud = refFrame.cloud;
  PCCPointSet3& predPointCloud = refPointCloud;
  MSOctree& mSOctree = interPredParams.mSOctreeRef;
  auto& mvField = interPredParams.mvField;
  std::unique_ptr<MotionEntropyEncoder> motionEncoder;
  const bool isInter = gbh.slice_inter_prediction_flag;

  int log2MotionBlockSize = 0;
  if (isInter) {
    int log2MinPUSize = ilog2(uint32_t(gps.motion.motion_min_pu_size - 1)) + 1;
    mSOctree = MSOctree(&predPointCloud, -gbh.geomBoxOrigin, std::min(2, log2MinPUSize));
    motionEncoder.reset(new MotionEntropyEncoder(ctxtMemMotion, arithmeticEncoders[0].get()));

    // local motion prediction structure -> LPUs from predPointCloud
    log2MotionBlockSize = ilog2(uint32_t(gps.motion.motion_block_size - 1)) + 1;
    // LPU is bigger than root node, must adjust
    if (gbh.rootNodeSizeLog2 < log2MotionBlockSize)
      log2MotionBlockSize = gbh.rootNodeSizeLog2;
  }

  auto arithmeticEncoderIt = arithmeticEncoders.begin();
  GeometryOctreeEncoder encoder(gps, gbh, ctxtMem, arithmeticEncoderIt->get());

  // local QU
  if (forTrisoup) {
    if (gbh.qu_size_log2 > 0) {
      int quSize = (1 << gbh.qu_size_log2) * gbh.trisoupNodeSize(gps);
      std::cout << "QU size = " << quSize << "\n";
    }
    ctxtMem.quLastIndex = 0;
    ctxtMem.listOfQUs.clear();
  }

  // saved state for use with parallel bistream coding.
  // the saved state is restored at the start of each parallel octree level
  std::unique_ptr<GeometryOctreeContexts> savedState;

  const int S2 = gbh.trisoupNodeSizeLog2(gps);
  const int S = gbh.trisoupNodeSize(gps);
  const int factorS = (1 << S2) - S;
  const bool flagNonPow2 = factorS != 0;

  // process point cloud positions depending on TriSoup node size
  if (flagNonPow2) {
    for (int np = 0; np < pointCloud.getPointCount(); np++) {
      const int temp0 = pointCloud[np][0] / S;
      const int temp1 = pointCloud[np][1] / S;
      const int temp2 = pointCloud[np][2] / S;
      pointCloud[np][0] += temp0 * factorS;
      pointCloud[np][1] += temp1 * factorS;
      pointCloud[np][2] += temp2 * factorS;
    }
  }

  // init main fifo and reserve for max size
  int reservedBufferSize = pointCloud.getPointCount();
  if (gps.trisoup_enabled_flag && gbh.trisoupNodeSize(gps) > 1)
    reservedBufferSize = std::max(1000, reservedBufferSize >> 2 * gbh.trisoupNodeSizeLog2(gps) - 1);
  std::vector<PCCOctree3Node> fifo;
  std::vector<PCCOctree3Node> fifoNext;
  fifo.reserve(reservedBufferSize);
  fifoNext.reserve(reservedBufferSize);

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(pointCloud.getPointCount());
  node00.pos = int32_t(0);
  node00.mSOctreeNodeIdx = uint32_t(0);
  node00.predEnd = isInter ? predPointCloud.getPointCount() : uint32_t(0);
  node00.mSONodeIdx = isInter ? 0 : -1;
  node00.predStart = uint32_t(0);
  node00.hasMotion = 0;
  node00.isCompensated = 0;
  node00.isSkiped = false;
  node00.quIndex = -1;

  RasterScanContext rsc(fifo);

  // the termination depth of the octree phase
  int minNodeSizeLog2 = gbh.trisoupNodeSizeLog2(gps);
  int maxDepth = gbh.rootNodeSizeLog2 - minNodeSizeLog2;
  gbh.tree_depth_minus1 = maxDepth - 1;

  // Derived parameter used by local attributes
  if (!forTrisoup) {
    gbh.maxRootNodeSize = 1 << maxDepth;
    gbh.uniqueSlabBlock = gbh.uniqueSlabBlock && gbh.maxRootNodeSize <= std::min(
      sps.localized_attributes_slab_thickness_minus1 + 1, sps.localized_attributes_slab_block_size_minus1 + 1);
  }

  if (gps.octree_point_count_list_present_flag)
    gbh.footer.octree_lvl_num_points_minus1.reserve(maxDepth);

  if (!(isInter && gps.gof_geom_entropy_continuation_enabled_flag) && !gbh.entropy_continuation_flag) {
    encoder.clearMap();
    encoder.resetMap(forTrisoup);
  }

  // localized attributes point indexes
  std::vector<std::vector<int> > laPointIdx;
  point_t localSlabBlockIdxStart = 0;
  point_t numSlabBlocksPerDim = -1;
  point_t localSlabBlockStart = 0;
  point_t slabBlockSize;
  bool useLocalAttr = !forTrisoup && !gbh.uniqueSlabBlock;
  PCCPointSet3 localPointCloud;
  std::vector<int32_t> numPointsSlabBlocks;
  if (useLocalAttr) {
    slabBlockSize = {
      sps.localized_attributes_slab_thickness_minus1 + 1,
      sps.localized_attributes_slab_block_size_minus1 + 1,
      sps.localized_attributes_slab_block_size_minus1 + 1
    };
    numSlabBlocksPerDim = {
      (gbh.maxRootNodeSize + slabBlockSize[0] - 1) / slabBlockSize[0],
      (gbh.maxRootNodeSize + slabBlockSize[1] - 1) / slabBlockSize[1],
      (gbh.maxRootNodeSize + slabBlockSize[2] - 1) / slabBlockSize[2]
    };
    // in case we would have IDCM or thing like that we need to buffer all
    // potential slab blocks
    int numSlabBlocksPerSlab = numSlabBlocksPerDim[1] * numSlabBlocksPerDim[2];
    int numSlabBlocks = numSlabBlocksPerDim[0] * numSlabBlocksPerSlab;
    laPointIdx.resize(numSlabBlocks);
    localPointCloud.reserve(pointCloud.getPointCount());
    numPointsSlabBlocks.reserve(numSlabBlocksPerSlab);
  }

  PCCPointSet3 recPointCloud;
  recPointCloud.addRemoveAttributes(
    pointCloud.hasColors(), pointCloud.hasReflectances());
  recPointCloud.resize(pointCloud.getPointCount());
  int nRecPoints = 0;
  int nOutputPoints = 0;

  MSOctree mSOctreeCurr;
  if (isInter)
    mSOctreeCurr = MSOctree(
      &pointCloud, {}, gbh.trisoupNodeSizeLog2(gps),
      std::max(0, gbh.rootNodeSizeLog2 - gbh.trisoupNodeSizeLog2(gps)));

  // Note: the below hashtable is not really welcome, even if encoder only
  // TODO: try to avoid it
  std::unordered_map<uint32_t, std::array<int, 2>> DepthNodeNum;
  if (forTrisoup && isInter && gps.trisoup_early_skip_coding_mode_enabled_flag) {
    computeSkipNodesStatsEncoder(mSOctreeCurr, gps.trisoup_early_skip_coding_mode_node_size, DepthNodeNum);
  }

  // ----------  loop on depth ----------------
  int lastPos0 = INT_MIN; // used to detect slice change and call for TriSoup
  int nodeSizeLog2 = gbh.rootNodeSizeLog2;
  for (int depth = 0; depth < maxDepth; depth++, nodeSizeLog2--) {
    // set at the start of each level
    auto fifoCurrLvlEnd = fifo.end();
    int numNodesNextLvl = 0;
    const int childSizeLog2 = nodeSizeLog2 - 1;
    const int nodeSizePow2 = 1 << nodeSizeLog2;

    if (forTrisoup && isInter && gps.trisoup_early_skip_coding_mode_enabled_flag
        && nodeSizePow2 == gps.trisoup_early_skip_coding_mode_node_size) {
      trisoup->enableSkip = true;
    }

    // save context state for parallel coding
    if (depth == maxDepth - 1 - gbh.geom_stream_cnt_minus1)
      if (gbh.geom_stream_cnt_minus1)
        savedState.reset(new GeometryOctreeContexts(encoder.ctx));

    // load context state for parallel coding starting one level later
    if (depth > maxDepth - 1 - gbh.geom_stream_cnt_minus1) {
      encoder.ctx = *savedState;
      encoder._arithmeticEncoder = (++arithmeticEncoderIt)->get();
    }

    int currPUIdx = 0;
    if (isInter) {
      if (nodeSizeLog2 == log2MotionBlockSize) {
        // TODO: allocate motion field according to number of nodes
        uint32_t numNodes = fifo.end() - fifo.begin();
        mvField.puNodes.reserve(numNodes * 8); // arbitrary value.
        mvField.mvPool.reserve(numNodes * 8); // arbitrary value.
        // allocate all the PU roots
        mvField.numRoots = numNodes;
        mvField.puNodes.resize(numNodes);
      }
    }

    rsc.initializeNextDepth();
    bool isLastDepth = depth == maxDepth - 1;
    if (forTrisoup && isLastDepth)
      trisoup->leaves.reserve(8 * fifo.size());

    // process all nodes within a single level
    IterOneLevelSubnodesRSO(
      fifo.begin(),
      fifoCurrLvlEnd,
      [&](const decltype(fifo.begin())& itNode) -> const point_t& {
        return itNode->pos;
      },
      [&](const decltype(fifo.begin())& fifoCurrNode, int childIdx) {

      if (childIdx & 1 || isLeafNode(nodeSizeLog2))
        return;

      const int tubeIndex = (childIdx >> 1) & 1;
      const int nodeSliceIndex = (childIdx >> 2) & 1;
      bool firstVisit = !tubeIndex && !nodeSliceIndex;

      PCCOctree3Node& node0 = *fifoCurrNode;

      if (firstVisit) {
        //local motion : determine PU tree by motion search and RDO
        const int nodeSize = S * (1 << nodeSizeLog2 - S2);
        if (isInter && nodeSizeLog2 == log2MotionBlockSize) {
          node0.hasMotion = motionSearchForNode(
            mSOctreeCurr, mSOctree, encParams.motion, encParams.gps.motion,
            node0.start, node0.end, node0.pos * nodeSize, nodeSize,
            *motionEncoder, mvField, nullptr, currPUIdx, flagNonPow2, S, S2);
          node0.mvFieldNodeIdx = currPUIdx++;
        }

        // code split PU flag. If not split, code  MV and apply MC
        // results of MC are stacked in compensatedPointCloud that starts empty
        if (node0.mvFieldNodeIdx != -1) {
          uint64_t before = motionEncoder->getArithmeticEncoder()->getNumBitsEstimate();

          encode_splitPU_MV_MC<false>(mSOctree, nullptr,
            &node0, mvField, node0.mvFieldNodeIdx, gps.motion, nodeSize,
            *motionEncoder, &compensatedPointCloud,
            flagNonPow2, S, S2);

          uint64_t after = motionEncoder->getArithmeticEncoder()->getNumBitsEstimate();

          fp16IOSizeMotionBits += after - before;
        }

        if (forTrisoup) {
          // ---------- skip Mode ---------------
          if (gps.trisoup_early_skip_coding_mode_enabled_flag && node0.isCompensated
              && nodeSizePow2 == gps.trisoup_early_skip_coding_mode_node_size
              && node0.predStart < node0.predEnd) {
            // test early termination Skip mode pred inter here
            node0.isSkiped = activateSkipMode(
              pointCloud, node0, compensatedPointCloud,
              ctxtMem.ctxSkipMode, DepthNodeNum, encParams);

            encoder._arithmeticEncoder->encode(node0.isSkiped ? 1 : 0, ctxtMem.ctxSkipMode);
          }
        }

        // split the current node into 8 children
        //  - perform an 8-way counting sort of the current node's points
        //  - (later) map to child nodes
        if (isInter) {
          if (node0.mSOctreeNodeIdx || depth == 0) {
            auto& msoNode = mSOctreeCurr.nodes[node0.mSOctreeNodeIdx];
            // 8-way counting sort already made by mSOctreeCurr
            for (int childIndex = 0; childIndex < 8; ++childIndex) {
              auto msoChildIdx = msoNode.child[childIndex];
              if (msoChildIdx) {
                auto& msoChildNode = mSOctreeCurr.nodes[msoChildIdx];
                node0.childCounts[childIndex] = msoChildNode.end - msoChildNode.start;
              }
            }
          }
        } else {
          countingSort(
            PCCPointSet3::iterator(&pointCloud, node0.start),
            PCCPointSet3::iterator(&pointCloud, node0.end), node0.childCounts,
            [=](const PCCPointSet3::Proxy& proxy) {
              const auto& point = *proxy;
              return ((point[2] >> childSizeLog2) & 1)
                | (((point[1] >> childSizeLog2) & 1) << 1)
                | (((point[0] >> childSizeLog2) & 1) << 2);
            });
        }

        /// sort and partition the predictor for local motion
        // TODO: check we can remove
        if (isInter) {
          if (node0.isCompensated) {
            countingSort(
              PCCPointSet3::iterator(&compensatedPointCloud, node0.predStart),
              PCCPointSet3::iterator(&compensatedPointCloud, node0.predEnd),
              node0.predCounts, [=](const PCCPointSet3::Proxy& proxy) {
                const auto& point = *proxy;
                return ((point[2] >> childSizeLog2) & 1)
                  | (((point[1] >> childSizeLog2) & 1) << 1)
                  | (((point[0] >> childSizeLog2) & 1) << 2);
              });
          }
          else {
            if (!forTrisoup && depth < mSOctree.depth && node0.mSONodeIdx >= 0) {
              const auto& msoNode = mSOctree.nodes[node0.mSONodeIdx];
              for (int i = 0; i < 8; ++i) {
                uint32_t msoChildIdx = msoNode.child[i];
                if (msoChildIdx) {
                  const auto& msoChild = mSOctree.nodes[msoChildIdx];
                  node0.predCounts[i] = msoChild.end - msoChild.start;
                }
              }
            }
          }
        }

        node0.predPointsStartIdx = node0.predStart;

        // local QU
        if (forTrisoup) {
          if (gbh.qu_size_log2 && nodeSizeLog2 == S2 + gbh.qu_size_log2) {
            ctxtMem.listOfQUs.emplace_back();
            localQU& qu = ctxtMem.listOfQUs.back();
            encoder.createQU(
              qu, (node0.pos << nodeSizeLog2 - S2) * S,
              (1 << nodeSizeLog2 - S2) * S, gbh.trisoup_QP);
            encoder.encodeQU(qu, gbh.trisoup_QP);
            node0.quIndex = ctxtMem.quLastIndex++;
          }
        }

        // occupancy of child nodes
        if (!forTrisoup || !node0.isSkiped) {
          // generate the bitmap of child occupancy and count
          // the number of occupied children in node0.
          int occupancy = 0;
          for (int i = 0; i < 8; i++) {
            occupancy |= (node0.childCounts[i] > 0) << i;
          }
          node0.childOccupancy = occupancy;

          // occupancy inter predictor
          int predOccupancy[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
          for (int i = 0; i < 8; i++) {
            predOccupancy[i] += node0.predCounts[i] > std::max(0, nodeSizePow2 >> 8);
            predOccupancy[i] += node0.predCounts[i] > std::max(2, nodeSizePow2 >> 4);
            predOccupancy[i] += node0.predCounts[i] >
              (forTrisoup ? std::max(16, nodeSizePow2) : std::max(8, nodeSizePow2 >> 2));
          }

          // update contexts
          RasterScanContext::occupancy contextualOccupancy;
          GeometryNeighPattern gnp{};
          rsc.nextNode(&*fifoCurrNode, contextualOccupancy);
          node0.neighPattern = gnp.neighPattern = contextualOccupancy.neighPattern;
          gnp.adjNeighOcc[0] = contextualOccupancy.childOccupancyContext[4];
          gnp.adjNeighOcc[1] = contextualOccupancy.childOccupancyContext[10];
          gnp.adjNeighOcc[2] = contextualOccupancy.childOccupancyContext[12];

          // encode child occupancy map
          encoder.encodeOccupancyFullNeihbourgs(
            contextualOccupancy, node0.childOccupancy, predOccupancy, isInter, forTrisoup);
        }
        else {
          int predOcc = 0;
          for (int i = 0; i < 8; i++)
            predOcc |= (node0.predCounts[i] > 0) << i;
          node0.childOccupancy = predOcc;
        }

      } // end firstVisit

      // Leaf nodes are immediately coded.  No further splitting occurs.
      if (!forTrisoup && tubeIndex && nodeSliceIndex && isLeafNode(childSizeLog2)) {
        if (!useLocalAttr) {
          for (auto idx = node0.start; idx < node0.end; idx++, nOutputPoints++) {
            recPointCloud[nOutputPoints] = pointCloud[idx];
            if (pointCloud.hasColors())
              recPointCloud.setColor(nOutputPoints, pointCloud.getColor(idx));
            if (pointCloud.hasReflectances())
              recPointCloud.setReflectance(nOutputPoints, pointCloud.getReflectance(idx));
          }
        } else {
          int slabBlockIdx;
          int nodeposX = node0.pos[0] << 1 + childSizeLog2;
          int nodeposY = node0.pos[1] << 1 + childSizeLog2;
          int nodeposZ = node0.pos[2] << 1 + childSizeLog2;
          if (isLastDepth) {
            // selection of the right Slab Block
            // and rendering attributes of any finished Slab
            while (nodeposZ < localSlabBlockStart[2]
                || nodeposZ >= localSlabBlockStart[2] + slabBlockSize[2]
                || nodeposY < localSlabBlockStart[1]
                || nodeposY >= localSlabBlockStart[1] + slabBlockSize[1]
                || nodeposX >= localSlabBlockStart[0] + slabBlockSize[0]) {
              // TODO: better implementation:
              // - no need to walk through all non existing "slab nodes"
              // - see if we cannot avoid some intermediate copies

              // keep slabs aligned on a regular grid
              // end of a block
              ++localSlabBlockIdxStart[2];
              localSlabBlockStart[2] += slabBlockSize[2];
              if (localSlabBlockStart[2] >= numSlabBlocksPerDim[2] * slabBlockSize[2]) {
                // end of tube
                localSlabBlockStart[2] = 0;
                if (nodeposY < localSlabBlockStart[1]
                    || nodeposY >= localSlabBlockStart[1] + slabBlockSize[1]
                    || nodeposX >= localSlabBlockStart[0] + slabBlockSize[0]) {
                  // next tube
                  localSlabBlockIdxStart[1] += numSlabBlocksPerDim[2]; // = localSlabBlockIdxStart[2];
                  localSlabBlockStart[1] += slabBlockSize[1];
                  if (localSlabBlockStart[1] >= numSlabBlocksPerDim[1] * slabBlockSize[1]) {
                    localSlabBlockStart[1] = 0;
                    if (nodeposX >= localSlabBlockStart[0] + slabBlockSize[0]) {
                      // output finished Slab
                      // encoder Hack:
                      //  pointcloud for complete slab is passed to attributes handling,
                      //  to be processed by blocks, then result is obtained from same
                      //  intermediate buffer
                      numPointsSlabBlocks.clear();
                      slabBlockIdx = localSlabBlockIdxStart[0];
                      int numPointsLocal = 0;
                      for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
                        for (int z = 0; z < numSlabBlocksPerDim[0]; ++z, ++slabBlockIdx) {
                          int numPoints = laPointIdx[slabBlockIdx].size();
                          numPointsSlabBlocks.emplace_back(numPoints);
                          if (numPoints) {
                            localPointCloud.appendPartition(pointCloud, laPointIdx[slabBlockIdx]);
                            laPointIdx[slabBlockIdx] = std::vector<int>(); // just release memory
                            numPointsLocal += numPoints;
                          }
                        }
                      rootEncoder.processNextSlabAttributes(localPointCloud,
                        {localSlabBlockStart[0], 0, 0}, numPointsSlabBlocks, false);
                      recPointCloud.setFromPartition(localPointCloud, 0, numPointsLocal, nRecPoints);
                      localPointCloud.clear();
                      nRecPoints += numPointsLocal;
                      // next slice
                      localSlabBlockIdxStart[0] += numSlabBlocksPerDim[1] * numSlabBlocksPerDim[2]; // = localSlabBlockIdxStart[2];
                      localSlabBlockStart[0] += slabBlockSize[0];
                    }
                    // restart from beginning of the slice
                    localSlabBlockIdxStart[2] = localSlabBlockIdxStart[0];
                    localSlabBlockIdxStart[1] = localSlabBlockIdxStart[0];
                  }
                }
                // restart from beginning of the tube
                localSlabBlockIdxStart[2] = localSlabBlockIdxStart[1];
              }
            }
            slabBlockIdx = localSlabBlockIdxStart[2];
          }

          int idxLaPoint = laPointIdx[slabBlockIdx].size();
          laPointIdx[slabBlockIdx].resize(idxLaPoint + node0.end - node0.start);
          for (auto idx = node0.start; idx < node0.end; idx++)
            laPointIdx[slabBlockIdx][idxLaPoint++] = idx;
        }

        for (int i = 0; i < 8; i++) {
          if (!node0.childCounts[i])
            continue;

          // if the bitstream is configured to represent unique points,
          // no point count is sent.
          if (!gps.geom_duplicated_points_flag) {
            assert(node0.childCounts[i] == 1);
            continue;
          }

          encoder.encodePositionLeafNumPoints(node0.childCounts[i]);
        }

        // leaf nodes do not get split
        return;
      }

      if (!isLeafNode(childSizeLog2)) {
        // push child nodes to fifo
        for (int i = 0; i < 2; ++i) {
          int childIndex = (nodeSliceIndex << 2) + (tubeIndex << 1) + i;

          bool occupiedChild = node0.childOccupancy & 1 << childIndex;
          if (!occupiedChild) {
            // child is empty: skip
            node0.predPointsStartIdx += node0.predCounts[childIndex];
          }
          else {
            // create new child
            fifoNext.emplace_back();
            auto& child = fifoNext.back();

            child.pos[0] = (node0.pos[0] << 1) + nodeSliceIndex;
            child.pos[1] = (node0.pos[1] << 1) + tubeIndex;
            child.pos[2] = (node0.pos[2] << 1) + i;

            int childPointsStartIdx = node0.start;
            for (int j = 0; j < childIndex; ++j)
              childPointsStartIdx += node0.childCounts[j];

            child.start = childPointsStartIdx;
            childPointsStartIdx += node0.childCounts[childIndex];
            child.end = childPointsStartIdx;

            child.predStart = node0.predPointsStartIdx;
            node0.predPointsStartIdx += node0.predCounts[childIndex];
            child.predEnd = node0.predPointsStartIdx;
            if (node0.mSONodeIdx >= 0) {
              child.mSONodeIdx = mSOctree.nodes[node0.mSONodeIdx].child[childIndex];
            }

            //status PU inheritance
            child.hasMotion = node0.hasMotion;
            child.isCompensated = node0.isCompensated;
            child.isSkiped = node0.isSkiped;
            child.quIndex = node0.quIndex;

            if (node0.hasMotion && !node0.isCompensated) {
              auto& puNode = mvField.puNodes[node0.mvFieldNodeIdx];
              assert(puNode._childsMask);
              assert(puNode._firstChildIdx != MVField::kNotSetMVIdx);
              child.mvFieldNodeIdx = puNode._firstChildIdx;
              for (int j = 0; j < childIndex; ++j)
                if (puNode._childsMask & (1 << j))
                  ++child.mvFieldNodeIdx;
            }

            if (isInter)
              if (!child.isSkiped) {
                child.mSOctreeNodeIdx = mSOctreeCurr.nodes[node0.mSOctreeNodeIdx].child[childIndex];
                assert(child.mSOctreeNodeIdx); // if not, mismatch between octrees
              } else {
                if (node0.mSOctreeNodeIdx || nodeSizeLog2 == gbh.rootNodeSizeLog2)
                  child.mSOctreeNodeIdx = mSOctreeCurr.nodes[node0.mSOctreeNodeIdx].child[childIndex];
              }

            if (forTrisoup) {
              if (isLastDepth) {
                trisoup->leaves.emplace_back(child);
                trisoup->leaves.back().pos *= S;
                if (flagNonPow2) {
                  int maskS = (1 << S2) - 1;
                  for (int np = child.start; np < child.end; np++) {
                    pointCloud[np][0] = ((pointCloud[np][0] >> S2) * S) + (pointCloud[np][0] & maskS);
                    pointCloud[np][1] = ((pointCloud[np][1] >> S2) * S) + (pointCloud[np][1] & maskS);
                    pointCloud[np][2] = ((pointCloud[np][2] >> S2) * S) + (pointCloud[np][2] & maskS);
                  }

                  for (int np = child.predStart; np < child.predEnd; np++) {
                    compensatedPointCloud[np][0] = ((compensatedPointCloud[np][0] >> S2) * S) + (compensatedPointCloud[np][0] & maskS);
                    compensatedPointCloud[np][1] = ((compensatedPointCloud[np][1] >> S2) * S) + (compensatedPointCloud[np][1] & maskS);
                    compensatedPointCloud[np][2] = ((compensatedPointCloud[np][2] >> S2) * S) + (compensatedPointCloud[np][2] & maskS);
                  }
                }

                if (lastPos0 != INT_MIN && child.pos[0] != lastPos0)
                  trisoup->sliceEncoder(false); // TriSoup unpile slices (not final = false)
                lastPos0 = child.pos[0];
              }
            }

            numNodesNextLvl++;
          }
        }
      }
    });

    fifo.resize(0);
    fifo.swap(fifoNext);

    // calculate the number of points that would be decoded if decoding were
    // to stop at this point.
    if (gps.octree_point_count_list_present_flag) {
      int numPtsAtLvl = numNodesNextLvl + nOutputPoints - 1;
      gbh.footer.octree_lvl_num_points_minus1.push_back(numPtsAtLvl);
    }
  }

  if (!(gps.inter_prediction_enabled_flag
    && gps.gof_geom_entropy_continuation_enabled_flag)
    && !(gps.trisoup_enabled_flag || gbh.entropy_continuation_flag)) {
    encoder.clearMap();
  }

  // the last element is the number of decoded points
  if (!gbh.footer.octree_lvl_num_points_minus1.empty())
    gbh.footer.octree_lvl_num_points_minus1.pop_back();

  if (motionEncoder.get()) {
    ctxtMemMotion = motionEncoder->getCtx();
  }

  rootEncoder.setMotionBits(fp16IOSizeMotionBits);

  if (forTrisoup) {
    // TriSoup final pass (true)
    trisoup->sliceEncoder(true);
    trisoup->finishSliceEncoder();
    return;
  }

  if (!useLocalAttr) {
    nRecPoints = nOutputPoints;
  } else {
  // The following is to render the remaining Slabs
    numPointsSlabBlocks.clear();
    int slabBlockIdx = localSlabBlockIdxStart[0];
    int numPointsLocal = 0;
    for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
      for (int z = 0; z < numSlabBlocksPerDim[2]; ++z, ++slabBlockIdx) {
        int numPoints = laPointIdx[slabBlockIdx].size();
        numPointsSlabBlocks.emplace_back(numPoints);
        if (numPoints) {
          localPointCloud.appendPartition(pointCloud, laPointIdx[slabBlockIdx]);
          laPointIdx[slabBlockIdx] = std::vector<int>(); // just release memory
          numPointsLocal += numPoints;
        }
      }
    rootEncoder.processNextSlabAttributes(localPointCloud,
      {localSlabBlockStart[0], 0, 0}, numPointsSlabBlocks, true);
    recPointCloud.setFromPartition(localPointCloud, 0, numPointsLocal, nRecPoints);
    localPointCloud.clear();
    nRecPoints += numPointsLocal;
  }
  recPointCloud.resize(nRecPoints);
  swap(pointCloud, recPointCloud);

  // The following is to render single slab-block for full slice
  if(!forTrisoup && !useLocalAttr) {
    std::vector<int32_t> numPoints {nRecPoints};
    rootEncoder.processNextSlabAttributes(pointCloud, {0, 0, 0}, numPoints, true);
  }
}

// instanciate for Trisoup
template void
encodeGeometryOctree<true>(
  const EncoderParams& encParams,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& rootEncoder,
  TrisoupEncoder* trisoup);

// instanciate for Octree
template void
encodeGeometryOctree<false>(
  const EncoderParams& encParams,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& rootEncoder,
  TrisoupEncoder* trisoup);

//============================================================================

void
encodeGeometryOctree(
  const EncoderParams& opt,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& encoder)
{
  encodeGeometryOctree<false>(
    opt, gps, gbh, pointCloud, ctxtMem, ctxtMemMotion, arithmeticEncoders,
    nullptr, refFrame, sps, interPredParams, encoder);
}

//-------------------------------------------------------------------------

//============================================================================
}  // namespace pcc
