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
#include "geometry_trisoup_decoder.h"
#include "io_hls.h"
#include "tables.h"
#include "PCCTMC3Decoder.h"
#include "motionWip.h"

namespace pcc {

//============================================================================

class GeometryOctreeDecoder {
public:
  GeometryOctreeContexts& ctx;

  GeometryOctreeDecoder(
    const GeometryParameterSet& gps,
    const GeometryBrickHeader& gbh,
    GeometryOctreeContexts& ctxMem,
    EntropyDecoder* arithmeticDecoder);

  GeometryOctreeDecoder(const GeometryOctreeDecoder&) = default;
  GeometryOctreeDecoder(GeometryOctreeDecoder&&) = default;

  // dynamic OBUF
  void clearMap() { ctx.clearMap(); };
  void resetMap(bool forTrisoup) { ctx.resetMap(forTrisoup); }

  int decodePositionLeafNumPoints();

  uint32_t decodeOccupancyFullNeihbourgs(
    const RasterScanContext::occupancy& occ,
    int predOcc[8],
    bool isInter,
    bool trisoup
  );

  // local QU
  void decodeQU(localQU& qu, int baseQP);

  const GeometryOctreeContexts& getCtx() const { return ctx; }

public:
  EntropyDecoder* _arithmeticDecoder;
};

//============================================================================

GeometryOctreeDecoder::GeometryOctreeDecoder(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  GeometryOctreeContexts& ctxtMem,
  EntropyDecoder* arithmeticDecoder)
  : ctx(ctxtMem)
  , _arithmeticDecoder(arithmeticDecoder)
{
}

//============================================================================
// Decode the number of points in a leaf node of the octree.
int
GeometryOctreeDecoder::decodePositionLeafNumPoints()
{
  int val = _arithmeticDecoder->decode(ctx._ctxDupPointCntGt0);
  if (val)
    val += _arithmeticDecoder->decodeExpGolomb(0, ctx._ctxDupPointCntEgl);

  return val + 1;
}

//-------------------------------------------------------------------------
// decode node occupancy bits

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
// decode node occupancy bits

uint32_t
GeometryOctreeDecoder::decodeOccupancyFullNeihbourgs(
  const RasterScanContext::occupancy& occ,
  int predOcc[8],
  bool isInter,
  bool trisoup)
{
  uint32_t occupancy = 0;
  int coded0 = 0;
  OctreeNeighours octreeNeighours;
  prepareGeometryAdvancedNeighPattern(occ, octreeNeighours);

  // loop on occupancy bits from occupancy map
  for (int i = 0; i < 8; i++) {
    if (coded0 >= 7) { // bit is 1
      occupancy += 1 << i;
      continue;
    }

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

    // decode
    int bit;
    if (Sparse) {
      if (trisoup)
        ctx2 >>= 4;

      bit = ctx._MapOccupancySparse[isInter][i].decodeEvolve(
        _arithmeticDecoder, ctx._CtxMapDynamicOBUF[ctxTable], ctx2, ctx1,
        &ctx._OBUFleafNumber, ctx._BufferOBUFleaves);
    }
    else {
      if (trisoup)
        ctx2 >>= 2;

      bit = ctx._MapOccupancy[isInter][i].decodeEvolve(
        _arithmeticDecoder, ctx._CtxMapDynamicOBUF[ctxTable], ctx2, ctx1,
        &ctx._OBUFleafNumber, ctx._BufferOBUFleaves);
    }

    // update partial occupancy of current node
    occupancy += bit << i;
    coded0 += !bit;
  }

  return occupancy;
}

//-------------------------------------------------------------------------
void
GeometryOctreeDecoder::decodeQU(localQU& qu, int baseQP)
{
  qu.isBaseParameters = _arithmeticDecoder->decode(ctx._ctxQUflag);
  if (qu.isBaseParameters) {
    qu.localQP = baseQP;
    return;
  }

  // sign
  bool sign = _arithmeticDecoder->decode(ctx._ctxQUSign);

  // magnitude
  int mag = 1
    + _arithmeticDecoder->decodeExpGolomb(1, ctx._ctxQUQPpref, ctx._ctxQUQPsuf);
  qu.localQP = sign ? mag : -mag;
  qu.localQP += baseQP;
}


//-------------------------------------------------------------------------

template <bool forTrisoup>
void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  EntropyDecoder& arithmeticDecoder,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  const Vec3<int> minimum_position,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& rootDecoder,
  TrisoupDecoder* trisoup)
{
  PCCPointSet3& compensatedPointCloud = interPredParams.compensatedPointCloud;
  PCCPointSet3& predPointCloud = interPredParams.referencePointCloud;
  MSOctree& mSOctree = interPredParams.mSOctreeRef;
  auto& mvField = interPredParams.mvField;
  std::unique_ptr<MotionEntropyDecoder> motionDecoder;
  const bool isInter = gbh.slice_inter_prediction_flag;

  int log2MotionBlockSize = 0;
  if (isInter) {
    int log2MinPUSize = ilog2(uint32_t(gps.motion.motion_min_pu_size - 1)) + 1;
    predPointCloud = refFrame->cloud;
    // for recoloring, need same depth as in encoder
    mSOctree = MSOctree(&predPointCloud, -gbh.geomBoxOrigin, std::min(2,log2MinPUSize));
    motionDecoder.reset(new MotionEntropyDecoder(ctxtMemMotion, &arithmeticDecoder));

    // local motion prediction structure -> LPUs from predPointCloud
    log2MotionBlockSize = ilog2(uint32_t(gps.motion.motion_block_size - 1)) + 1;
    // LPU is bigger than root node, must adjust
    if (gbh.rootNodeSizeLog2 < log2MotionBlockSize)
      log2MotionBlockSize = gbh.rootNodeSizeLog2;
  }

  // local QU
  if (forTrisoup) {
    if (gbh.qu_size_log2 > 0) {
      int quSize = (1 << gbh.qu_size_log2) * gbh.trisoupNodeSize(gps);
      std::cout << "QU size = " << quSize << "\n";
    }
    ctxtMem.quLastIndex = 0;
    ctxtMem.listOfQUs.clear();
  }

  const int S = gbh.trisoupNodeSize(gps);
  const int S2 = gbh.trisoupNodeSizeLog2(gps);
  const int factorS = (1 << S2) - S;
  const bool flagNonPow2 = factorS != 0;

  // init main fifo and reserve for max size
  int ringBufferSize = gbh.footer.geom_num_points_minus1 + 1;
  if (gps.trisoup_enabled_flag && gbh.trisoupNodeSize(gps) > 1)
     ringBufferSize = std::max(1000,ringBufferSize >> 2 * gbh.trisoupNodeSizeLog2(gps) - 2);
  std::vector<PCCOctree3Node> fifo;
  std::vector<PCCOctree3Node> fifoNext;
  fifo.reserve(ringBufferSize);
  fifoNext.reserve(ringBufferSize);

  // saved state for use with parallel bistream coding.
  // the saved state is restored at the start of each parallel octree level
  std::unique_ptr<GeometryOctreeContexts> savedState;

  // push the first node
  fifo.emplace_back();
  PCCOctree3Node& node00 = fifo.back();
  node00.start = uint32_t(0);
  node00.end = uint32_t(0);
  node00.pos = int32_t(0);
  node00.predStart = uint32_t(0);
  node00.predEnd = isInter ? predPointCloud.getPointCount() : uint32_t(0);
  node00.mSONodeIdx = isInter ? 0 : -1;
  node00.hasMotion = 0;
  node00.isCompensated = 0;
  node00.isSkiped = false;
  node00.quIndex = -1;
  size_t processedPointCount = 0;

  RasterScanContext rsc(fifo);

  // the termination depth of the octree phase
  int maxDepth = gbh.tree_depth_minus1 + 1;

  // Derived parameter used by local attributes
  if (!forTrisoup) {
    gbh.maxRootNodeSize = 1 << maxDepth;
    gbh.uniqueSlabBlock = gbh.uniqueSlabBlock && gbh.maxRootNodeSize <= std::min(
      sps.localized_attributes_slab_thickness_minus1 + 1, sps.localized_attributes_slab_block_size_minus1 + 1);
  }

  // NB: this needs to be after the root node size is determined to
  GeometryOctreeDecoder decoder(gps, gbh, ctxtMem, &arithmeticDecoder);

  if (!(isInter && gps.gof_geom_entropy_continuation_enabled_flag) && !gbh.entropy_continuation_flag) {
    decoder.clearMap();
    decoder.resetMap(forTrisoup);
  }

  // localized attributes point indexes
  std::vector<std::vector<point_t> > laPoints;
  point_t localSlabBlockIdxStart = 0;
  point_t numSlabBlocksPerDim = -1;
  point_t localSlabBlockStart = 0;
  point_t slabBlockSize;
  bool useLocalAttr = !forTrisoup && !gbh.uniqueSlabBlock;
  PCCPointSet3 localPointCloud;
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
    // memory allocation for worst case
    int numSlabBlocksPerSlab = numSlabBlocksPerDim[1] * numSlabBlocksPerDim[2];
    int numSlabBlocks = numSlabBlocksPerDim[0] * numSlabBlocksPerSlab;
    laPoints.resize(numSlabBlocks);
    localPointCloud.addRemoveAttributes(pointCloud);
    localPointCloud.reserve(ringBufferSize);
  }

  // ----------  loop on depth ----------------
  int lastPos0 = INT_MIN; // used to detect slice change and call for TriSoup
  int nodeSizeLog2 = gbh.rootNodeSizeLog2;
  for (int depth = 0; depth < maxDepth; depth++, nodeSizeLog2--) {
    // setup at the start of each level
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
        savedState.reset(new GeometryOctreeContexts(decoder.ctx));

    // a new entropy stream starts one level after the context state is saved.
    // restore the saved state and flush the arithmetic decoder
    if (depth > maxDepth - 1 - gbh.geom_stream_cnt_minus1) {
      decoder.ctx = *savedState;
      arithmeticDecoder.flushAndRestart();
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
        // decode local motion PU tree
        if (isInter) {
          if (nodeSizeLog2 == log2MotionBlockSize) {
            node0.mvFieldNodeIdx = currPUIdx++;
            node0.hasMotion = true;
          }

          // decode LPU/PU/MV
          if (node0.hasMotion && !node0.isCompensated) {
            const int  nodeSize = S * (1 << nodeSizeLog2 - S2);
            decode_splitPU_MV_MC<false>(mSOctree, nullptr,
              &node0, mvField, node0.mvFieldNodeIdx,
              gps.motion, nodeSize, *motionDecoder,
              &compensatedPointCloud, flagNonPow2, S, S2);
          }
        }

        // ---------- skip Mode ---------------
        if (gps.trisoup_early_skip_coding_mode_enabled_flag && node0.isCompensated
          && nodeSizePow2 == gps.trisoup_early_skip_coding_mode_node_size
          && node0.predStart < node0.predEnd) {
          node0.isSkiped = arithmeticDecoder.decode(ctxtMem.ctxSkipMode) == 1;
        }

        // ...for local motion
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
            decoder.decodeQU(qu, gbh.trisoup_QP);
            node0.quIndex = ctxtMem.quLastIndex++;
          }
        }

        // occupancy of child nodes
        if (!forTrisoup || !node0.isSkiped) {
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

          // decode child occupancy map
          node0.childOccupancy = decoder.decodeOccupancyFullNeihbourgs(
            contextualOccupancy, predOccupancy, isInter, forTrisoup);
        }
        else {
          int predOcc = 0;
          for (int i = 0; i < 8; i++)
            predOcc |= (node0.predCounts[i] > 0) << i;
          node0.childOccupancy = predOcc;
        }

        // create motion field child nodes if needed
        if (isInter && node0.hasMotion && !node0.isCompensated) {
          auto& puNode = mvField.puNodes[node0.mvFieldNodeIdx];
          uint8_t occupancy = node0.childOccupancy;
          puNode._childsMask = occupancy;
          puNode._firstChildIdx = mvField.puNodes.size();
          int numOccupied = popcnt(occupancy);
          for (int puChildIdx = 0; puChildIdx < numOccupied; ++puChildIdx) {
            mvField.puNodes.emplace_back();
          }
        }

      } //end firstVisit

      uint8_t occupancy = node0.childOccupancy;

      if (!forTrisoup && tubeIndex && nodeSliceIndex && isLeafNode(childSizeLog2)) {
        int slabBlockIdx;
        if (useLocalAttr) {
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
                      slabBlockIdx = localSlabBlockIdxStart[0];
                      for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
                        for (int z = 0; z < numSlabBlocksPerDim[0]; ++z, ++slabBlockIdx) {
                          int numPoints = laPoints[slabBlockIdx].size();
                          if (numPoints) {
                            localPointCloud.resize(numPoints);
                            for (int i = 0; i < numPoints; ++i)
                              localPointCloud[i] = laPoints[slabBlockIdx][i];
                            laPoints[slabBlockIdx] = std::vector<point_t>(); // just release memory
                            rootDecoder.processNextSlabBlockAttributes(localPointCloud,
                              {localSlabBlockStart[0], y * slabBlockSize[1], z * slabBlockSize[2]});
                            pointCloud.setFromPartition(localPointCloud, 0, numPoints, processedPointCount);
                            localPointCloud.clear();
                            processedPointCount += numPoints;
                          }
                        }
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
        }

        for (int i = 0; i < 8; i++) {
          uint32_t mask = 1 << i;
          if (!(occupancy & mask))
            continue;

          // point counts for leaf nodes are coded immediately
          int numPoints = 1;
          if (gps.geom_duplicated_points_flag) {
            numPoints = decoder.decodePositionLeafNumPoints();
          }

          // the final bits from the leaf:
          Vec3<int32_t> point{
            (node0.pos[0] << 1) + !!(i & 4),
            (node0.pos[1] << 1) + !!(i & 2),
            (node0.pos[2] << 1) + !!(i & 1) };

          if (useLocalAttr) {
            int idxLaPoint = laPoints[slabBlockIdx].size();
            laPoints[slabBlockIdx].resize(idxLaPoint + numPoints);
            for (int i = 0; i < numPoints; ++i)
              laPoints[slabBlockIdx][idxLaPoint++] = point;
          } else {
            for (int i = 0; i < numPoints; ++i)
              pointCloud[processedPointCount++] = point;
          }
        }

        // do not recurse into leaf nodes
        return;
      }


      if (!isLeafNode(childSizeLog2)) {
        // push child nodes to fifo
        for (int i = 0; i < 2; ++i) {
          int childIndex = (nodeSliceIndex << 2) + (tubeIndex << 1) + i;
          uint32_t mask = 1 << childIndex;
          bool occupiedChild = occupancy & mask;
          if (!occupiedChild) {
            // child is empty: skip
            node0.predPointsStartIdx += node0.predCounts[childIndex];
          }
          else {
            // create & enqueue new child.
            fifoNext.emplace_back();
            auto& child = fifoNext.back();

            child.pos[0] = (node0.pos[0] << 1) + nodeSliceIndex;
            child.pos[1] = (node0.pos[1] << 1) + tubeIndex;
            child.pos[2] = (node0.pos[2] << 1) + i;

            child.predStart = node0.predPointsStartIdx;
            node0.predPointsStartIdx += node0.predCounts[childIndex];
            child.predEnd = node0.predPointsStartIdx;
            if (node0.mSONodeIdx >= 0) {
              child.mSONodeIdx = mSOctree.nodes[node0.mSONodeIdx].child[childIndex];
            }

            //status inheritance
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

            if (forTrisoup) {
              if (isLastDepth) {
                trisoup->leaves.emplace_back(child);
                trisoup->leaves.back().pos *= S;

                if (flagNonPow2) {
                  int maskS = (1 << S2) - 1;
                  for (int np = child.predStart; np < child.predEnd; np++) {
                    compensatedPointCloud[np][0] = ((compensatedPointCloud[np][0] >> S2) * S) + (compensatedPointCloud[np][0] & maskS);
                    compensatedPointCloud[np][1] = ((compensatedPointCloud[np][1] >> S2) * S) + (compensatedPointCloud[np][1] & maskS);
                    compensatedPointCloud[np][2] = ((compensatedPointCloud[np][2] >> S2) * S) + (compensatedPointCloud[np][2] & maskS);
                  }
                }

                if (lastPos0 != INT_MIN && child.pos[0] != lastPos0)
                  trisoup->sliceDecoder(false); // TriSoup unpile slices (not final = false)
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
  }
  if (!(gps.inter_prediction_enabled_flag
        && gps.gof_geom_entropy_continuation_enabled_flag)
        && !(gps.trisoup_enabled_flag || gbh.entropy_continuation_flag))
    decoder.clearMap();

  if (motionDecoder.get()) {
    ctxtMemMotion = motionDecoder->getCtx();
  }

  // The following is to render the remaining Slabs
  if (useLocalAttr) {
    int slabBlockIdx = localSlabBlockIdxStart[0];
    for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
      for (int z = 0; z < numSlabBlocksPerDim[2]; ++z, ++slabBlockIdx) {
        int numPoints = laPoints[slabBlockIdx].size();
        if (numPoints) {
          localPointCloud.resize(numPoints);
          for (int i = 0; i < numPoints; ++i)
            localPointCloud[i] = laPoints[slabBlockIdx][i];
          laPoints[slabBlockIdx] = std::vector<point_t>(); // just release memory
          rootDecoder.processNextSlabBlockAttributes(localPointCloud,
            {localSlabBlockStart[0], y * slabBlockSize[1], z * slabBlockSize[2]});
          pointCloud.setFromPartition(localPointCloud, 0, numPoints, processedPointCount);
          localPointCloud.clear();
          processedPointCount += numPoints;
        }
      }
  }

  // resize point cloud (not needed for TriSoup)
  if (!forTrisoup)
    pointCloud.resize(processedPointCount);

  // TriSoup final pass (true)
  if (forTrisoup) {
    trisoup->sliceDecoder(true);
    trisoup->finishSliceDecoder();
  }

  // The following is to render single slab-block for full slice
  if(!forTrisoup && !useLocalAttr) {
    rootDecoder.processNextSlabBlockAttributes(pointCloud, {0, 0, 0});
  }
}

// instanciate for Trisoup
template void
decodeGeometryOctree<true>(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  EntropyDecoder& arithmeticDecoder,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  const Vec3<int> minimum_position,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& rootDecoder,
  TrisoupDecoder* trisoup);

// instanciate for Octree
template void
decodeGeometryOctree<false>(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int skipLastLayers,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  EntropyDecoder& arithmeticDecoder,
  std::vector<PCCOctree3Node>* nodesRemaining,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  const Vec3<int> minimum_position,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& rootDecoder,
  TrisoupDecoder* trisoup);

//-------------------------------------------------------------------------

void
decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  const Vec3<int> minimum_position,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& decoder
)
{
  decodeGeometryOctree<false>(
    gps, gbh, 0, pointCloud, ctxtMem, ctxtMemMotion, arithmeticDecoder, nullptr,
    refFrame, sps, minimum_position, interPredParams, decoder, nullptr);
}

//-------------------------------------------------------------------------

void
decodeGeometryOctreeScalable(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  int minGeomNodeSizeLog2,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMem,
  MotionEntropy& ctxtMemMotion,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  PCCTMC3Decoder3& decoder
)
{
  std::vector<PCCOctree3Node> nodes;
  InterPredParams interPredParams;
  decodeGeometryOctree<false>(
    gps, gbh, minGeomNodeSizeLog2, pointCloud, ctxtMem, ctxtMemMotion, arithmeticDecoder,
    &nodes, refFrame, sps, { 0, 0, 0 }, interPredParams, decoder);

  if (minGeomNodeSizeLog2 > 0) {
    size_t size =
      pointCloud.removeDuplicatePointInQuantizedPoint(minGeomNodeSizeLog2);

    pointCloud.resize(size + nodes.size());
    size_t processedPointCount = size;

    if (minGeomNodeSizeLog2 > 1) {
      uint32_t mask = uint32_t(-1) << minGeomNodeSizeLog2;
      for (auto node0 : nodes) {
        for (int k = 0; k < 3; k++)
          node0.pos[k] &= mask;
        node0.pos += 1 << (minGeomNodeSizeLog2 - 1);
        pointCloud[processedPointCount++] = node0.pos;
      }
    } else {
      for (const auto& node0 : nodes)
        pointCloud[processedPointCount++] = node0.pos;
    }
  }
}

//============================================================================

}  // namespace pcc
