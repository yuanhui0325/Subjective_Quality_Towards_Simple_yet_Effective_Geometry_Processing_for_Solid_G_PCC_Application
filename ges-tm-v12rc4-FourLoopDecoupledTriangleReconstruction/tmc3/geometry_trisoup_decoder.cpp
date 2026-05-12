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

#include <cstdio>
#include <queue>
#include "geometry_trisoup_decoder.h"
#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"
#include "PCCTMC3Decoder.h"

namespace pcc {

//============================================================================

void
decodeGeometryTrisoup(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  MotionEntropy& ctxtMemMotion,
  EntropyDecoder& arithmeticDecoder,
  const CloudFrame* refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Decoder3& decoder)
{
  const Vec3<int> minimum_position = sps.seqBoundingBoxOrigin;
  bool isInter = gbh.slice_inter_prediction_flag;

  // prepare TriSoup parameters
  int blockWidth = gbh.trisoupNodeSize(gps);
  std::cout << "TriSoup QP = " << gbh.trisoup_QP << "\n";

  point_t slabBlockSize = {
    sps.localized_attributes_slab_thickness_minus1 + 1,
    sps.localized_attributes_slab_block_size_minus1 + 1,
    sps.localized_attributes_slab_block_size_minus1 + 1
  };

  // used to determine if full frame optimized implementation is used or not
  {
    // Derived parameter used by trisoup
    gbh.maxRootNodeSize = gbh.trisoupNodeSize(gps) << gbh.tree_depth_minus1 + 1; // get bounding box of point cloud geometry
    gbh.uniqueSlabBlock = gbh.uniqueSlabBlock && gbh.maxRootNodeSize <= std::min(
      sps.localized_attributes_slab_thickness_minus1 + 1, sps.localized_attributes_slab_block_size_minus1 + 1);
  }

  // trisoup uses octree coding until reaching the triangulation level.
  TrisoupDecoder trisoup(blockWidth, pointCloud,
    1 /*distanceSearchEncoder*/, isInter, interPredParams.compensatedPointCloud,
    gps, gbh, &arithmeticDecoder, ctxtMemOctree,
    !gbh.uniqueSlabBlock, slabBlockSize, &decoder);

  trisoup.initDecoder();

  trisoup.thVertexDeterminationSkip = gbh.trisoup_early_skip_vertex_determination_threshold;

  // octree
  decodeGeometryOctree<true>(
    gps, gbh, 0, pointCloud, ctxtMemOctree, ctxtMemMotion, arithmeticDecoder,
    nullptr, refFrame, sps, minimum_position, interPredParams, decoder, &trisoup);

  std::cout << "\nSize compensatedPointCloud for TriSoup = "
    << interPredParams.compensatedPointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes for TriSoup = " << trisoup.leaves.size() << "\n";
  std::cout << "TriSoup gives " << pointCloud.getPointCount() << " points \n";
}

//---------------------------------------------------------------------------
//  decoding of centroid residual in TriSOup node
int
decodeCentroidResidual(
  pcc::EntropyDecoder* arithmeticDecoder,
  GeometryOctreeContexts& ctxtMemOctree,
  int resCentroQPred,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound)
{
  int resCentroQ = 1;

  if (resCentroQPred == -100) //intra
    resCentroQ = arithmeticDecoder->decode(ctxtMemOctree.ctxRes0[ctxMinMax][0]) ? 0 : 1;
  else {//inter
    resCentroQ = arithmeticDecoder->decode(ctxtMemOctree.ctxRes0[ctxMinMax][1 + std::min(3, std::abs(resCentroQPred))]) ? 0 : 1;
  }
  if (resCentroQ == 0)
    return 0;

  // if not 0, residual in [-lowBound; highBound]
  // code sign
  int sign = 1;
  if (highBound && lowBound) {// otherwise sign is knwow
    int lowS = std::min(7, lowBoundSurface);
    int highS = std::min(7, highBoundSurface);
    sign = arithmeticDecoder->decode(ctxtMemOctree.ctxResSign[lowBound == highBound ? 0 : 1 + (lowBound < highBound)][lowS][highS][(resCentroQPred && resCentroQPred != -100) ? 1 + (resCentroQPred > 0) : 0]);
  }
  else if (!highBound) // highbound is 0 , so sign is negative; otherwise sign is already set to positive
    sign = 0;

  // code remaining bits 1 to 7 at most
  int magBound = (sign ? highBound : lowBound) - 1;
  bool sameSignPred = resCentroQPred != -100 && (resCentroQPred > 0 && sign) || (resCentroQPred < 0 && !sign);

  int ctx = 0;
  int ctx2 = resCentroQPred != -100 ? 1 + std::min(8, sameSignPred * std::abs(resCentroQPred)) : 0;
  while (magBound > 0) {
    int bit;
    if (ctx < 4)
      bit = arithmeticDecoder->decode(ctxtMemOctree.ctxResMag[ctx][ctx2]);
    else
      bit = arithmeticDecoder->decode();

    if (bit) // magRes==0 and magnitude coding is finished
      break;

    resCentroQ++;
    magBound--;
    ctx++;
  }

  return sign ? resCentroQ : -resCentroQ;
}

//============================================================================

void
TrisoupDecoder
::decodeOneTriSoupVertexRasterScan(
  pcc::EntropyDecoder* const arithmeticDecoder,
  GeometryOctreeContexts& ctxtMemOctree,
  std::vector<int8_t>& TriSoupVerticesQP,
  std::vector<int8_t>& TriSoupVertices2bits,
  int neigh,
  std::array<int, 18>& patternIdx,
  int8_t interPredictor,
  int8_t colocatedVertex,
  int i,
  int qpEdge)
{
  int maxMag = LUT_maxMagnitude[qpEdge];
  int nBitMag = LUT_nBitMagnitude[qpEdge];

  codeVertexCtxInfo ctxInfo;
  constructCtxInfo(ctxInfo, neigh, patternIdx, TriSoupVertices2bits);

  // decode vertex presence
  int ctxMap1, ctxMap2, ctxMap2Base, ctxInter;
  int shiftMag = std::max(0, nBitMag - 1);
  int interPredictor2bit = 0;
  if (nBitMag >= 1) {
    int bit1 = interPredictor > 0;
    int bit2 = ((std::abs(interPredictor) - 1) >> nBitMag - 1) & 1;
    interPredictor2bit = (bit1 << 1) + (bit1 ? bit2 : !bit2);
  }

  constructCtxPresence(ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, colocatedVertex);
  bool vertexPresent = ctxtMemOctree.MapOBUFTriSoup[0].decodeEvolve(
    arithmeticDecoder, ctxtMemOctree.ctxTriSoup[0][ctxInter], ctxMap2,
    ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
    ctxtMemOctree._BufferOBUFleavesTrisoup);

  // decode vertex position
  if (vertexPresent) {
    uint8_t partialDec = 0;
    int bitRemaining = nBitMag - 1;
    int blockWidthLog2 = ilog2(uint32_t(blockWidth));

    // first position bit is for left vs right  -> sign of vertexQP
    constructCtxPos1(
      ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor,
      colocatedVertex, interPredictor2bit);

    int bit = ctxtMemOctree.MapOBUFTriSoup[1].decodeEvolve(
      arithmeticDecoder, ctxtMemOctree.ctxTriSoup[1][ctxInter], ctxMap2,
      ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
      ctxtMemOctree._BufferOBUFleavesTrisoup);
    bool vertexSign = bit;
    partialDec = bit;

    // second position bit
    int magnitude = 0;
    int goodColo = 0;
    if (bitRemaining >= 0) {
      constructCtxPos2(
        ctxMap2Base, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor,
        shiftMag, partialDec, colocatedVertex, blockWidthLog2, goodColo);

      bit = 0;
      if (isCodedBit(magnitude, bitRemaining, maxMag))
        bit = ctxtMemOctree.MapOBUFTriSoup[2].decodeEvolve(
          arithmeticDecoder, ctxtMemOctree.ctxTriSoup[2][ctxInter], ctxMap2,
          (ctxMap1 << 1) + partialDec, &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);

      partialDec = (partialDec << 1) | bit;
      magnitude = (magnitude << 1) | bit;
      bitRemaining--;
    }

    // third bit
    if (bitRemaining >= 0) {
      constructCtxPos3(ctxMap2Base, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, shiftMag, partialDec, colocatedVertex, goodColo);

      bit = 0;
      if (isCodedBit(magnitude, bitRemaining, maxMag))
        bit = ctxtMemOctree.MapOBUFTriSoup[3].decodeEvolve(
          arithmeticDecoder, ctxtMemOctree.ctxTriSoup[3][ctxInter], ctxMap2,
          (ctxMap1 << 2) + partialDec, &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);

      partialDec = (partialDec << 1) | bit;
      magnitude = (magnitude << 1) | bit;
      bitRemaining--;
    }

    // forth bit
    if (bitRemaining >= 0) {
      constructCtxPos4(ctxMap2Base, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, shiftMag, partialDec, colocatedVertex, goodColo);

      bit = 0;
      if (isCodedBit(magnitude, bitRemaining, maxMag))
        bit = ctxtMemOctree.MapOBUFTriSoup[4].decodeEvolve(
          arithmeticDecoder, ctxtMemOctree.ctxTriSoup[4][ctxInter], ctxMap2,
          (ctxMap1 << 3) + partialDec, &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);

      magnitude = (magnitude << 1) | bit;
      bitRemaining--;
    }

    // remaining bits are bypassed
    for (; bitRemaining >= 0; bitRemaining--) {
      bit = 0;
      if (isCodedBit(magnitude, bitRemaining, maxMag))
        bit = arithmeticDecoder->decode();
      magnitude = (magnitude << 1) | bit;
    }

    TriSoupVerticesQP[firstVertexToCode] = vertexSign ? 1 + magnitude : -1 - magnitude;
  }
}

//----------------------------------------------------------------------------

void
TrisoupDecoder
::decodeFlags(
  pcc::EntropyDecoder* const decoder,
  uint8_t& axiFlag,
  Vec3<bool>& posFlag,
  AdaptiveBitModel (&axiFlagctx)[3][2],
  AdaptiveBitModel (&posFlagctx)[3][2],
  Vec3<int>& axiFlag_ctx,
  Vec3<int>& posFlag_ctx,
  int idx)
{
  if (decoder->decode(axiFlagctx[idx][axiFlag_ctx[idx]])) {
    axiFlag |= 1 << idx;
    posFlag[idx] = decoder->decode(posFlagctx[idx][posFlag_ctx[idx]]);
  }
}

//----------------------------------------------------------------------------

void
TrisoupDecoder
::prepareVertexConsistencyDecoder(
  pcc::EntropyDecoder* const arithmeticDecoder,
  TrisoupOriginalPCinfo& oriPCinfo)
{
  int cnt = 0;
  Vec3<int> poscount = 0;
  for (int j = 0; j < 12; j++) {
    int uniqueIndex = segmentUniqueIndex[idxSegment + j];
    int vertexQ = TriSoupVerticesQP[uniqueIndex];
    if (vertexQ) {// skip segments that do not intersect the surface
      cnt++;
      if (vertexQ < 0) {
        poscount[0] += posTableQNegative[0][j];
        poscount[1] += posTableQNegative[1][j];
        poscount[2] += posTableQNegative[2][j];
      }
      poscount[0] += posTableQPositive[0][j];
      poscount[1] += posTableQPositive[1][j];
      poscount[2] += posTableQPositive[2][j];
    }
  }
  int eligible = 0;
  const int numth = blockWidth >= 10 ? 6 : 5;
  if (cnt > numth) {
    for (int k = 0; k < 3; k++) {
      if (poscount[k] && poscount[k] != cnt
        && (poscount[k] >= cnt - eligible_th[cnt - 6] || poscount[k] <= eligible_th[cnt - 6]))
        eligible |= 1 << k;
    }
  }
  oriPCinfo.eligible = eligible;

  if (!eligible)
    return;

  uint8_t axisFlagMask = 0;
  Vec3<bool> posFlag = false;
  if (arithmeticDecoder->decode(multiFlagctx)) {

    Vec3<int> axiFlag_ctx;
    Vec3<int> posFlag_ctx;
    for (int k = 0; k < 3; k++) {
      axiFlag_ctx[k] = (poscount[k] >= (cnt - 1) || poscount[k] <= 1) ? 1 : 0;
      posFlag_ctx[k] = poscount[k] >= (cnt - 1) ? 1 : 0;
    }

    if (eligible > 4) {
      decodeFlags(arithmeticDecoder, axisFlagMask, posFlag, axiFlagctx, posFlagctx, axiFlag_ctx, posFlag_ctx, 2);
    }
    if ((eligible & 3) > 2) {
      decodeFlags(arithmeticDecoder, axisFlagMask, posFlag, axiFlagctx, posFlagctx, axiFlag_ctx, posFlag_ctx, 1);
    }
    int b = !(eligible & 1) + (eligible == 4);
    axisFlagMask |= (!axisFlagMask || arithmeticDecoder->decode(axiFlagctx[b][axiFlag_ctx[b]])) << b;
    if (axisFlagMask & 1 << b)
      posFlag[b] = arithmeticDecoder->decode(posFlagctx[b][posFlag_ctx[b]]);
  }
  oriPCinfo.axisFlagMask = axisFlagMask;
  oriPCinfo.posFlag = posFlag;
}

//----------------------------------------------------------------------------
void
TrisoupDecoder
::CentroidsResidualDecoder(
  int stepQcentro,
  int dominantAxis,
  int scaleQ,
  Vec3<int32_t>& blockCentroid,
  int vtxCount,
  const TrisoupNodeDecoder& leaf,
  std::vector<Vec3<int32_t>>& leafVertices,
  Vec3<int32_t>& neiCentroid,
  bool eligDoubleCentro,
  bool& isDoubleCentro,
  TrisoupNodeEdgeVertex& neVertex,
  uint16_t subCentroMask,
  std::array<Vec3<int32_t>, 2>& subCentroid
  )
{
  int half = stepQcentro >> 1;
  int DZ = 682 * half >> 10; // 2 * half / 3;

  int lowBound, highBound, lowBoundSurface, highBoundSurface, ctxMinMax;
  Vec3<int32_t> normalV =
    determineCentroidNormalAndBounds(
      lowBound, highBound, lowBoundSurface, highBoundSurface, ctxMinMax,
      vtxCount, blockCentroid, dominantAxis,
      leafVertices, blockWidth, stepQcentro, scaleQ, neiCentroid);

  bool flagCentroCoding = lowBound != 0 || highBound != 0;
  int resCentroQ = 0;

  if (flagCentroCoding) {
    int resCentroQPred = 0;
    determineCentroidPredictor(
      resCentroQPred, normalV, blockCentroid, leaf.pos,
      compensatedPointCloud, leaf.predStart, leaf.predEnd, lowBound,
      highBound, stepQcentro, scaleQ, vtxCount);

    // decode centroid residual
    resCentroQ = decodeCentroidResidual(arithmeticDecoder, ctxtMemOctree, resCentroQPred, ctxMinMax,
        lowBoundSurface, highBoundSurface, lowBound, highBound);
    if (eligDoubleCentro && resCentroQ == 0)
      isDoubleCentro = arithmeticDecoder->decode(ctxtMemOctree.ctxDoubleCentroid);

    if (isDoubleCentro) {
      std::vector<Vec3<int32_t>> vertices0;
      std::vector<Vec3<int32_t>> vertices1;
      for (int j = 0; j < vtxCount; j++) {
        if (subCentroMask & 1 << j)
          vertices0.push_back(neVertex.vertices[j]);
        else
          vertices1.push_back(neVertex.vertices[j]);
      }
      subCentroid[0] = getSubCentroid(vertices0);
      subCentroid[1] = getSubCentroid(vertices1);

      return;
    }
  }

  // dequantize and apply residual
  int resCentroDQ = 0;
  if (resCentroQ)
    resCentroDQ = std::abs(resCentroQ) * stepQcentro + DZ - half;

  if (resCentroQ < 0)
    resCentroDQ = -resCentroDQ;

  resCentroDQ = resCentroDQ > 0 ? (resCentroDQ * scaleQ >> 8) : -((-resCentroDQ) * scaleQ >> 8);
  blockCentroid += (resCentroDQ * normalV) >> 8;

  int boundW = ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1;
  blockCentroid[0] = std::min(boundW, std::max(-kTrisoupFpHalf, blockCentroid[0]));
  blockCentroid[1] = std::min(boundW, std::max(-kTrisoupFpHalf, blockCentroid[1]));
  blockCentroid[2] = std::min(boundW, std::max(-kTrisoupFpHalf, blockCentroid[2]));
}


//----------------------------------------------------------------------------
void
TrisoupDecoder
::GenerateCentroidsInNodeRasterScanDecoder(
  const TrisoupNodeDecoder& leaf,
  const std::vector<int8_t>& TriSoupVerticesQP,
  const std::vector<int8_t>& TriSoupEdgeLocalQP,
  int& idxSegment,
  const bool isCentroidResActivated,
  std::vector<int>& segmentUniqueIndex,
  TrisoupOriginalPCinfo* poriPCinfo,
  int nodeIdx)
{
  // Find up to 12 vertices for this leaf in neVertex.vertices
  eVerts.emplace_back();
  TrisoupNodeEdgeVertex& neVertex = eVerts.back();
  cVerts.emplace_back(TrisoupCentroidVertex{ false, { 0, 0, 0 }, false, {}, 0, true, { 0, 0, 0 } });
  TrisoupCentroidVertex& ncVertex = cVerts.back();
  isNonClosedForwardNode.push_back(false);

  TrisoupOriginalPCinfo& oriPCinfo = *poriPCinfo;
  std::array<bool, 8> eligibleMergedCorner{};
  std::array<bool, 8> isMergedCorner{};
  constexpr int M = 4;
  constexpr int N = 6;
  int vertexCnt = 0;
  int32_t th = std::min(8 << kTrisoupFpBits, blockWidth << (kTrisoupFpBits - 2));
  int32_t bmax = blockWidth << kTrisoupFpBits;
  if (vertexMergeFlag) {
    for (int n = 0; n < 8; n++) {
      Vec3< int64_t> posKey = leaf.pos + corner[n];
      const int64_t shift = 1 << 18;
      int64_t keyEdge = (posKey[0] + shift << 40) + (posKey[1] + shift << 20) + posKey[2] + shift;
      auto it = triNodeEnable.find(keyEdge);
      eligibleMergedCorner[n] = it != triNodeEnable.end() && it->second >= M;
    }

    int cnt = idxSegment;
    for (int j = 0; j < 12; j++)
      vertexCnt += TriSoupVerticesQP[segmentUniqueIndex[cnt++]] != 0;
  }
  uint8_t per_oriPointCloud = 255;
  if (oriPCinfo.eligible) {
    for (int k = 0; k < 3; k++) {
      if (oriPCinfo.axisFlagMask & 1 << k) {
        per_oriPointCloud = per_oriPointCloud & oriPointCloudMask[oriPCinfo.posFlag[k]][k];
      }
    }
  }
  for (int j = 0; j < 12; j++) {
    int uniqueIndex = segmentUniqueIndex[idxSegment++];
    int vertexQ = TriSoupVerticesQP[uniqueIndex];

    bool flagContinue = vertexQ == 0;
    flagContinue = flagContinue || (vertexConsistencyFlag && oriPCinfo.eligible
      && !(per_oriPointCloud & faceIdxToSubBlockMask[LUTSegmentExtToFaceIdx[j][0]])
      || !(per_oriPointCloud & faceIdxToSubBlockMask[LUTSegmentExtToFaceIdx[j][1]]));
    if (flagContinue)
      continue;  // skip segments that do not intersect the surface

    int qpEdge = TriSoupEdgeLocalQP[uniqueIndex];
    int deQVert = dequantizeQP(vertexQ, qpEdge);
    if (vertexMergeFlag && vertexCnt >= N) {
      int32_t deQVert1 = deQVert + kTrisoupFpHalf;
      if (eligibleMergedCorner[vertexIdx1[j]] && deQVert1 < th) {
        isMergedCorner[vertexIdx1[j]] = true;
        continue;
      }
      int32_t deQVert2 = bmax - deQVert1;
      if (eligibleMergedCorner[vertexIdx2[j]] && deQVert2 < th) {
        isMergedCorner[vertexIdx2[j]] = true;
        continue;
      }
    }
    // Get 3D position of point of intersection.
    Vec3<int32_t> point = (corner[startCorner[j]] << kTrisoupFpBits) - kTrisoupFpHalf; // the volume is [-0.5; B-0.5]^3
    point[LUTsegmentDirection[j]] = deQVert;
    // Add vertex to list of vertices.
    neVertex.vertices.push_back(point);
  }

  if (vertexMergeFlag)
    for (int k = 0; k < 8; k++)
      if (isMergedCorner[k])
        neVertex.vertices.push_back((corner[k] << kTrisoupFpBits) - kTrisoupFpHalf);

  int vtxCount = (int)neVertex.vertices.size();
  Vec3<int32_t> gCenter = 0;
  for (int j = 0; j < vtxCount; j++)
    gCenter += neVertex.vertices[j];
  gCenter = gCenter * LUTdiv32[vtxCount] >> 12;
  ncVertex.gravityCenter = gCenter;
  ncVertex.pos = gCenter;

  if (vtxCount >= 3)
    neVertex.dominantAxis = findDominantAxis(neVertex.vertices, blockWidth, gCenter);

  if (leaf.isSkiped || vtxCount < 2)
    return;

  Vec3<int32_t> neiCentroid;
  bool backwardNei = false;
  bool isCentroidNeiPresent = false;
  if (vtxCount == 2 && isFaceVertexActivated && eligibleNonClosed(neVertex.vertices[0], neVertex.vertices[1], blockWidth))
    isCentroidNeiPresent = getNeighboringCentroid(
      neiCentroid, neVertex.dominantAxis, nodeIdx, neVertex.vertices[0], neVertex.vertices[1], blockWidth, backwardNei);

  // Skip leaves that have fewer than 3 vertices and not non closed surface
  bool backwardNonclosed = backwardNei && isCentroidNeiPresent;
  if (!backwardNonclosed && vtxCount == 2) {
    if (!backwardNei)
      isNonClosedForwardNode.back() = true;
    return;
  }

  // refine gravity center
  int scaleQ = 256; // scaled on 8 bits; 256 is one
  int quIndex = leaf.quIndex;
  int qpNode = quIndex != -1 ? ctxtMemOctree.listOfQUs[quIndex].localQP : gbh.trisoup_QP;
  int stepQcentro = LUT_QP_Trisoup[qpNode];
  if (vtxCount > 2)
    determineRefinedGravityCenter(gCenter, neVertex.vertices, blockWidth, scaleQ, stepQcentro);
  Vec3<int32_t> blockCentroid = gCenter;

  ncVertex.pos = blockCentroid;
  ncVertex.gravityCenter = gCenter;

  // Refinement of the centroid
  if (isCentroidResActivated && vtxCount >= 2) {
    // eligible doube centro
    uint16_t subCentroMask  = 0;
    bool eligDoubleCentro = false;
    if (vtxCount >= 6)
      eligibleDoubleCentro(
        gCenter, neVertex.dominantAxis, neVertex.vertices, blockWidth, subCentroMask,  eligDoubleCentro);

    // centroid residual
    bool isDoubleCentro = false;
    std::array<Vec3<int32_t>, 2> subCentroid;
    CentroidsResidualDecoder(stepQcentro, neVertex.dominantAxis, scaleQ, blockCentroid, vtxCount, leaf, neVertex.vertices,
      neiCentroid, eligDoubleCentro, isDoubleCentro, neVertex,subCentroMask, subCentroid);

    bool boundaryInside = checkIfPointIsInsideNode(blockCentroid, blockWidth << kTrisoupFpBits);

    // store centroid residual for next frame
    ncVertex ={ true, blockCentroid, isDoubleCentro == 1, subCentroid, subCentroMask, boundaryInside , gCenter };

    if (isFaceVertexActivated)
      generateFaceVerticesInNodeRasterScanDecoder(nodeIdx);
  }
}

//----------------------------------------------------------------------------
void
TrisoupDecoder
::FinalizeNonCloseNodesDecoder(
  const TrisoupNodeDecoder& leaf,
  const bool isCentroidResActivated,
  int nodeIdx)
{
  TrisoupNodeEdgeVertex& neVertex = eVerts[nodeIdx];
  //int vtxCount = (int)neVertex.vertices.size();  // equal to 2 for secodn pass
  int dominantAxis = neVertex.dominantAxis;

  int neiCentroIdx = neiNodeIdxVec[nodeIdx][dominantAxis + 3];
  if (neiCentroIdx == -1 || !cVerts[neiCentroIdx].valid)
    return;

  Vec3<int32_t> neiCentroid = cVerts[neiCentroIdx].pos;
  neiCentroid[dominantAxis] += blockWidth << kTrisoupFpBits;

  // get gravity center
  Vec3<int32_t> gCenter = cVerts[nodeIdx].gravityCenter;
  Vec3<int32_t> blockCentroid = gCenter;

  if (!isCentroidResActivated) { // isDoubleCentro is -1
    cVerts[nodeIdx] = { false, blockCentroid, false, {}, 0, true , gCenter };
  }
  else { // Refinement of the centroid along the dominant axis
    uint16_t subCentroMask  = 0;
    bool isDoubleCentro = false;

    int scaleQ = 256; // scaled on 8 bits; 256 is one

    int quIndex = leaf.quIndex;
    int qpNode = quIndex != -1 ? ctxtMemOctree.listOfQUs[quIndex].localQP : gbh.trisoup_QP;
    int stepQcentro = LUT_QP_Trisoup[qpNode];

    std::array<Vec3<int32_t>, 2> subCentroid;
    CentroidsResidualDecoder(stepQcentro, dominantAxis, scaleQ, blockCentroid, 2, leaf, neVertex.vertices,
      neiCentroid, false, isDoubleCentro, neVertex, subCentroMask, subCentroid);

    bool boundaryInside = checkIfPointIsInsideNode(blockCentroid, blockWidth << kTrisoupFpBits);

    // store centroid residual for next frame
    cVerts[nodeIdx] = { true, blockCentroid, false,  {}, 0, boundaryInside , gCenter };

    if (isFaceVertexActivated)
      generateFaceVerticesInNodeRasterScanDecoder(nodeIdx);
  }
}

//----------------------------------------------------------------------------

void
TrisoupDecoder
::generateFaceVerticesInNodeRasterScanDecoder(const int nodeIdx)
{
  int i = nodeIdx;
  // For 6-neighbour nodes of three which have smaller coordinates than current node,
  // if current node and its neighbour node have refined centroid vertex each other,
  // and when the centroids are connected, if there is no bite on the current surface,
  // then the intersection of centroid connection segment and
  // node boundary face is defined as the temporary face vertex.
  // And if original points are distributed around the temporary face vertex,
  // it is defined as a true face vertex and determine to connect these centroids,
  // and then face-flag becomes true.
  // to generate 3 faces per node, neighbour direction loop must be placed at the most outer loop.
  // x,y,z-axis order 0,1,2
  if (!cVerts[i].boundaryInside)
    return;

  for (int axis = 0; axis < 3; axis++) {
    int ii = neiNodeIdxVec[i][axis];
    // neighbour-node exists on this direction
    // centroid of the neighbour-node exists and inside of the boundary
    if (-1 == ii || !cVerts[ii].valid || !cVerts[ii].boundaryInside|| leaves[ii].isSkiped)
      continue;

    int cntVtxFace = 0; // count edge vertices included in the face
    for (int k = 0; k < eVerts[i].vertices.size(); k++)
      cntVtxFace += eVerts[i].vertices[k][axis] + kTrisoupFpHalf == 0;

    int eIdxCurr[2] = { -1, -1 };
    int eIdxNei[2] = { -1, -1 };
    Vec3<int32_t > fVert[2];

    if (1 == cntVtxFace && eVerts[i].vertices.size() == 2 && eVerts[ii].vertices.size() == 2) {
      findTrisoupFaceVertex(cVerts[i].pos, cVerts[ii].pos, axis, fVert);
      placeFaceVertex(eVerts[i], axis, fVert[0], eIdxCurr);
      placeFaceVertex(eVerts[ii], axis, fVert[1], eIdxNei);
      fVerts[i].EdgeVerIdxBeforeFaceVer.push_back(eIdxCurr[0]);
      fVerts[i].vertices.push_back(fVert[0]);
      fVerts[ii].EdgeVerIdxBeforeFaceVer.push_back(eIdxNei[0]);
      fVerts[ii].vertices.push_back(fVert[1]);
      continue;
    }
    if (2 != cntVtxFace && 3 != cntVtxFace)
      continue;

    findTrisoupFaceVertex(cVerts[i].pos, cVerts[ii].pos, axis, fVert);
    placeFaceVertex(eVerts[i], axis, fVert[0], eIdxCurr);
    placeFaceVertex(eVerts[ii], axis, fVert[1], eIdxNei);

    if (-1 == eIdxCurr[0] || -1 == eIdxCurr[1])
      continue;

    // c0, c1, and face vertex is on the same side of the surface
    bool valid = validateFaceVertex(eVerts[i], cVerts[i], cVerts[ii], eIdxCurr[0], eIdxCurr[1], fVert[0]);
    if (valid) {
      bool faceConnected = arithmeticDecoder->decode(ctxFaces);
      if (faceConnected) {
        fVerts[i].EdgeVerIdxBeforeFaceVer.push_back(eIdxCurr[0]);
        fVerts[i].vertices.push_back(fVert[0]);
        fVerts[ii].EdgeVerIdxBeforeFaceVer.push_back(eIdxNei[0]);
        fVerts[ii].vertices.push_back(fVert[1]);
      }
    } // end valid
  } // end loop on axis
}

//----------------------------------------------------------------------------

void
TrisoupDecoder
::initDecoder()
{
  BBorig = gbh.geomBoxOrigin;
  keyshift = BBorig + (1 << 18);

  lastWedgex = currWedgePos[0];

  if (useLocalAttr) {
    localPointCloud.reserve(pointCloud.getPointCount());
    numSlabBlocksPerDim = {
      (gbh.maxRootNodeSize + slabBlockSize[0] - 1) / slabBlockSize[0],
      (gbh.maxRootNodeSize + slabBlockSize[1] - 1) / slabBlockSize[1],
      (gbh.maxRootNodeSize + slabBlockSize[2] - 1) / slabBlockSize[2]
    };
    localSlabBlockStart = 0;
    localSlabBlockIdxStart = 0;
    // we only need to buffer slab blocks for two successive slabs
    // when one is finished it can be be processed, the next one starts with
    // partial trisoup nodes' content
    int numSlabBlocksPerSlab = numSlabBlocksPerDim[1] * numSlabBlocksPerDim[2];
    laPoints[0].resize(numSlabBlocksPerSlab);
    laPoints[1].resize(numSlabBlocksPerSlab);
  }

  if (useLocalAttr) {
    localPointCloud.addRemoveAttributes(pointCloud);
  }

  // determine vertex quantization step
  midBlock = (blockWidth - 1) << 7; // << 8) >> 1;
  maxBlock = (blockWidth - 1) << 8;

  // quantization mapping par half segment
  for (int QP = 0; QP < 55; QP++) {
    int stepInv = LUT_QP_Trisoup_inv[QP];
    LUT_maxMagnitude[QP] = midBlock * stepInv >> 21;

    // max number of bits of quantization
    int nBitMag = 0;
    while ((1 << nBitMag) < LUT_maxMagnitude[QP] + 1)
      nBitMag++;
    LUT_nBitMagnitude[QP] = nBitMag;
  }

  // rendering
  renderedBlock.resize(blockWidth * blockWidth * 16, 0);

  haloFlag = gbh.trisoup_halo_flag;
  vertexMergeFlag = gbh.trisoup_vertex_merge_flag;
  vertexConsistencyFlag = gbh.trisoup_vertex_consistency_flag;
  thickness = gbh.trisoup_thickness;
  isCentroidResActivated = gbh.trisoup_centroid_vertex_residual_flag;
  isFaceVertexActivated = gbh.trisoup_face_vertex_flag;

  if (haloFlag && blockWidth > 4) {
    for (int QP = 0; QP < 55; QP++) {
      int step = LUT_QP_Trisoup[QP];
      int haloTriangle = std::max(0, step - 256) * LUTdiv32[blockWidth] >> 12;
      haloTriangle = (haloTriangle * 28) >> 5;
      haloTriangle = haloTriangle > 36 ? 36 : haloTriangle;
      haloTriangleQP[QP] = haloTriangle;
    }
  }
}

//----------------------------------------------------------------------------

void
TrisoupDecoder
::sliceDecoder(bool isFinalPass)
{
  neighbNodes.reserve(leaves.size() * 12); // at most 12 edges per node (to avoid reallocations)
  segmentUniqueIndex.resize(12 * leaves.size(), -1); // temporarily set to -1 to check everything is working
  // TODO: set to -1 could be removed when everything will work properly

  if (isFaceVertexActivated) {
    fVerts.resize(leaves.size());
  }

  if (currWedgePos == Vec3<int32_t>{INT32_MIN, INT32_MIN, INT32_MIN}) {
    currWedgePos = leaves[0].pos;
  }

  int lastAcceptable = INT_MAX;
  if (!isFinalPass)
    lastAcceptable = leaves.back().pos[0];

  while (nextIsAvailable() && currWedgePos[0] < lastAcceptable) { // this a loop on start position of edges; 3 edges along x,y and z per start position
    // process current wedge position

    std::array<bool, 8> isNeigbourSane;
    for (int i = 0; i < 8; ++i) // sanity of neighbouring nodes
      isNeigbourSane[i] = edgesNeighNodes[i] < leaves.size() && currWedgePos + offsets[i] == leaves[edgesNeighNodes[i]].pos;

    if (isFaceVertexActivated && isNeigbourSane[6]) {
      std::array<int, 6> nei3idx = { -1, -1, -1, -1, -1, -1 };
      if (isNeigbourSane[5]) {
        nei3idx[2] = edgesNeighNodes[5];
        neiNodeIdxVec[nei3idx[2]][5] = neiNodeIdxVec.size();
      }
      if (isNeigbourSane[4]) {
        nei3idx[1] = edgesNeighNodes[4];
        neiNodeIdxVec[nei3idx[1]][4] = neiNodeIdxVec.size();
      }
      if (isNeigbourSane[2]) {
        nei3idx[0] = edgesNeighNodes[2];
        neiNodeIdxVec[nei3idx[0]][3] = neiNodeIdxVec.size();
      }
      neiNodeIdxVec.push_back(nei3idx);
    }

    for (int dir = 0; dir < 3; ++dir) { // this the loop on the 3 edges along z, then y, then x
       // verify there is an edge
      bool processedEdge = false;
      for (int neighIdx = 0; neighIdx < 4; ++neighIdx) {
        int edgeNeighNodeIdx = edgesNeighNodesIdx[dir][neighIdx]; // gives the neighbour index in the list of 8 neighbours
        processedEdge = processedEdge || isNeigbourSane[edgeNeighNodeIdx]; // at least one neighbor intersects edge
      }
      if (!processedEdge)
        continue;

      int dir0 = axisdirection[dir][0];
      int dir1 = axisdirection[dir][1];
      int dir2 = axisdirection[dir][2];

      // for TriSoup Vertex inter prediction
      int countNearPointsPred = 0;
      int distanceSumPred = 0;

      bool anyNodeSkip = false;
      bool noGenVertex = true;
      int nNodeNotSkip = 0;

      if (enableSkip) {
        for (int neighIdx = 0; neighIdx < 4; ++neighIdx) { // this the loop on the 4 nodes to interesect the edge
          int edgeNeighNodeIdx = edgesNeighNodesIdx[dir][neighIdx]; // gives the neighbour index in the list of 8 neighbours
          if (isNeigbourSane[edgeNeighNodeIdx]) { // test for sanity of the neighbour, process only if sane
            int neighbNodeIndex = edgesNeighNodes[edgeNeighNodeIdx];
            bool nodeSkipStatus = leaves[neighbNodeIndex].isSkiped;
            nNodeNotSkip += !nodeSkipStatus;
            anyNodeSkip = anyNodeSkip || nodeSkipStatus;
            noGenVertex = noGenVertex && nodeSkipStatus; // If all neighboring (shared-edge) nodes are skipped, do not generate vertices.
          }
        }
      }
      anyNodeSkip = anyNodeSkip && (nNodeNotSkip <= 2);

      // local QU
      int countNodesQU = 0;
      int QPedge = 0;

      uint16_t neighboursMask = 0;
      std::array<int, 18> pattern{ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
      for (int neighIdx = 0; neighIdx < 4; ++neighIdx) { // this the loop on the 4 nodes to interesect the edge
        int edgeNeighNodeIdx = edgesNeighNodesIdx[dir][neighIdx]; // gives the neighbour inde in the list of 8 neighbours

        if (isNeigbourSane[edgeNeighNodeIdx]) { // test for sanity of the neighbour, process only if sane
          int neighbNodeIndex = edgesNeighNodes[edgeNeighNodeIdx];
          int idx = neighbNodeIndex * 12 + edgeIdx[dir][neighIdx];
          segmentUniqueIndex[idx] = uniqueIndex;

          // update mask from nodes touching edge
          neighboursMask |= neighMask[dir][neighIdx];

          int indexLow = edgeIdx[dir][neighIdx];
          for (int v = 0; v < 11; v++) {
            if (localEdgeindex[indexLow][v] == -1)
              break;

            int indexV = neighbNodeIndex * 12 + localEdgeindex[indexLow][v]; // index of segment
            int Vidx = segmentUniqueIndex[indexV];
            //assert(Vidx != -1); // check if already coded
            pattern[patternIndex[indexLow][v]] = Vidx;
          }

          // determine TriSoup Vertex
          auto& leafNei = leaves[neighbNodeIndex];
          const int offset1 = leafNei.pos[dir1] < currWedgePos[dir1];
          const int offset2 = leafNei.pos[dir2] < currWedgePos[dir2];
          const int pos1 = currWedgePos[dir1] - offset1;
          const int pos2 = currWedgePos[dir2] - offset2;

          // determine TriSoup Vertex inter prediction
          if (isInter && !(enableSkip && noGenVertex)) {
            for (int j = leafNei.predStart; j < leafNei.predEnd; j++) {
              Vec3<int> voxel = compensatedPointCloud[j];
              int pointOK = voxel[dir1] == pos1 && voxel[dir2] == pos2;
              countNearPointsPred += pointOK;
              distanceSumPred += pointOK * voxel[dir0];
            }
          }

          // determine QP from QU
          countNodesQU++;
          QPedge += leafNei.quIndex != -1
            ? ctxtMemOctree.listOfQUs[leafNei.quIndex].localQP : gbh.trisoup_QP;
        }
      } // end loop on 4 neighbouring nodes
      distanceSumPred -= countNearPointsPred * currWedgePos[dir0];

      int segmentUniqueIdxPrevEdge = -1;
      for (int prevNeighIdx = 0; prevNeighIdx < 4; ++prevNeighIdx) {
        int wedgeNeighNodeIdx = wedgeNeighNodesIdx[dir][prevNeighIdx];
        if (isNeigbourSane[wedgeNeighNodeIdx]) {
          // update current mask from nodes touching wedge
          neighboursMask |= wedgeNeighMask[dir][prevNeighIdx];
          if (segmentUniqueIdxPrevEdge == -1) {
            int idx = edgesNeighNodes[wedgeNeighNodeIdx] * 12 + wedgeNeighNodesEdgeIdx[dir][prevNeighIdx];
            segmentUniqueIdxPrevEdge = segmentUniqueIndex[idx];
            pattern[0] = segmentUniqueIdxPrevEdge;
            //assert(segmentUniqueIdxPrevEdge != -1);
            for (int neighIdx = 0; neighIdx < 4; ++neighIdx) {
              int edgeNeighNodeIdx = edgesNeighNodesIdx[dir][neighIdx];
              if (isNeigbourSane[edgeNeighNodeIdx]) {
                neighbNodes[segmentUniqueIdxPrevEdge] |= toPrevEdgeNeighMask[dir][neighIdx];
              }
            }
          }
        }
      }

      ++uniqueIndex;
      neighbNodes.push_back(neighboursMask);
      edgePattern.push(pattern);
      posForEdgeOfVertex.push(currWedgePos);

      // for colocated edges
      int64_t key = dir
        + (int64_t(currWedgePos[0] + keyshift[0]) << 42)
        + (int64_t(currWedgePos[1] + keyshift[1]) << 22)
        + (int64_t(currWedgePos[2] + keyshift[2]) << 2);
      currentFrameEdgeKeys.push_back(key);

      // QP for edge
      QPedge = QPedge * LUTdiv32[countNodesQU] + 2048 >> 12;

      // determine TriSoup Vertex inter prediction and skip
      int QPedgeWithSkip = anyNodeSkip ? 0 : QPedge;
      TriSoupEdgeLocalQP.push_back(QPedgeWithSkip);
      int8_t vertexQ = 0;
      if (isInter) {
        bool isSkippedVertex = false;
        int8_t vertexPos = 0;
        if (countNearPointsPred > 0) {
          vertexPos = quantizeQP(distanceSumPred, countNearPointsPred, QPedgeWithSkip);
        }
        TriSoupVerticesPred.push(vertexPos);

        if (anyNodeSkip) {
          isSkippedVertex = true;
          vertexQ = countNearPointsPred > thVertexDeterminationSkip ? vertexPos : 0;
        }
        TriSoupVerticesIsSkipped.push(isSkippedVertex); //skip status
      }
      TriSoupVerticesQP.push_back(vertexQ);

    } // end loop on three directions

    // move to next wedge
    goNextWedge(isNeigbourSane);


    // code vertices and rendering of preceding slices in case the loop has moved up one slice or if finished
    if (changeSlice()) {

      // clean merge map
      if (vertexMergeFlag && nodeIdxC < leaves.size()) {
        Vec3< int64_t> posKey = leaves[nodeIdxC].pos;
        const int64_t shift = 1 << 18;
        int64_t keyEdge = (posKey[0] + shift << 40) + (posKey[1] + shift << 20) + posKey[2] + shift;
        auto it = triNodeEnable.lower_bound(keyEdge);
        triNodeEnable.erase(triNodeEnable.begin(), it);
      }

      // coding vertices
      int upperxForCoding = !nextIsAvailable() ? INT32_MAX : currWedgePos[0] - blockWidth;
      while (!posForEdgeOfVertex.empty() && posForEdgeOfVertex.front()[0] < upperxForCoding) {
        // spatial neighbour and inter comp predictors
        int8_t  interPredictor = isInter ? TriSoupVerticesPred.front() : 0;
        auto pattern = edgePattern.front();

        // colocated edge predictor
        int8_t colocatedVertex = 0;
        int8_t colocatedVertexQP = gbh.trisoup_QP;
        auto keyCurrent = currentFrameEdgeKeys[firstVertexToCode];
        if (interSkipEnabled) {
          while (colocatedEdgeIdx < ctxtMemOctree.refFrameEdgeKeys.size() - 1
              && ctxtMemOctree.refFrameEdgeKeys[colocatedEdgeIdx] < keyCurrent)
            colocatedEdgeIdx++;

          if (ctxtMemOctree.refFrameEdgeKeys[colocatedEdgeIdx] == keyCurrent) {
            colocatedVertex = ctxtMemOctree.refFrameEdgeValue[colocatedEdgeIdx];
            colocatedVertexQP = ctxtMemOctree.refFrameEdgeQP[colocatedEdgeIdx];
          }
        }

        // code edge
        int8_t vertex;
        int QPedge = TriSoupEdgeLocalQP[firstVertexToCode];
        if (colocatedVertex && QPedge != colocatedVertexQP)
          // requantize colocated edge according to local QP
          colocatedVertex = quantizeQP(
            dequantizeQP(colocatedVertex, colocatedVertexQP), QPedge);

        // decode vertex
        if (!(isInter && TriSoupVerticesIsSkipped.front()))
          decodeOneTriSoupVertexRasterScan(
            arithmeticDecoder, ctxtMemOctree, TriSoupVerticesQP,TriSoupVertices2bits,
            neighbNodes[firstVertexToCode], pattern, interPredictor,
            colocatedVertex, firstVertexToCode, QPedge);

        vertex = TriSoupVerticesQP[firstVertexToCode];
        TriSoupVertices2bits.push_back(map2bits(vertex, QPedge));

        if (vertexMergeFlag && vertex) {
          const int32_t thMerge = std::min(8 << kTrisoupFpBits, blockWidth << kTrisoupFpBits - 2);
          Vec3< int64_t> posKey = posForEdgeOfVertex.front();
          int32_t deQVert = dequantizeQP(vertex, QPedge) + kTrisoupFpHalf;

          if (deQVert < thMerge) {
            const int64_t shift = 1 << 18;
            int64_t keyEdge = (posKey[0] + shift << 40) + (posKey[1] + shift << 20) + posKey[2] + shift;
            triNodeEnable[keyEdge]++;
          }
          if ((blockWidth << kTrisoupFpBits) - deQVert < thMerge) {
            posKey[2 - (keyCurrent & 3)] += blockWidth;
            const int64_t shift = 1 << 18;
            int64_t keyEdge = (posKey[0] + shift << 40) + (posKey[1] + shift << 20) + posKey[2] + shift;
            triNodeEnable[keyEdge]++;
          }
        }
        posForEdgeOfVertex.pop();
        edgePattern.pop();
        if (isInter) {
          TriSoupVerticesPred.pop();
          TriSoupVerticesIsSkipped.pop();
        }
        firstVertexToCode++;
      }

      // centroid processing
      int upperxForVertex = !nextIsAvailable() ? INT32_MAX : currWedgePos[0] - 2 * blockWidth;
      while (nodeIdxC < leaves.size() && leaves[nodeIdxC].pos[0] < upperxForVertex) {
        auto leaf = leaves[nodeIdxC];

        TrisoupOriginalPCinfo oriPCinfo;
        if (vertexConsistencyFlag)
          prepareVertexConsistencyDecoder(arithmeticDecoder, oriPCinfo);

        GenerateCentroidsInNodeRasterScanDecoder(   // dequantization of vertices is here
          leaf, TriSoupVerticesQP, TriSoupEdgeLocalQP, idxSegment,
          isCentroidResActivated, segmentUniqueIndex, &oriPCinfo, nodeIdxC);

        nodeIdxC++;
      }

      // rendering by TriSoup triangles + local attributes rendering
      int upperxForRendering = !nextIsAvailable() ?  INT32_MAX  : currWedgePos[0] - 3 * blockWidth;
      int nonClosedIdx = firstNodeToRender;
      bool finishedNonClose = false;
      while (firstNodeToRender < leaves.size() && leaves[firstNodeToRender].pos[0] < upperxForRendering) {

        if(isFaceVertexActivated && nonClosedIdx < leaves.size() && leaves[nonClosedIdx].pos[0] < upperxForRendering) {
          if (isNonClosedForwardNode[nonClosedIdx])
            FinalizeNonCloseNodesDecoder(leaves[nonClosedIdx], isCentroidResActivated, nonClosedIdx);
          nonClosedIdx++;
        } else
          finishedNonClose = true;

        bool eligibleRendering = finishedNonClose || nonClosedIdx < leaves.size() && (leaves[nonClosedIdx].pos[0] > leaves[firstNodeToRender].pos[0]
          || leaves[nonClosedIdx].pos[1] - leaves[firstNodeToRender].pos[1] >= 2 * blockWidth);
        if (eligibleRendering) {
          auto leaf = leaves[firstNodeToRender];

          // selection of the right Slab Block
          // and rendering attributes of any finished Slab
          int slabBlockIdx;
          if (useLocalAttr) {
            int nodeposX = leaves[firstNodeToRender].pos[0];
            int nodeposY = leaves[firstNodeToRender].pos[1];
            int nodeposZ = leaves[firstNodeToRender].pos[2];

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
                      slabBlockIdx = localSlabBlockIdxStart[0];
                      for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
                        for (int z = 0; z < numSlabBlocksPerDim[2]; ++z, ++slabBlockIdx) {
                          int numPoints = laPoints[slabIdxMod2][slabBlockIdx].size();
                          if (numPoints) {
                            int size = localPointCloud.size();
                            localPointCloud.resize(size + numPoints);
                            for (int i = 0; i < numPoints; ++i)
                              localPointCloud[size++] = laPoints[slabIdxMod2][slabBlockIdx][i];
                            laPoints[slabIdxMod2][slabBlockIdx].clear(); // only clear, may be reused
                            localSlabBlockStart[2] = z * slabBlockSize[2];
                            localSlabBlockStart[1] = y * slabBlockSize[1];
                            processLocalAttributesDecoder(pointCloud, false);
                            localPointCloud.clear();
                          }
                        }
                      // next slice
                      slabIdxMod2 = slabIdxMod2 ^ 1;
                      localSlabBlockIdxStart[0] = 0; //numSlabBlocksPerDim[1] * numSlabBlocksPerDim[2]; // = localSlabBlockIdxStart[2];
                      localSlabBlockStart[0] += slabBlockSize[0];
                      localSlabBlockStart[2] = localSlabBlockStart[1] = 0;
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

          int nPointsInBlock = 0;
          if (leaf.isSkiped)
            nPointsInBlock = renderSkip(renderedBlock, nPointsInBlock, leaf);
          else
            nPointsInBlock = generateTrianglesInNodeRasterScan(
              leaf, firstNodeToRender, renderedBlock, thickness, isFaceVertexActivated);

          if (nPointsInBlock) {
            if (useLocalAttr) {
              flush2SlabBuff(renderedBlock.begin(), nPointsInBlock,
                slabIdxMod2, slabBlockIdx, leaf.isSkiped);
            } else {
              flush2PointCloud(
                nRecPoints, renderedBlock.begin(), nPointsInBlock,
                pointCloud, leaf.isSkiped);
            }
          }
          firstNodeToRender++;
        }
      }
    } // end if on slice chnage

    lastWedgex = currWedgePos[0];
  } // end while loop on wedges
}

//----------------------------------------------------------------------------

void
TrisoupDecoder
::finishSliceDecoder()
{
  // store edges for colocated next frame
  ctxtMemOctree.refFrameEdgeKeys = currentFrameEdgeKeys;
  ctxtMemOctree.refFrameEdgeValue = TriSoupVerticesQP;
  ctxtMemOctree.refFrameEdgeQP = TriSoupEdgeLocalQP;

  if (useLocalAttr) {
    ////
    // The following is to render the remaining Slabs
    int numSlabs = 1 + (localSlabBlockStart[0] + slabBlockSize[0] < leaves.back().pos[0] + blockWidth);
    for (int x = 0; x < numSlabs; ++x) {
      int slabBlockIdx = localSlabBlockIdxStart[0];
      for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
        for (int z = 0; z < numSlabBlocksPerDim[2]; ++z, ++slabBlockIdx) {
          int numPoints = laPoints[slabIdxMod2][slabBlockIdx].size();
          if (numPoints) {
            int size = localPointCloud.size();
            localPointCloud.resize(size + numPoints);
            for (int i = 0; i < numPoints; ++i)
              localPointCloud[size++] = laPoints[slabIdxMod2][slabBlockIdx][i];
            laPoints[slabIdxMod2][slabBlockIdx] = std::vector<point_t>(); // just release memory
            localSlabBlockStart[2] = z * slabBlockSize[2];
            localSlabBlockStart[1] = y * slabBlockSize[1];
            processLocalAttributesDecoder(pointCloud, true);
            localPointCloud.clear();
          }
        }
      // next slice/slab
      slabIdxMod2 = slabIdxMod2 ^ 1;
      localSlabBlockIdxStart[0] = 0; //numSlabBlocksPerDim[1] * numSlabBlocksPerDim[2]; // = localSlabBlockIdxStart[2];
      localSlabBlockStart[0] += slabBlockSize[0];
      localSlabBlockStart[2] = localSlabBlockStart[1] = 0;
    }
  }

  assert(nRecPoints == pointCloud.getPointCount());

  ////
  // The following is to render single slab-block for full slice
  if (!useLocalAttr) {
    processGlobalAttributesDecoder(pointCloud);
  }
  clearTrisoupElements();
}

//----------------------------------------------------------------------------

void
TrisoupDecoder::processLocalAttributesDecoder(
  PCCPointSet3& recPointCloud, bool isLast)
{
  //auto allocatedSizeLocal = localPointCloud.size();
  auto nRecPointsLocal = localPointCloud.size();
  decoder->processNextSlabBlockAttributes(localPointCloud, localSlabBlockStart);

  if (recPointCloud.getPointCount() < nRecPoints + nRecPointsLocal)
    recPointCloud.resize(nRecPoints + nRecPointsLocal + PC_PREALLOCATION_SIZE);

  recPointCloud.setFromPartition(localPointCloud, 0, nRecPointsLocal, nRecPoints);
  nRecPoints += nRecPointsLocal;
}

//----------------------------------------------------------------------------

void
TrisoupDecoder::processGlobalAttributesDecoder(
  PCCPointSet3& recPointCloud)
{
  decoder->processNextSlabBlockAttributes(pointCloud, {0, 0, 0});
}

//============================================================================

}  // namespace pcc