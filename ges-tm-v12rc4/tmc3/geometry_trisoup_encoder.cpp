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

#include "geometry_trisoup_encoder.h"

#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"
#include "PCCTMC3Encoder.h"

namespace pcc {

//============================================================================
void
encodeGeometryTrisoup(
  const EncoderParams& encParams,
  const GeometryParameterSet& gps,
  GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  GeometryOctreeContexts& ctxtMemOctree,
  MotionEntropy& ctxtMemMotion,
  std::vector<std::unique_ptr<EntropyEncoder>>& arithmeticEncoders,
  const CloudFrame& refFrame,
  const SequenceParameterSet& sps,
  InterPredParams& interPredParams,
  PCCTMC3Encoder3& encoder)
{
  bool isInter = gbh.slice_inter_prediction_flag;

  // prepare TriSoup parameters
  int blockWidth = gbh.trisoupNodeSize(gps);

  std::cout << "TriSoup QP = " << gbh.trisoup_QP << "\n";

  // get first encoder
  pcc::EntropyEncoder* arithmeticEncoder = arithmeticEncoders.begin()->get();

  point_t slabBlockSize = {
    sps.localized_attributes_slab_thickness_minus1 + 1,
    sps.localized_attributes_slab_block_size_minus1 + 1,
    sps.localized_attributes_slab_block_size_minus1 + 1
  };

  // used to determine if full frame optimized implementation is used or not
  {
    // Derived parameter used by trisoup
    gbh.maxRootNodeSize = gbh.trisoupNodeSize(gps) << gbh.rootNodeSizeLog2 - gbh.trisoupNodeSizeLog2(gps); // get bounding box of point cloud geometry
    gbh.uniqueSlabBlock = gbh.uniqueSlabBlock && gbh.maxRootNodeSize <= std::min(
      sps.localized_attributes_slab_thickness_minus1 + 1, sps.localized_attributes_slab_block_size_minus1 + 1);
  }

  // trisoup uses octree coding until reaching the triangulation level.
  TrisoupEncoder trisoup(blockWidth, pointCloud,
    1, isInter, interPredParams.compensatedPointCloud,
    gps, gbh, arithmeticEncoder, ctxtMemOctree,
    !gbh.uniqueSlabBlock, slabBlockSize, &encoder);

  trisoup.initEncoder();
  trisoup.thVertexDeterminationEnc = encParams.trisoup.thVertexDetermination;
  trisoup.thVertexDeterminationSkip = gbh.trisoup_early_skip_vertex_determination_threshold;

  // octree
  encodeGeometryOctree<true>(
    encParams, gps, gbh, pointCloud, ctxtMemOctree, ctxtMemMotion, arithmeticEncoders,
    nullptr, refFrame, sps, interPredParams, encoder, &trisoup);

  std::cout << "Size compensatedPointCloud for TriSoup = "
    << interPredParams.compensatedPointCloud.getPointCount() << "\n";
  std::cout << "Number of nodes for TriSoup = " << trisoup.leaves.size() << "\n";
  std::cout << "TriSoup gives " << pointCloud.getPointCount() << " points \n";

}

//============================================================================

int
determineCentroidResidual(
  Vec3<int32_t> normalV,
  Vec3<int32_t> blockCentroid,
  Vec3<int32_t> nodepos,
  PCCPointSet3& pointCloud,
  int start,
  int end,
  int lowBound,
  int  highBound,
  int stepQcentro,
  int  scaleQ,
  int& res,
  int vtxCount,
  int blockWidth,
  bool eligDoubleCentro,
  int dominantAxis,
  bool& isDoubleCentro)
{
  // determine quantized residual
  int counter = 0;
  int maxD = std::max(1, stepQcentro >> 9);
  bool posNonclosed[64] = { false };

  // determine residual
  int numPointsAroundCentroid = 0;
  for (int p = start; p < end; p++) {
    auto point = (pointCloud[p] - nodepos) << kTrisoupFpBits;

    Vec3<int64_t> CP = crossProduct(normalV, point - blockCentroid) >> kTrisoupFpBits;
    int64_t dist = isqrt(CP[0] * CP[0] + CP[1] * CP[1] + CP[2] * CP[2]);
    dist >>= kTrisoupFpBits;

    if ((dist << 10) <= 1774 * maxD) {
      if (vtxCount != 2) {
        int32_t w = (1 << 10) + 4 * (1774 * maxD - ((1 << 10) * dist));
        counter += w >> 10;
        res += (w >> 10) * ((normalV * (point - blockCentroid)) >> kTrisoupFpBits);
      } else {
        int dists = (normalV * (point - blockCentroid)) >> kTrisoupFpBits;
        posNonclosed[32 + (dists >> kTrisoupFpBits)] = true;
      }
    }
    if (eligDoubleCentro) {
      int distance = point[dominantAxis] > blockCentroid[dominantAxis] ? point[dominantAxis] - blockCentroid[dominantAxis] : blockCentroid[dominantAxis] - point[dominantAxis];
      if (distance <= (1 << kTrisoupFpBits))
        numPointsAroundCentroid++;
    }
  }

  if (vtxCount == 2) {
    res = divApprox((int64_t(findOffsetForNonClosedSurface(posNonclosed)) << 8), scaleQ, 0);
    res = res < 0 ? 0 : res;
  } else if (counter) {
    res = divApprox((int64_t(res) << 8), counter * scaleQ, 0);// res is shift by 8 bit // res is shift by +8 due to scaleQ
  }

  int half = stepQcentro >> 1;
  int DZ = 682 * half >> 10; //2 * half / 3;

  int resCentroQ = 0;
  if (abs(res) >= DZ) {
    resCentroQ = 1 + (abs(res) ) / stepQcentro;
    if (res < 0)
      resCentroQ = -resCentroQ;
  }
  resCentroQ = std::min(std::max(resCentroQ, -lowBound), highBound);  // res in [-lowBound; highBound]

  isDoubleCentro = eligDoubleCentro ? (numPointsAroundCentroid > 1 ? false : true) : false;
  return resCentroQ;
}

// --------------------------------------------------------------------------
//  encoding of centroid residual in TriSOup node
void
encodeCentroidResidual(
  int resCentroQ,
  pcc::EntropyEncoder* arithmeticEncoder,
  GeometryOctreeContexts & ctxtMemOctree,
  int resCentroQPred,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound)
{
  if (resCentroQPred == -100) //intra
     arithmeticEncoder->encode(resCentroQ == 0, ctxtMemOctree.ctxRes0[ctxMinMax][0]);
  else { //inter
    arithmeticEncoder->encode(resCentroQ == 0, ctxtMemOctree.ctxRes0[ctxMinMax][1 + std::min(3, std::abs(resCentroQPred))]);
  }

  // if not 0, residual in [-lowBound; highBound]
  if (resCentroQ) {
    // code sign
    if (highBound && lowBound) {  // otherwise sign is known
      int lowS = std::min(7, lowBoundSurface);
      int highS = std::min(7, highBoundSurface);
      arithmeticEncoder->encode(resCentroQ > 0, ctxtMemOctree.ctxResSign[lowBound == highBound ? 0 : 1 + (lowBound < highBound)][lowS][highS][(resCentroQPred && resCentroQPred != -100) ? 1 + (resCentroQPred > 0) : 0]);
    }

    // code remaining bits 1 to 7 at most
    int magBound = (resCentroQ > 0 ? highBound : lowBound) - 1;
    bool sameSignPred = resCentroQPred != -100 && (resCentroQPred > 0 && resCentroQ > 0) || (resCentroQPred < 0 && resCentroQ < 0);

    int magRes = std::abs(resCentroQ) - 1;
    int ctx = 0;
    int ctx2 = resCentroQPred != -100 ? 1 + std::min(8, sameSignPred * std::abs(resCentroQPred)) : 0;
    while (magBound > 0 && magRes >= 0) {
      if (ctx < 4)
        arithmeticEncoder->encode(magRes == 0, ctxtMemOctree.ctxResMag[ctx][ctx2]);
      else
        arithmeticEncoder->encode(magRes == 0);

      magRes--;
      magBound--;
      ctx++;
    }
  }  // end if not 0
}

//============================================================================

void
TrisoupEncoder
::encodeOneTriSoupVertexRasterScan(
  int8_t vertexQP,
  pcc::EntropyEncoder* arithmeticEncoder,
  GeometryOctreeContexts& ctxtMemOctree,
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

  // encode vertex presence
  int ctxMap1, ctxMap2, ctxMap2Base, ctxInter;
  int shiftMag = std::max(0, nBitMag - 1);
  int interPredictor2bit = 0;
  if (nBitMag >= 1)
  {
    int bit1 = interPredictor > 0;
    int bit2 = ((std::abs(interPredictor) - 1) >> nBitMag - 1) & 1;
    interPredictor2bit = (bit1 << 1) + (bit1 ? bit2 : !bit2);
  }

  constructCtxPresence(ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, colocatedVertex);
  int ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[0].getEvolve(
    vertexQP != 0, ctxMap2, ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
    ctxtMemOctree._BufferOBUFleavesTrisoup);

  arithmeticEncoder->encode(
    (int)(vertexQP != 0), ctxTrisoup >> 3,
    ctxtMemOctree.ctxTriSoup[0][ctxInter][ctxTrisoup],
    ctxtMemOctree.ctxTriSoup[0][ctxInter].obufSingleBound);

  int magnitude = std::abs(vertexQP) - 1;

  // encode  vertex position
  if (vertexQP != 0) {
    int partialDec = 0;
    int bitRemaining = nBitMag - 1;
    int partialMag = 0;
    int blockWidthLog2 = ilog2(uint32_t(blockWidth));

    // first position bit is for left vs right  -> sign of vertexQP
    constructCtxPos1(
      ctxMap1, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor,
      colocatedVertex, interPredictor2bit);
    int bit = vertexQP > 0;

    ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[1].getEvolve(
      bit, ctxMap2, ctxMap1, &ctxtMemOctree._OBUFleafNumberTrisoup,
      ctxtMemOctree._BufferOBUFleavesTrisoup);
    arithmeticEncoder->encode(
      bit, ctxTrisoup >> 3,
      ctxtMemOctree.ctxTriSoup[1][ctxInter][ctxTrisoup],
      ctxtMemOctree.ctxTriSoup[1][ctxInter].obufSingleBound);
    partialDec = bit;

    // second position bit
    int goodColo = 0;
    if (bitRemaining >= 0) {
      constructCtxPos2(
        ctxMap2Base, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor,
        shiftMag, partialDec, colocatedVertex, blockWidthLog2, goodColo);
      bit = (magnitude >> bitRemaining) & 1;

      if (isCodedBit(partialMag, bitRemaining, maxMag)) {
        ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[2].getEvolve(
          bit, ctxMap2, (ctxMap1 << 1) + partialDec,
          &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);
        arithmeticEncoder->encode(
          bit, ctxTrisoup >> 3,
          ctxtMemOctree.ctxTriSoup[2][ctxInter][ctxTrisoup],
          ctxtMemOctree.ctxTriSoup[2][ctxInter].obufSingleBound);
      }
      partialDec = (partialDec << 1) | bit;
      partialMag = (partialMag << 1) | bit;
      bitRemaining--;
    }

    // third bit
    if (bitRemaining >= 0 ) {
      constructCtxPos3(ctxMap2Base, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, shiftMag, partialDec, colocatedVertex, goodColo);
      bit = (magnitude >> bitRemaining) & 1;

      if (isCodedBit(partialMag, bitRemaining, maxMag)) {
        ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[3].getEvolve(
          bit, ctxMap2, (ctxMap1 << 2) + partialDec,
          &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);
        arithmeticEncoder->encode(
          bit, ctxTrisoup >> 3,
          ctxtMemOctree.ctxTriSoup[3][ctxInter][ctxTrisoup],
          ctxtMemOctree.ctxTriSoup[3][ctxInter].obufSingleBound);
      }
      partialDec = (partialDec << 1) | bit;
      partialMag = (partialMag << 1) | bit;
      bitRemaining--;
    }

    // fourth bit
    if (bitRemaining >= 0) {
      constructCtxPos4(ctxMap2Base, ctxMap2, ctxInter, ctxInfo, isInter, interPredictor, shiftMag, partialDec, colocatedVertex, goodColo);
      bit = (magnitude >> bitRemaining) & 1;

      if (isCodedBit(partialMag, bitRemaining, maxMag)) {
        ctxTrisoup = ctxtMemOctree.MapOBUFTriSoup[4].getEvolve(
          bit, ctxMap2, (ctxMap1 << 3) + partialDec,
          &ctxtMemOctree._OBUFleafNumberTrisoup,
          ctxtMemOctree._BufferOBUFleavesTrisoup);
        arithmeticEncoder->encode(
          bit, ctxTrisoup >> 3,
          ctxtMemOctree.ctxTriSoup[4][ctxInter][ctxTrisoup],
          ctxtMemOctree.ctxTriSoup[4][ctxInter].obufSingleBound);
      }
      partialMag = (partialMag << 1) | bit;
      bitRemaining--;
    }

    // remaining bits are bypassed
    for (; bitRemaining >= 0; bitRemaining--) {
      bit = (magnitude >> bitRemaining) & 1;
      if (isCodedBit(partialMag, bitRemaining, maxMag))
        arithmeticEncoder->encode(bit);
      partialMag = (partialMag << 1) | bit;
    }
  }
}

//----------------------------------------------------------------------------

void
TrisoupEncoder
::encodeAxiAndPosFlags(
  pcc::EntropyEncoder* const arithmeticEncoder,
  const TrisoupOriginalPCinfo oriPCinfo,
  AdaptiveBitModel (&axiFlagctx)[3][2],
  AdaptiveBitModel (&posFlagctx)[3][2],
  Vec3<int>& axiFlag_ctx,
  Vec3<int>& posFlag_ctx,
  int index)
{
  bool flag = oriPCinfo.axisFlagMask & 1 << index;
  arithmeticEncoder->encode(flag, axiFlagctx[index][axiFlag_ctx[index]]);
  if (flag) {
    arithmeticEncoder->encode(oriPCinfo.posFlag[index], posFlagctx[index][posFlag_ctx[index]]);
  }
}

//----------------------------------------------------------------------------

void
TrisoupEncoder
::prepareVertexConsistencyEncoder(
  pcc::EntropyEncoder* const arithmeticEncoder,
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

  //axiFlagctx perNode 3
  Vec3<int> axiFlag_ctx;
  Vec3<int> posFlag_ctx;
  for (int k = 0; k < 3; k++) {
    axiFlag_ctx[k] = (poscount[k] >= (cnt - 1) || poscount[k] <= 1) ? 1 : 0;
    posFlag_ctx[k] = poscount[k] >= (cnt - 1) ? 1 : 0;
  }

  uint8_t perNodeoriPCinfo = 0;
  std::array<int, 8> perNodeoriPCinfo_num = {};
  int perWidth = blockWidth >> 1;
  int twoPerWidth = perWidth << 1;
  for (int j = leaves[nodeIdxC].start(); j < leaves[nodeIdxC].end(); j++) {
    auto voxel = pointCloud[j] - leaves[nodeIdxC].pos;
    if (voxel[0] >= twoPerWidth || voxel[1] >= twoPerWidth || voxel[2] >= twoPerWidth)
      continue;
    int idx = (voxel[0] >= perWidth) << 2 | (voxel[1] >= perWidth) << 1 | (voxel[2] >= perWidth);
    perNodeoriPCinfo_num[idx]++;
  }
  for (int k = 0; k < 8; k++) {
    if ((perNodeoriPCinfo_num[k]) >= (blockWidth >> 2))
      perNodeoriPCinfo |= 1 << k;
  }
  uint8_t isFix = 0;
  for (int j = 0; j < 12; j++) {
    int uniqueIndex = segmentUniqueIndex[idxSegment + j];
    int vertexQ = TriSoupVerticesQP[uniqueIndex];
    if (vertexQ == 0)
      continue;  // skip segments that do not intersect the surface
    isFix |= !(perNodeoriPCinfo & faceIdxToSubBlockMask[LUTSegmentExtToFaceIdx[j][0]])
      << (LUTSegmentExtToFaceIdx[j][0] >> 1);
    isFix |= !(perNodeoriPCinfo & faceIdxToSubBlockMask[LUTSegmentExtToFaceIdx[j][1]])
      << (LUTSegmentExtToFaceIdx[j][1] >> 1);
  }
  const uint8_t per_axiFlag = isFix & eligible;
  Vec3<bool> per_posFlag;
  for (int k = 0; k < 3; k++) {
    int axisIdx = k << 1;  //0 2 4
    per_posFlag[k] = perNodeoriPCinfo & faceIdxToSubBlockMask[axisIdx];
  }
  oriPCinfo.axisFlagMask = per_axiFlag;
  oriPCinfo.posFlag = per_posFlag;
  bool multiFlag = 0;
  multiFlag = per_axiFlag;
  oriPCinfo.nodeMultiFlag = multiFlag;
  arithmeticEncoder->encode(multiFlag, multiFlagctx);
  if (multiFlag) {
    if (eligible > 4) {
      encodeAxiAndPosFlags(arithmeticEncoder, oriPCinfo, axiFlagctx, posFlagctx, axiFlag_ctx, posFlag_ctx, 2);
    }
    if ((eligible & 3) > 2) {
      encodeAxiAndPosFlags(arithmeticEncoder, oriPCinfo, axiFlagctx, posFlagctx, axiFlag_ctx, posFlag_ctx, 1);
    }
    int b = !(eligible & 1) + (eligible == 4);

    if (per_axiFlag >> b + 1)
      encodeAxiAndPosFlags(arithmeticEncoder, oriPCinfo, axiFlagctx, posFlagctx, axiFlag_ctx, posFlag_ctx, b);
    else
      arithmeticEncoder->encode(oriPCinfo.posFlag[b], posFlagctx[b][posFlag_ctx[b]]);
  }
}
//----------------------------------------------------------------------------
void
TrisoupEncoder
::CentroidsResidualEncoder(
  int stepQcentro,
  int dominantAxis,
  int scaleQ,
  Vec3<int32_t>& blockCentroid,
  int vtxCount,
  const TrisoupNodeEncoder& leaf,
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

    // encode centroid residual
    int res = 0;
    resCentroQ = determineCentroidResidual(normalV, blockCentroid, leaf.pos, pointCloud,
        leaf.start(), leaf.end(), lowBound, highBound, stepQcentro,
        scaleQ, res, vtxCount, blockWidth, eligDoubleCentro, dominantAxis, isDoubleCentro);
    if (eligDoubleCentro && isDoubleCentro) resCentroQ = 0;

    encodeCentroidResidual(resCentroQ, arithmeticEncoder, ctxtMemOctree, resCentroQPred, ctxMinMax,
      lowBoundSurface, highBoundSurface, lowBound, highBound);
    if (eligDoubleCentro && resCentroQ == 0)
      arithmeticEncoder->encode(isDoubleCentro, ctxtMemOctree.ctxDoubleCentroid);
    else
      isDoubleCentro = false;
  }

  if (isDoubleCentro) {
    std::vector<Vec3<int32_t>> vertices0;
    std::vector<Vec3<int32_t>> vertices1;
    for (int j = 0; j < vtxCount; j++) {
      if (subCentroMask  & 1 << j)
        vertices0.push_back(neVertex.vertices[j]);
      else
        vertices1.push_back(neVertex.vertices[j]);
    }
    subCentroid[0] = getSubCentroid(vertices0);
    subCentroid[1] = getSubCentroid(vertices1);
  }

  // dequantize and apply residual
  int resCentroDQ = 0;
  if (resCentroQ) {
    resCentroDQ = std::abs(resCentroQ) * stepQcentro + DZ - half;
  }

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
TrisoupEncoder
::GenerateCentroidsInNodeRasterScanEncoder(
  const TrisoupNodeEncoder& leaf,
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
    if (vtxCount >= 6 && isCentroidResActivated)
      eligibleDoubleCentro(
        gCenter, neVertex.dominantAxis, neVertex.vertices, blockWidth, subCentroMask, eligDoubleCentro);

    // centroid residual
    bool isDoubleCentro = false;
    std::array<Vec3<int32_t>, 2> subCentroid;
    CentroidsResidualEncoder(stepQcentro, neVertex.dominantAxis, scaleQ, blockCentroid, vtxCount, leaf, neVertex.vertices,
      neiCentroid, eligDoubleCentro, isDoubleCentro, neVertex, subCentroMask, subCentroid);

    bool boundaryInside = checkIfPointIsInsideNode(blockCentroid, blockWidth << kTrisoupFpBits);

  // store centroid residual for next frame
    ncVertex = { true, blockCentroid, isDoubleCentro == 1, subCentroid, subCentroMask, boundaryInside, gCenter };

    if (isFaceVertexActivated)
      generateFaceVerticesInNodeRasterScanEncoder(nodeIdx);
  }
}

//----------------------------------------------------------------------------
void
TrisoupEncoder
::FinalizeNonCloseNodesEncoder(
  const TrisoupNodeEncoder& leaf,
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
    CentroidsResidualEncoder(stepQcentro, dominantAxis, scaleQ, blockCentroid, 2, leaf, neVertex.vertices,
      neiCentroid, false, isDoubleCentro, neVertex, subCentroMask, subCentroid);

    bool boundaryInside = checkIfPointIsInsideNode(blockCentroid, blockWidth << kTrisoupFpBits);

    // store centroid residual for next frame
    cVerts[nodeIdx] = { true, blockCentroid, false,  {}, 0, boundaryInside , gCenter };

    if (isFaceVertexActivated)
      generateFaceVerticesInNodeRasterScanEncoder(nodeIdx);
  }
}

//----------------------------------------------------------------------------

void
TrisoupEncoder
::generateFaceVerticesInNodeRasterScanEncoder(const int nodeIdx)
{
  const int32_t tmin1 = 2 * 4;
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

    int neVtxBoundaryFace = 0; // count edge vertices included in the face
    for (int k = 0; k < eVerts[i].vertices.size(); k++)
      neVtxBoundaryFace += eVerts[i].vertices[k][axis] + kTrisoupFpHalf == 0;

    int eIdxCurr[2] = { -1, -1 };
    int eIdxNei[2] = { -1, -1 };
    Vec3<int32_t > fVert[2];

    if (1 == neVtxBoundaryFace && eVerts[i].vertices.size() == 2 && eVerts[ii].vertices.size() == 2) {
      findTrisoupFaceVertex(cVerts[i].pos, cVerts[ii].pos, axis, fVert);
      placeFaceVertex(eVerts[i], axis, fVert[0], eIdxCurr);
      placeFaceVertex(eVerts[ii], axis, fVert[1], eIdxNei);
      fVerts[i].EdgeVerIdxBeforeFaceVer.push_back(eIdxCurr[0]);
      fVerts[i].vertices.push_back(fVert[0]);
      fVerts[ii].EdgeVerIdxBeforeFaceVer.push_back(eIdxNei[0]);
      fVerts[ii].vertices.push_back(fVert[1]);
      continue;
    }
    if (2 != neVtxBoundaryFace && 3 != neVtxBoundaryFace)
      continue;

    findTrisoupFaceVertex(cVerts[i].pos, cVerts[ii].pos, axis, fVert);
    placeFaceVertex(eVerts[i], axis, fVert[0], eIdxCurr);
    placeFaceVertex(eVerts[ii], axis, fVert[1], eIdxNei);

    if (-1 == eIdxCurr[0] || -1 == eIdxCurr[1])
      continue;

    // c0, c1, and face vertex is on the same side of the surface
    bool valid = validateFaceVertex(eVerts[i], cVerts[i], cVerts[ii], eIdxCurr[0], eIdxCurr[1], fVert[0]);
    if (valid) {
      bool faceConnected;
      int32_t weight1 = 0;

      // current-node
      for (int k = leaves[i].start(); k < leaves[i].end(); k++) {
        Vec3<int32_t> dist = fVert[0] - (pointCloud[k] - leaves[i].pos << kTrisoupFpBits);
        int32_t d = dist.abs().max() + kTrisoupFpHalf >> kTrisoupFpBits;
        weight1 += d < tmin1;
      }

      // nei-node
      for (int k = leaves[ii].start(); k < leaves[ii].end(); k++) {
        Vec3<int32_t> dist = fVert[1] - (pointCloud[k] - leaves[ii].pos << kTrisoupFpBits);
        int32_t d = dist.abs().max() + kTrisoupFpHalf >> kTrisoupFpBits;
        weight1 += d < tmin1;
      }

      faceConnected = weight1 > 0;
      arithmeticEncoder->encode((int)faceConnected, ctxFaces);

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
TrisoupEncoder
::initEncoder()
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
    numPointsSlabBlocks.reserve(numSlabBlocksPerSlab);
  }

  recPointCloud.resize(pointCloud.getPointCount());

  if (useLocalAttr) {
    localPointCloud.addRemoveAttributes(pointCloud);
    recPointCloud.addRemoveAttributes(pointCloud);
  }

  // determine vertex quantization step
  midBlock = (blockWidth - 1) << 7; // << 8) >> 1;
  maxBlock = (blockWidth - 1) << 8;

  // quantization mapping par half segment
  for (int QP = 0; QP < 55; QP++) {
    int stepInv = LUT_QP_Trisoup_inv[QP];
    LUT_maxMagnitude[QP] = (midBlock)*stepInv >> 21;

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
TrisoupEncoder
::sliceEncoder(bool isFinalPass)
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

      // for TriSoup Vertex by the encoder
      int countNearPoints = 0;
      int distanceSum = 0;

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

          // determine TriSoup Vertex current point cloud
          if (!(isInter && anyNodeSkip))
            for (int j = leafNei.start(); j < leafNei.end(); j++) {
              Vec3<int> voxel = pointCloud[j];
              int pointOK = voxel[dir1] == pos1 && voxel[dir2] == pos2;
              countNearPoints += pointOK;
              distanceSum += pointOK * voxel[dir0];
            }

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
      distanceSum -= countNearPoints * currWedgePos[dir0];
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

      // determine TriSoup Vertex by the encoder
      int QPedgeWithSkip = anyNodeSkip ? 0 : QPedge;
      TriSoupEdgeLocalQP.push_back(QPedgeWithSkip);
      int8_t vertexQ = 0;
      bool isSkippedVertex = false;
      if (anyNodeSkip) {
        if (countNearPointsPred > thVertexDeterminationSkip)
          vertexQ = quantizeQP(distanceSumPred, countNearPointsPred, QPedgeWithSkip);
        isSkippedVertex = true;
      }
      else if (countNearPoints > thVertexDeterminationEnc) {
        vertexQ = quantizeQP(distanceSum, countNearPoints, QPedgeWithSkip);
      }
      TriSoupVerticesQP.push_back(vertexQ);
      TriSoupVertices2bits.push_back(map2bits(vertexQ, QPedgeWithSkip));  // map edge to two bits

      // determine TriSoup Vertex inter prediction
      if (isInter) {
        int8_t vertexPos = 0;
        if (countNearPointsPred > 0) {
          vertexPos =
            quantizeQP(distanceSumPred, countNearPointsPred, QPedgeWithSkip);
        }
        TriSoupVerticesPred.push(vertexPos);
        TriSoupVerticesIsSkipped.push(isSkippedVertex); //skip status
      }

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

        // encode vertex if not skipped
        vertex = TriSoupVerticesQP[firstVertexToCode];
        if (!(isInter && TriSoupVerticesIsSkipped.front()))
          encodeOneTriSoupVertexRasterScan(
            vertex, arithmeticEncoder, ctxtMemOctree, TriSoupVertices2bits,
            neighbNodes[firstVertexToCode], pattern, interPredictor,
            colocatedVertex, firstVertexToCode, QPedge);

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
          prepareVertexConsistencyEncoder(arithmeticEncoder, oriPCinfo);

        GenerateCentroidsInNodeRasterScanEncoder(   // dequantization of vertices is here
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
            FinalizeNonCloseNodesEncoder(leaves[nonClosedIdx], isCentroidResActivated, nonClosedIdx);
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
                      numPointsSlabBlocks.clear();
                      slabBlockIdx = localSlabBlockIdxStart[0];
                      for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
                        for (int z = 0; z < numSlabBlocksPerDim[2]; ++z, ++slabBlockIdx) {
                          int numPoints = laPoints[slabIdxMod2][slabBlockIdx].size();
                          numPointsSlabBlocks.emplace_back(numPoints);
                          if (numPoints) {
                            int size = localPointCloud.size();
                            localPointCloud.resize(size + numPoints);
                            for (int i = 0; i < numPoints; ++i)
                              localPointCloud[size++] = laPoints[slabIdxMod2][slabBlockIdx][i];
                            laPoints[slabIdxMod2][slabBlockIdx].clear(); // only clear, may be reused
                          }
                        }
                      processLocalAttributesEncoder(recPointCloud, false);
                      localPointCloud.clear();
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
                recPointCloud, leaf.isSkiped);
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
TrisoupEncoder
::finishSliceEncoder()
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
      numPointsSlabBlocks.clear();
      int slabBlockIdx = localSlabBlockIdxStart[0];
      for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
        for (int z = 0; z < numSlabBlocksPerDim[2]; ++z, ++slabBlockIdx) {
          int numPoints = laPoints[slabIdxMod2][slabBlockIdx].size();
          numPointsSlabBlocks.emplace_back(numPoints);
          if (numPoints) {
            int size = localPointCloud.size();
            localPointCloud.resize(size + numPoints);
            for (int i = 0; i < numPoints; ++i)
              localPointCloud[size++] = laPoints[slabIdxMod2][slabBlockIdx][i];
            laPoints[slabIdxMod2][slabBlockIdx] = std::vector<point_t>(); // just release memory
          }
        }
      processLocalAttributesEncoder(recPointCloud, true);
      localPointCloud.clear();
      // next slice/slab
      slabIdxMod2 = slabIdxMod2 ^ 1;
      localSlabBlockIdxStart[0] = 0; //numSlabBlocksPerDim[1] * numSlabBlocksPerDim[2]; // = localSlabBlockIdxStart[2];
      localSlabBlockStart[0] += slabBlockSize[0];
      localSlabBlockStart[2] = localSlabBlockStart[1] = 0;
    }
  }

  // copy reconstructed point cloud to point cloud
  recPointCloud.resize(nRecPoints);
  pointCloud.resize(0);
  pointCloud.swap(recPointCloud);

  ////
  // The following is to render single slab-block for full slice
  if (!useLocalAttr) {
    processGlobalAttributesEncoder(pointCloud);
  }
  clearTrisoupElements();
}

//----------------------------------------------------------------------------

void
TrisoupEncoder
::processLocalAttributesEncoder(
  PCCPointSet3& recPointCloud, bool isLast)
{
  //auto allocatedSizeLocal = localPointCloud.size();
  auto nRecPointsLocal = localPointCloud.size();
  encoder->processNextSlabAttributes(localPointCloud,
                      {localSlabBlockStart[0], 0, 0}, numPointsSlabBlocks, isLast);

  if (recPointCloud.getPointCount() < nRecPoints + nRecPointsLocal)
    recPointCloud.resize(nRecPoints + nRecPointsLocal + PC_PREALLOCATION_SIZE);

  recPointCloud.setFromPartition(localPointCloud, 0, nRecPointsLocal, nRecPoints);
  nRecPoints += nRecPointsLocal;
}

//----------------------------------------------------------------------------

void
TrisoupEncoder
::processGlobalAttributesEncoder(
  PCCPointSet3& recPointCloud)
{
  std::vector<int32_t> numPoints {nRecPoints};
  encoder->processNextSlabAttributes(pointCloud,
                      {0, 0, 0}, numPoints, true);
}

//============================================================================

}  // namespace pcc
