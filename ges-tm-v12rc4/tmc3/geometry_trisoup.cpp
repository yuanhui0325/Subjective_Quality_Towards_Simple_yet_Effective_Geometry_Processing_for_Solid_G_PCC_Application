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

#include "geometry_trisoup.h"

#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"

#include "PCCTMC3Encoder.h"
#include "PCCTMC3Decoder.h"

namespace pcc {

//============================================================================

static const int LUTsqrt[13] = { 0, 256, 362, 443, 512, 572, 627, 677, 724, 768, 810, 849, 887 };

void
determineRefinedGravityCenter(
  Vec3<int32_t>& gCenter,
  std::vector<Vec3<int32_t>>& leafVertices,
  int blockWidth,
  int& scaleQ,
  int stepQcentro)

{
  // refine gravity center
  int triCount = (int)leafVertices.size();

  std::vector<int> Weights(leafVertices.size(), 0);
  int Wtotal = 0;
  for (int k = 0; k < triCount; k++) {
    int k2 = k + 1;
    if (k2 >= triCount)
      k2 -= triCount;
    Vec3<int32_t> segment = (leafVertices[k] - leafVertices[k2]).abs();
    int weight = segment[0] + segment[1] + segment[2];

    Weights[k] += weight;
    Weights[k2] += weight;
    Wtotal += 2 * weight;
  }

  Vec3<int64_t> refinedGravity = 0;
  for (int j = 0; j < triCount; j++) {
    refinedGravity += int64_t(Weights[j]) * leafVertices[j];
  }
  int64_t divisor = divApprox(int64_t(1) << 32, uint64_t(Wtotal), 0);
  for (int i = 0; i < 3; i++)
    refinedGravity[i] = fpReduce<32>(refinedGravity[i] * divisor);
  gCenter = { int(refinedGravity[0]),int(refinedGravity[1]), int(refinedGravity[2]) };

  // compute quantization step for centroid residual
  // >> bitDropped  is like / (step/256)
  int64_t a = int64_t(7680) * blockWidth * 256;  // cst * node size
  int64_t b = int64_t(LUTsqrt[triCount]) * Wtotal * stepQcentro;  // sqrt(#vertices) * triangle area * Qstep
  int ratio = divApprox(a << 16, b, 0); // cst * node size / sqrt(#vertices) / triangle area / Qstep
  scaleQ = std::max(190, ratio); // precision on 8 bits
}

// --------------------------------------------------------------------------

void
eligibleDoubleCentro(
  Vec3<int32_t>& gCenter,
  int dominantAxis,
  std::vector<Vec3<int32_t>>& leafVertices,
  int blockWidth,
  uint16_t& subCentroMask,
  bool& eligDoubleCentro)

{
  int triCount = (int)leafVertices.size();

  int32_t dominantCenter = gCenter[dominantAxis];
  int32_t maxNegative = 0;
  int32_t minPositive = blockWidth << kTrisoupFpBits;
  int negativeNum = 0;
  int positiveNum = 0;
  for (int j = 0; j < triCount; j++) {
    int32_t dominantVertex = leafVertices[j][dominantAxis];
    if (dominantVertex > dominantCenter) {
      positiveNum++;
      if (dominantVertex < minPositive)
        minPositive = dominantVertex;
    }
    else {
      subCentroMask  |= 1 << j;
      negativeNum++;
      if (dominantVertex > maxNegative)
        maxNegative = dominantVertex;
    }
  }
  int32_t distNP = minPositive - maxNegative;
  eligDoubleCentro = positiveNum >= 3 && negativeNum >= 3 && distNP >= (blockWidth << (kTrisoupFpBits - 2));
}

// --------------------------------------------------------------------------

Vec3<int32_t>
getSubCentroid(const std::vector<Vec3<int32_t>>& vertices)
{
  int triCount = (int)vertices.size();
  std::vector<int> Weigths(triCount, 0);
  int Wtotal = 0;
  for (int k = 0; k < triCount; k++) {
    int k2 = k + 1;
    if (k2 >= triCount)
      k2 -= triCount;
    Vec3<int32_t> segment = (vertices[k] - vertices[k2]).abs();
    int weight = segment[0] + segment[1] + segment[2];

    Weigths[k] += weight;
    Weigths[k2] += weight;
    Wtotal += 2 * weight;
  }
  Vec3<int64_t> blockCentroid2 = 0;
  for (int j = 0; j < triCount; j++) {
    blockCentroid2 += int64_t(Weigths[j]) * vertices[j];
  }
  int64_t divisor = divApprox(int64_t(1) << 32, uint64_t(Wtotal), 0);
  for (int i = 0; i < 3; i++)
    blockCentroid2[i] = fpReduce<32>(blockCentroid2[i] * divisor);
  Vec3<int32_t> pos = {int(blockCentroid2[0]), int(blockCentroid2[1]), int(blockCentroid2[2])};
  return pos;
}

// --------------------------------------------------------------------------

int
findOffsetForNonClosedSurface(bool posNonclosed[64]) {
  int maxIdx = 0;
  for (int i = 0; i < 63; ++i) {
    if (posNonclosed[i]) {
      if (!posNonclosed[i + 1])
        return (i - 32) << kTrisoupFpBits;
      maxIdx = i;
    }
  }
  return std::max(maxIdx - 32, 0) << kTrisoupFpBits;
}

// --------------------------------------------------------------------------

Vec3<int32_t>
determineCentroidNormalAndBounds(
  int& lowBound,
  int& highBound,
  int& lowBoundSurface,
  int& highBoundSurface,
  int & ctxMinMax,
  int triCount,
  Vec3<int32_t> blockCentroid,
  int dominantAxis,
  std::vector<Vec3<int32_t>>& leafVertices,
  int blockWidth,
  int stepQcentro,
  int scaleQ,
  Vec3<int32_t> neiCentroid)
{
  // contextual information for residual coding
  int minPos = leafVertices[0][dominantAxis];
  int maxPos = minPos;
  for (int k = 1; k < triCount; k++) {
    int vDom = leafVertices[k][dominantAxis];
    if (vDom < minPos)
      minPos = vDom;
    if (vDom > maxPos)
      maxPos = vDom;
  }

  // find normal vector
  Vec3<int64_t> accuNormal = 0;
  if (triCount != 2) {
    for (int k = 0; k < triCount; k++) {
      int k2 = k + 1;
      if (k2 >= triCount)
        k2 -= triCount;
      accuNormal += crossProduct(leafVertices[k] - blockCentroid, leafVertices[k2] - blockCentroid);
    }
  } else
    accuNormal += blockCentroid - neiCentroid;

  int64_t invNormN = irsqrt(accuNormal[0] * accuNormal[0] + accuNormal[1] * accuNormal[1] + accuNormal[2] * accuNormal[2]);
  Vec3<int32_t> normalV = accuNormal  * invNormN >> 40 - kTrisoupFpBits;

  // residual bounds
  ctxMinMax = std::min(8, int(divApprox(maxPos - minPos, stepQcentro, 0)));
  int boundL = -kTrisoupFpHalf;
  int boundH = ((blockWidth - 1) << kTrisoupFpBits) + kTrisoupFpHalf - 1;
  int m = 1;

  int half = stepQcentro >> 1;
  int DZ = 682 * half >> 10; // 2 * half / 3;

  for (; m < blockWidth; m++) {
    int resDm = m * stepQcentro;
    resDm += DZ - half;
    resDm = resDm * scaleQ  >> 8;

    Vec3<int32_t> temp = blockCentroid + (resDm * normalV >> 8); // >>6
    if (temp[0]<boundL || temp[1]<boundL || temp[2]<boundL || temp[0]>boundH || temp[1]>boundH || temp[2]> boundH)
      break;
  }
  highBound = m - 1;

  m = 1;
  for (; m < blockWidth; m++) {
    int resDm = m * stepQcentro;
    resDm += DZ - half;
    resDm = resDm * scaleQ >> 8;

    Vec3<int32_t> temp = blockCentroid + (-resDm * normalV >> 8); // >>6
    if (temp[0]<boundL || temp[1]<boundL || temp[2]<boundL || temp[0]>boundH || temp[1]>boundH || temp[2]> boundH)
      break;
  }
  lowBound = m - 1;
  lowBoundSurface = std::max(0, (blockCentroid[dominantAxis] - minPos) + kTrisoupFpHalf >> kTrisoupFpBits);
  highBoundSurface = std::max(0, (maxPos - blockCentroid[dominantAxis]) + kTrisoupFpHalf >> kTrisoupFpBits);

  return normalV;
}

//----------------------------------------------------------------------------

void
determineCentroidPredictor(
  int& resCentroQPred,
  Vec3<int32_t> normalV,
  Vec3<int32_t> blockCentroid,
  Vec3<int32_t> nodepos,
  const PCCPointSet3& compensatedPointCloud,
  int start,
  int end,
  int lowBound,
  int  highBound,
  int stepQcentro,
  int scaleQ,
  int vtxCount)
{
  resCentroQPred = -100;

  // determine quantized residual for predictor
  if (end > start) {
    resCentroQPred = 0;

    bool posNonclosed[64] = { false };

    int resPred = 0;
    int counter = 0;
    int maxD = std::max(1, stepQcentro >> 9);

    for (int p = start; p < end; p++) {
      auto point = (compensatedPointCloud[p] - nodepos) << kTrisoupFpBits;

      Vec3<int32_t> CP = crossProduct(normalV, point - blockCentroid) >> kTrisoupFpBits;
      int dist = std::max(std::max(std::abs(CP[0]), std::abs(CP[1])), std::abs(CP[2]));
      dist >>= kTrisoupFpBits;

      if ((dist << 10) <= 1774 * maxD) {
        if (vtxCount != 2) {
          int32_t w = (1 << 10) + 4 * (1774 * maxD - ((1 << 10) * dist));
          counter += w >> 10;
          resPred += (w >> 10) * ((normalV * (point - blockCentroid)) >> kTrisoupFpBits);
        } else {
          int dists = (normalV * (point - blockCentroid)) >> kTrisoupFpBits;
          posNonclosed[32 + (dists >> kTrisoupFpBits)] = true;
        }
      }
    }

    if (vtxCount == 2) {
      resPred = divApprox((int64_t(findOffsetForNonClosedSurface(posNonclosed)) << 8), scaleQ, 0);
      resPred = resPred < 0 ? 0 : resPred;
    } else if (counter) { // res is shift by kTrisoupFpBits
      resPred = divApprox((int64_t(resPred) << 8), counter * scaleQ, 0); // res is shift by +8 due to scaleQ, then  kTrisoupFpBits - 6 due to precision res and Q
    }

    int half = stepQcentro >> 1;
    int DZ = 682 * half >> 10; //2 * half / 3;

    if (abs(resPred) >= DZ) {
      resCentroQPred = divApprox(1 + abs(resPred - DZ) << 1 , stepQcentro , 0);

      if (resPred < 0) {
        resCentroQPred = -resCentroQPred;
      }
    }
    resCentroQPred = std::min(std::max(resCentroQPred, -2 * lowBound), 2 * highBound);  // res in [-lowBound; highBound] but quantization is twice better
  }
}

//============================================================================

void
constructCtxInfo(
  codeVertexCtxInfo& ctxInfo,
  int neigh,
  std::array<int, 18>& patternIdx,
  std::vector<int8_t>& TriSoupVertices2bits) {

  // node info
  ctxInfo.ctxE = (neigh & 1) + !!(neigh & 2) + !!(neigh & 4) + !!(neigh & 8) - 1; // at least one node is occupied
  ctxInfo.ctx0 = !!(neigh & 16) + !!(neigh & 32) + !!(neigh & 64) + !!(neigh & 128);
  ctxInfo.ctx1 = !!(neigh & 256) + !!(neigh & 512) + !!(neigh & 1024) + !!(neigh & 2048);
  int direction = neigh >> 13; // 0=x, 1=y, 2=z
  ctxInfo.direction = direction;

  // reorganize node neighbours of vertex independently on xyz
  ctxInfo.neighbEdge = neigh & 15;
  ctxInfo.neighbEnd = (neigh >> 4) & 15;
  if (direction == 2) {
    ctxInfo.neighbEdge = neigh & 1;
    ctxInfo.neighbEdge += ((neigh >> 3) & 1) << 1;
    ctxInfo.neighbEdge += ((neigh >> 1) & 1) << 2;
    ctxInfo.neighbEdge += ((neigh >> 2) & 1) << 3;

    ctxInfo.neighbEnd = (neigh >> 4) & 1;
    ctxInfo.neighbEnd += ((neigh >> 7) & 1) << 1;
    ctxInfo.neighbEnd += ((neigh >> 5) & 1) << 2;
    ctxInfo.neighbEnd += ((neigh >> 6) & 1) << 3;
  }

  //neighbours info
  for (int v = 0; v < 9; v++) {
    int v18 = mapping18to9[direction][v];

    if (patternIdx[v18] != -1) {
      int vertexPos2bits = TriSoupVertices2bits[patternIdx[v18]];
      if (vertexPos2bits >= 0) {
        ctxInfo.pattern |= 1 << v;
        if (towardOrAway[v18])
          vertexPos2bits = 3 - vertexPos2bits; // reverses for away
        if (vertexPos2bits >= 2)
          ctxInfo.patternClose |= 1 << v;
        if (vertexPos2bits >= 3)
          ctxInfo.patternClosest |= 1 << v;
        ctxInfo.nclosestPattern += vertexPos2bits >= 3 && v <= 4;
      }
    }
  }

  ctxInfo.missedCloseStart = !(ctxInfo.pattern & 2) + !(ctxInfo.pattern & 4);
  ctxInfo.nClosestStart = (ctxInfo.patternClosest & 1) + !!(ctxInfo.patternClosest & 2) + !!(ctxInfo.patternClosest & 4);
  if (direction == 0) {
    ctxInfo.missedCloseStart += !(ctxInfo.pattern & 8) + !(ctxInfo.pattern & 16);
    ctxInfo.nClosestStart += !!(ctxInfo.patternClosest & 8) + !!(ctxInfo.patternClosest & 16);
  }
  if (direction == 1) {
    ctxInfo.missedCloseStart += !(ctxInfo.pattern & 8);
    ctxInfo.nClosestStart += !!(ctxInfo.patternClosest & 8) - !!(ctxInfo.patternClosest & 16);
  }
  if (direction == 2) {
    ctxInfo.nClosestStart += -!!(ctxInfo.patternClosest & 8) - !!(ctxInfo.patternClosest & 16);
  }

  ctxInfo.orderedPclosePar = (((ctxInfo.pattern >> 5) & 3) << 1) + !!(ctxInfo.pattern & 128);
  ctxInfo.orderedPcloseParPos = (((ctxInfo.patternClose >> 5) & 3) << 1) + !!(ctxInfo.patternClose & 128);

  // neighbour info
  for (int v = 0; v < 18; v++)
    if (patternIdx[v] != -1)
      ctxInfo.nNeiEdge++;
}

// ---------------------------------------------------------------------------

void
constructCtxPresence(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int8_t colocatedVertex) {

  ctxMap1 = std::min(ctxInfo.nclosestPattern, 2) * 30 + (ctxInfo.neighbEdge - 1) * 2 + (ctxInfo.ctx1 == 4);    // 2*15*3 = 90 -> 7 bits
  ctxMap2 = ctxInfo.neighbEnd << 10;
  ctxMap2 |= (ctxInfo.patternClose & (0b00000110)) << 7; // perp that do not depend on direction = to start
  ctxMap2 |= ctxInfo.direction << 6;
  ctxMap2 |= (ctxInfo.patternClose & (0b00011000)) << 1; // perp that  depend on direction = to start or to end
  ctxMap2 |= (ctxInfo.patternClose & (0b00000001)) << 3;  // before
  ctxMap2 |= ctxInfo.orderedPclosePar;

  bool isInterGood = isInter &&  ctxInfo.nNeiEdge <= 3;

  ctxInter = 0;
  if (isInterGood) {
    ctxInter = 1 + (TriSoupVerticesPred != 0) * 3;

    if (ctxInfo.nNeiEdge <= 0)
      ctxInter += 1 + (colocatedVertex != 0);
  }
}

// ---------------------------------------------------------------------------

void
constructCtxPos1(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int8_t colocatedVertex,
  int interPredictor2bit) {

  int ctxFullNbounds = (4 * (ctxInfo.ctx0 <= 1 ? 0 : (ctxInfo.ctx0 >= 3 ? 2 : 1)) + (std::max(1, ctxInfo.ctx1) - 1)) * 2 + (ctxInfo.ctxE == 3);
  ctxMap1 = ctxFullNbounds * 2 + (ctxInfo.nClosestStart > 0);
  ctxMap2 = ctxInfo.missedCloseStart << 8;
  ctxMap2 |= (ctxInfo.patternClosest & 1) << 7;
  ctxMap2 |= ctxInfo.direction << 5;
  ctxMap2 |= ctxInfo.patternClose & 0b00011111;
  ctxMap2 = (ctxMap2 << 3) + ctxInfo.orderedPcloseParPos;

  bool isGoodRef = isInter && colocatedVertex != 0;
  bool isInterGood = isInter && ((isGoodRef && ctxInfo.nNeiEdge <= 0) || ctxInfo.nNeiEdge <= 4);

  ctxInter = 0;
  if (isInterGood) {
    ctxInter =
      TriSoupVerticesPred != 0 ? 1 + (TriSoupVerticesPred > 0 ? 1 : 0) : 0;
    if (ctxInter > 0) {
      ctxInter += 2 * (interPredictor2bit == 1 || interPredictor2bit == 2);
      ctxInter = ctxInter * 3 - 2;

      int goodColo = colocatedVertex != 0 && ctxInfo.nNeiEdge <= 0;
      if (goodColo)
        ctxInter += colocatedVertex > 0 ? 2 : 1;
    }
  }
}

// ---------------------------------------------------------------------------

void
constructCtxPos2(
  int& ctxMap2Base,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int shiftMag,
  int partialDec,
  int8_t colocatedVertex,
  int blockWidthLog2,
  int& goodColo){

  ctxMap2Base = ctxInfo.missedCloseStart << 8;
  ctxMap2Base |= (ctxInfo.patternClose & 1) << 7;
  ctxMap2Base |= (ctxInfo.patternClosest & 1) << 6;
  ctxMap2Base |= ctxInfo.direction << 4;
  ctxMap2Base |= (ctxInfo.patternClose & 0b00011111) >> 1;
  ctxMap2 = (ctxMap2Base << 3) + ctxInfo.orderedPcloseParPos;

  ctxInter = 0;

  if (isInter) {
    ctxInter = TriSoupVerticesPred != 0
      ? 1 + (!!((std::abs(TriSoupVerticesPred) - 1) >> shiftMag)) : 0;

    goodColo = colocatedVertex != 0 ? 1 : 0;
    goodColo = goodColo && ((colocatedVertex > 0 ? 1 : 0) == partialDec)
      && ctxInfo.nNeiEdge <= blockWidthLog2 - 2;

    if (ctxInter > 0) {
      ctxInter = ctxInter * 3 - 2;
      if (goodColo)
        ctxInter += (!!((std::abs(colocatedVertex) - 1) >> shiftMag)) ? 2 : 1;
    }
  }
}

// ---------------------------------------------------------------------------

void
constructCtxPos3(
  int ctxMap2Base,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int shiftMag,
  int partialDec,
  int8_t colocatedVertex,
  int& goodColo) {

  ctxMap2 = ctxMap2Base;

  ctxInter = 0;
  if (isInter) {
    int temp = ((std::abs(TriSoupVerticesPred) - 1) >> std::max(0,shiftMag - 1)) & 1;
    ctxInter = TriSoupVerticesPred != 0 ? 1 + temp : 0;

    goodColo = goodColo && (!!((std::abs(colocatedVertex) - 1) >> shiftMag)) == (partialDec & 1);
    ctxMap2 |= goodColo << 12;
    if (goodColo) {
      int temp2 = ((std::abs(colocatedVertex) - 1) >> std::max(0, shiftMag - 1)) & 1;
      ctxMap2 |= temp2 << 11;
    }
  }
}

// ---------------------------------------------------------------------------

void
constructCtxPos4(
  int ctxMap2Base,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int shiftMag,
  int partialDec,
  int8_t colocatedVertex,
  int& goodColo) {

  ctxMap2 = ctxMap2Base >> 1;

  ctxInter = 0;
  if (isInter) {
    int temp = ((std::abs(TriSoupVerticesPred) - 1) >> std::max(0, shiftMag - 2)) & 1;
    ctxInter = TriSoupVerticesPred != 0 ? 1 + temp : 0;

    goodColo = goodColo && (((std::abs(colocatedVertex) - 1) >> std::max(0, shiftMag - 1)) & 1) == (partialDec & 1);
    ctxMap2 |= goodColo << 11;
    if (goodColo) {
      int temp2 = ((std::abs(colocatedVertex) - 1) >> std::max(0, shiftMag - 2)) & 1;
      ctxMap2 |= temp2 << 10;
    }
  }
  else {
    ctxMap2 <<= 1;
  }
}

//============================================================================

//============================================================================
// Representation for a vertex in preparation for sorting.
struct Vertex {
  Vec3<int32_t> pos;  // position of vertex
  int32_t theta;      // angle of vertex when projected along dominant axis
  int32_t tiebreaker;  // coordinate of vertex along dominant axis
};

// ---------------------------------------------------------------------------
// Project vertices along dominant axis (i.e., into YZ, XZ, or XY plane).
// Sort projected vertices by decreasing angle in [-pi,+pi] around center
// of block (i.e., clockwise) breaking ties in angle by
// increasing distance along the dominant axis.

int findDominantAxis(
  std::vector<Vec3<int32_t>>& leafVertices2,
  int blockWidth,
  Vec3<int32_t>& blockCentroid) {

  int dominantAxis = 0;
  int triCount = leafVertices2.size();

  std::vector<Vertex> leafVertices(triCount);
  for (int n = 0; n < triCount; n++)
    leafVertices[n].pos = leafVertices2[n];

  const int Width = blockWidth << kTrisoupFpBits;
  const int sIdx1[3] = { 2,2,1 };
  const int sIdx2[3] = { 1,0,0 };
  auto leafVerticesTemp = leafVertices;

  int maxNormTri = 0;
  for (int axis = 0; axis <= 2; axis++) {
    int axis1 = sIdx1[axis];
    int axis2 = sIdx2[axis];

    // order along axis
    for (int j = 0; j < triCount; j++) {
      // compute score closckwise
      int x = leafVerticesTemp[j].pos[axis1] + kTrisoupFpHalf; // back to [0,B]^3 for ordering
      int y = leafVerticesTemp[j].pos[axis2] + kTrisoupFpHalf; // back to [0,B]^3 for ordering

      int flag3 = x <= 0;
      int score3 = Width - flag3 * y + (!flag3) * x;
      int flag2 = y >= Width;
      int score2 = 2 * Width - flag2 * x + (!flag2) * score3;
      int flag1 = x >= Width;
      int score = flag1 * y + (!flag1) * score2;
      leafVerticesTemp[j].theta = score;
      leafVerticesTemp[j].tiebreaker = leafVerticesTemp[j].pos[axis]; // stable sort if same score
    }
    std::sort(leafVerticesTemp.begin(), leafVerticesTemp.end(),
      [](Vertex v1, Vertex v2) {
        if (v1.theta > v2.theta)
          return true;  // sort in decreasing order of theta
        if (v1.theta == v2.theta && v1.tiebreaker < v2.tiebreaker)
          return true;
        return false;
      });

    // compute sum normal
    int32_t accuN = 0;
    for (int k = 0; k < triCount; k++) {
      int k2 = k == triCount - 1 ? 0 : k + 1;
      int32_t h = (leafVerticesTemp[k].pos[axis1] - blockCentroid[axis1]) * (leafVerticesTemp[k2].pos[axis2] - blockCentroid[axis2]);
      h -= (leafVerticesTemp[k].pos[axis2] - blockCentroid[axis2]) * (leafVerticesTemp[k2].pos[axis1] - blockCentroid[axis1]);
      accuN += std::abs(h);
    }

    // if sumnormal is bigger , this is dominantAxis
    if (accuN > maxNormTri) {
      maxNormTri = accuN;
      dominantAxis = axis;
      leafVertices = leafVerticesTemp;
    }
  }

  for (int n = 0; n < triCount; n++)
    leafVertices2[n] = leafVertices[n].pos;

  return dominantAxis;
}

//============================================================================

}  // namespace pcc
