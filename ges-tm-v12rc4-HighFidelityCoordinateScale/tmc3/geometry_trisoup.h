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

#include <cstdint>
#include <vector>
#include <cstring>

#include "PCCPointSet.h"
#include "geometry_octree.h"
#include <map>

#define PC_PREALLOCATION_SIZE 200000

static const int precDivA = 30;
static const int LUTdiv32[33] = { // 12 bits precision
  32767, 4096, 2048, 1365, 1024, 819, 683, 585, 512, 455, 410, 372, 341, 315, 293, 273, 256, 241, 228, 216, 205, 195, 186, 178, 171, 164, 158, 152, 146, 141, 137, 132, 128 };

static const int LUT_QP_Trisoup[55] = { // 8 bit precision; QP in [0,54]
  64,    72,    81,    91,    102,   114,   128,   144,   161,   181,   203,
  228,   256,   287,   323,   362,   406,   456,   512,   575,   645,   724,
  813,   912,   1024,  1149,  1290,  1448,  1625,  1825,  2048,  2299,  2580,
  2896,  3251,  3649,  4096,  4598,  5161,  5793,  6502,  7298,  8192,  9195,
  10321, 11585, 13004, 14596, 16384, 18390, 20643, 23170, 26008, 29193, 32768 };

static const int64_t LUT_QP_Trisoup_inv[55] = { // 8 bit precision; QP in [0,54]
  32768, 29127, 25891, 23046, 20560, 18396, 16384, 14564, 13026, 11586, 10331,
  9198, 8192, 7307, 6493, 5793, 5165, 4599, 4096, 3647, 3251, 2897, 2580, 2300, 2048,
  1825, 1626, 1448, 1291, 1149, 1024, 912, 813, 724, 645, 575, 512, 456, 406, 362, 323,
  287, 256, 228, 203, 181, 161, 144, 128, 114, 102, 91, 81, 72, 64 };

namespace pcc {

//============================================================================
// The number of fractional bits used in trisoup triangle voxelisation
const int kTrisoupFpBits = 8;

// The value 1 in fixed-point representation
const int kTrisoupFpOne = 1 << (kTrisoupFpBits);
const int kTrisoupFpHalf = 1 << (kTrisoupFpBits - 1);

//============================================================================

struct codeVertexCtxInfo {
  int ctxE;
  int ctx0;
  int ctx1;
  int direction;
  int pattern = 0;
  int patternClose = 0;
  int patternClosest = 0;
  int nclosestPattern = 0;
  int missedCloseStart;
  int nClosestStart;
  int neighbEdge;
  int neighbEnd;
  int orderedPclosePar;
  int orderedPcloseParPos;
  int nNeiEdge = 0;
};

static const uint8_t faceIdxToSubBlockMask[6] = {
  1 << 0 | 1 << 1 | 1 << 2 | 1 << 3,
  1 << 4 | 1 << 5 | 1 << 6 | 1 << 7,
  1 << 0 | 1 << 1 | 1 << 4 | 1 << 5,
  1 << 2 | 1 << 3 | 1 << 6 | 1 << 7,
  1 << 0 | 1 << 2 | 1 << 4 | 1 << 6,
  1 << 1 | 1 << 3 | 1 << 5 | 1 << 7
};

static const int LUTSegmentExtToFaceIdx[12][2] = {
  {2, 4}, {0, 4}, {3, 4}, {1, 4}, {0, 2}, {0, 3},
  {1, 3}, {1, 2}, {2, 5}, {0, 5}, {3, 5}, {1, 5}
};

static const int towardOrAway[18] = { // 0 = toward; 1 = away
 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

static const int mapping18to9[3][9] = {
  { 0, 1, 2, 3,  4, 15, 14, 5,  7},
  { 0, 1, 2, 3,  9, 15, 14, 7, 12},
  { 0, 1, 2, 9, 10, 15, 14, 7, 12}
};

//============================================================================

inline bool isCodedBit(int partialMag, int bitRemaining, int maxMag)
{
  return ((partialMag << 1) | 1) << std::max(0, bitRemaining) <= maxMag;
}

//============================================================================

void constructCtxInfo(
  codeVertexCtxInfo& ctxInfo,
  int neigh,
  std::array<int, 18>& patternIdx,
  std::vector<int8_t>& TriSoupVertices);

void constructCtxPresence(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int8_t colocatedVertex
);

void constructCtxPos1(
  int& ctxMap1,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int8_t colocatedVertex,
  int interPredictor2bit);

void constructCtxPos2(
  int& ctxMap2Base,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int shiftMag,
  int v,
  int8_t colocatedVertex,
  int blockWidthLog2,
  int& goodColo);

void constructCtxPos3(
  int ctxMap2Base,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int shiftMag,
  int v,
  int8_t colocatedVertex,
  int& goodColo);

void constructCtxPos4(
  int ctxMap2Base,
  int& ctxMap2,
  int& ctxInter,
  codeVertexCtxInfo& ctxInfo,
  bool isInter,
  int8_t TriSoupVerticesPred,
  int shiftMag,
  int v,
  int8_t colocatedVertex,
  int& goodColo);

//============================================================================
// index in edgesNeighNodes of neighboring nodes for each edge
static const int edgesNeighNodesIdx[3][4] = {
  {6, 4, 0, 2}, // along z
  {6, 2, 5, 1}, // along y
  {6, 4, 5, 3}, // along x
};

// axis direction
static const int axisdirection[3][3] = { {2,0,1}, {1, 0,2}, {0, 1, 2 } }; // z, y ,x

// edge index in existing order
static const int edgeIdx[3][4] = {
  {4, 5, 6, 7},
  {1, 3, 9, 11},
  {0, 2, 8, 10},
};

static const uint16_t neighMask[3][4] = {
  {0x4001, 0x4002, 0x4004, 0x4008},
  {0x2001, 0x2002, 0x2004, 0x2008},
  {0x0001, 0x0002, 0x0004, 0x0008},
};

// index in edgesNeighNodes of other neighboring nodes for the wedge
static const int wedgeNeighNodesIdx[3][4] = {
  {1, 5, 3, 7}, // along z
  {0, 3, 4, 7}, // along y
  {0, 1, 2, 7}, // along x
};

// edge index in existing order
static const int wedgeNeighNodesEdgeIdx[3][4] = {
  {7, 4, 5,  6},
  {3, 9, 1, 11},
  {2, 8, 0, 10},
};

static const uint16_t wedgeNeighMask[3][4] = {
  {0x0800, 0x0100, 0x0200, 0x0400},
  {0x0200, 0x0400, 0x0100, 0x0800},
  {0x0200, 0x0400, 0x0100, 0x0800},
};

static const uint16_t toPrevEdgeNeighMask[3][4] = {
  {0x0010, 0x0020, 0x0040, 0x0080},
  {0x0010, 0x0020, 0x0040, 0x0080},
  {0x0010, 0x0020, 0x0040, 0x0080},
};

// neighbourhood staic tables
// ---------    8-bit pattern = 0 before, 1-4 perp, 5-12 others
static const int localEdgeindex[12][11] = {
  { 4,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 0
  { 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 1
  { 1,  5,  4,  9,  0,  8, -1, -1, -1, -1, -1}, // vertex 2
  { 0,  7,  4,  8,  2, 10,  1,  9, -1, -1, -1}, // vertex 3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 4
  { 1,  0,  9,  4, -1, -1, -1, -1, -1, -1, -1}, // vertex 5
  { 3,  2,  0, 10, 11,  9,  8,  7,  5,  4, -1}, // vertex 6
  { 0,  1,  2,  8, 10,  4,  5, -1, -1, -1, -1}, // vertex 7
  { 4,  9,  1,  0, -1, -1, -1, -1, -1, -1, -1}, // vertex 8
  { 4,  0,  1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 9
  { 5,  9,  1,  2,  8,  0, -1, -1, -1, -1, -1}, // vertex 10
  { 7,  8,  0, 10,  5,  2,  3,  9,  1, -1, -1}  // vertex 11
};
static const int patternIndex[12][11] = {
  { 3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 0
  { 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 1
  { 2,  3,  5,  8, 15, 17, -1, -1, -1, -1, -1}, // vertex 2
  { 2,  3,  5,  8,  9, 12, 15, 17, -1, -1, -1}, // vertex 3
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 4
  { 1,  7, 10, 14, -1, -1, -1, -1, -1, -1, -1}, // vertex 5
  { 1,  2,  6,  9, 10, 11, 13, 14, 15, 16, -1}, // vertex 6
  { 2,  5,  8,  9, 12, 15, 17, -1, -1, -1, -1}, // vertex 7
  { 1,  4,  7, 14, -1, -1, -1, -1, -1, -1, -1}, // vertex 8
  { 1,  7, 14, -1, -1, -1, -1, -1, -1, -1, -1}, // vertex 9
  { 1,  2,  6, 14, 15, 16, -1, -1, -1, -1, -1}, // vertex 10
  { 1,  2,  6,  9, 11, 13, 14, 15, 16, -1, -1}  // vertex 11
};

//============================================================================

int findDominantAxis(
  std::vector<Vec3<int32_t>>& leafVertices,
  int blockWidth,
  Vec3<int32_t>& blockCentroid);

//============================================================================

template<typename T>
Vec3<T>
crossProduct(const Vec3<T> a, const Vec3<T> b)
{
  Vec3<T> ret;
  ret[0] = a[1] * b[2] - a[2] * b[1];
  ret[1] = a[2] * b[0] - a[0] * b[2];
  ret[2] = a[0] * b[1] - a[1] * b[0];
  return ret;
}

//============================================================================

template<int axis>
void rayTracingAlongdirection_samp1(
  std::vector<int64_t>& renderedBlock,
  int& nPointsInBlock,
  int blockWidth,
  Vec3<int32_t>& nodepos,
  int minRange[3],
  int maxRange[3],
  Vec3<int32_t>& edge1,
  Vec3<int32_t>& edge2,
  Vec3<int32_t>& s0,
  int64_t inva,
  int haloTriangle,
  int thickness)
{
  constexpr int idx1[3] = { 1, 0, 0 };
  constexpr int idx2[3] = { 2, 2, 1 };
  constexpr int idx3[3] = { 1, 1, 2 };
  constexpr int sign1[3] = { 1, -1, 1 };
  constexpr int sign3[3] = { 1, 1, -1 };

  constexpr int i0 = axis;
  constexpr int i1 = idx1[axis];
  constexpr int i2 = idx2[axis];
  constexpr int i3 = idx3[axis];

  constexpr int s1 = sign1[axis];
  constexpr int s3 = sign3[axis];

  constexpr int shift0 = 20 * (2 - i0);
  constexpr int shift1 = 20 * (2 - i1);
  constexpr int shift2 = 20 * (2 - i2);

  int32_t u0 = (s1 * (-s0[i1] * edge2[i2] + s0[i2] * edge2[i1]) * inva) >> precDivA;
  Vec3<int32_t> q0 = crossProduct(s0, edge1);
  int32_t v0 = (q0[i0] * inva) >> precDivA;
  int32_t t0 = ((edge2 * (q0 >> kTrisoupFpBits)) * inva) >> precDivA;

  int32_t u1 = ((-s1 * edge2[i2] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t v1 = ((s1 * edge1[i2] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t t1 = (s1 * (edge2[i0] * edge1[i2] - edge2[i2] * edge1[i0]) * inva) >> precDivA;

  int32_t u2 = ((s1 * edge2[i1] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t v2 = ((-s1 * edge1[i1] << kTrisoupFpBits) * inva) >> precDivA;
  int32_t t2 = (s3 * (-edge2[0] * edge1[i3] + edge2[i3] * edge1[0]) * inva) >> precDivA;

  int64_t renderedPoint1D0 = (int64_t(nodepos[0]) << 40) + (int64_t(nodepos[1]) << 20) + int64_t(nodepos[2]);
  for (int32_t g1 = minRange[i1]; g1 <= maxRange[i1]; g1++, u0 += u1, v0 += v1, t0 += t1) {
    int32_t u = u0, v = v0, t = t0;
    for (int32_t g2 = minRange[i2]; g2 <= maxRange[i2]; g2++, u += u2, v += v2, t += t2) {
      int w = kTrisoupFpOne - u - v;
      if (u >= -haloTriangle && v >= -haloTriangle && w >= -haloTriangle) {
        int32_t foundvoxel = minRange[i0] + (t + kTrisoupFpHalf >> kTrisoupFpBits);
        if (foundvoxel >= 0 && foundvoxel < blockWidth)
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(foundvoxel) << shift0) + (int64_t(g1) << shift1) + (int64_t(g2) << shift2);
        int32_t foundvoxelUp = minRange[i0] + (t + thickness + kTrisoupFpHalf >> kTrisoupFpBits);
        if (foundvoxelUp != foundvoxel && foundvoxelUp >= 0 && foundvoxelUp < blockWidth)
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(foundvoxelUp) << shift0) + (int64_t(g1) << shift1) + (int64_t(g2) << shift2);
        int32_t foundvoxelDown = minRange[i0] + (t - thickness + kTrisoupFpHalf >> kTrisoupFpBits);
        if (foundvoxelDown != foundvoxel && foundvoxelDown >= 0 && foundvoxelDown < blockWidth)
          renderedBlock[nPointsInBlock++] = renderedPoint1D0 + (int64_t(foundvoxelDown) << shift0) + (int64_t(g1) << shift1) + (int64_t(g2) << shift2);
      }
    }// loop g2
  }//loop g1
}

//============================================================================

static inline bool
checkIfPointIsInsideNode(const Vec3<int32_t> point, const int width)
{
  return point[0] >= 0 && point[0] <= width
    && point[1] >= 0 && point[1] <= width
    && point[2] >= 0 && point[2] <= width;
}

//============================================================================

void
determineRefinedGravityCenter(
  Vec3<int32_t>& blockCentroid,
  std::vector<Vec3<int32_t>>& leafVertices,
  int blockWidth,
  int &scaleQ,
  int stepQcentro);

void
eligibleDoubleCentro(
  Vec3<int32_t>& blockCentroid,
  int dominantAxis,
  std::vector<Vec3<int32_t>>& leafVertices,
  int blockWidth,
  uint16_t& subCentroMask,
  bool& eligDoubleCentro);

Vec3<int32_t>
getSubCentroid(const std::vector<Vec3<int32_t>>& vertices);

int
findOffsetForNonClosedSurface(bool posNonclosed[64]);

Vec3<int32_t>
determineCentroidNormalAndBounds(
  int& lowBound,
  int& highBound,
  int& lowBoundSurface,
  int& highBoundSurface,
  int& ctxMinMax,
  int triCount,
  Vec3<int32_t> blockCentroid,
  int dominantAxis,
  std::vector<Vec3<int32_t>>& leafVertices,
  int blockWidth,
  int stepQcentro,
  int scaleQ,
  Vec3<int32_t> neiCentroid);

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
  int vtxCount);


struct TrisoupNodeEdgeVertex {
  uint8_t dominantAxis = 0;
  std::vector<Vec3<int32_t>> vertices;
};


struct TrisoupCentroidVertex {
  bool valid;  // this represents centroid existence
  Vec3<int32_t> pos;
  bool useSubPos;
  std::array<Vec3<int32_t>, 2> subPos;
  uint16_t subCentroMask;
  bool boundaryInside;  // true if pos is inside of node boundary
  Vec3<int32_t > gravityCenter;
};


struct TrisoupNodeFaceVertex {
  std::vector<Vec3<int32_t >> vertices;
  std::vector<int> EdgeVerIdxBeforeFaceVer;
};

//============================================================================

struct TrisoupOriginalPCinfo {
  bool nodeMultiFlag = false;
  bool eligible = false;
  uint8_t axisFlagMask;
  Vec3<bool> posFlag;
};

//============================================================================

enum{
  POS_000 = 0,
  POS_W00 = 1,
  POS_0W0 = 2,
  POS_WW0 = 3,
  POS_00W = 4,
  POS_W0W = 5,
  POS_0WW = 6,
  POS_WWW = 7
};

template <typename TrisoupNode>
struct Trisoup {

  const int32_t blockWidth;
  const Vec3<int32_t> corner[8] = {
    {       0,       0,       0 },
    { blockWidth,       0,       0 },
    {       0, blockWidth,       0 },
    { blockWidth, blockWidth,       0 },
    {       0,       0, blockWidth },
    { blockWidth,       0, blockWidth },
    {       0, blockWidth, blockWidth },
    { blockWidth, blockWidth, blockWidth }
  };

  const Vec3<int32_t> offsets[8] = {
    { -blockWidth, -blockWidth, 0 },
    { -blockWidth, 0, -blockWidth },
    { -blockWidth, 0, 0 },
    { 0, -blockWidth, -blockWidth },
    { 0, -blockWidth, 0 },
    { 0, 0, -blockWidth },
    { 0, 0, 0 },
    { -blockWidth, -blockWidth, -blockWidth },
  };

  const int startCorner[12] = { POS_000, POS_000, POS_0W0, POS_W00, POS_000, POS_0W0, POS_WW0, POS_W00, POS_00W, POS_00W, POS_0WW, POS_W0W };
  const int LUTsegmentDirection[12] = { 0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1 };

  const int posTableQNegative[3][12] = {
    {1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0},
    {0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1},
    {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0}
  };

  const int posTableQPositive[3][12] = {
    {0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0},
    {1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0},
    {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}
  };

  const int eligible_th[7] = { 1, 1, 2, 2, 3, 3, 4 };

  uint8_t oriPointCloudMask[2][3] = { {240, 204, 170}, {15, 51, 85} };
  const int vertexIdx1[12] = { 0, 0, 2, 1, 0, 2, 3, 1, 4, 4, 6, 5 };
  const int vertexIdx2[12] = { 1, 2, 3, 3, 4, 6, 7, 5, 5, 6, 7, 7 };

  std::array<int, 8> edgesNeighNodes; // neighboring nodes' index
  // The 7 firsts are used for unique segments generation/iteration
  // All the 8 are used for contextual information to be used by edge entropy coder
  Vec3<int32_t> currWedgePos = { INT32_MIN,INT32_MIN,INT32_MIN };
  int lastWedgex = 0;

  std::vector<TrisoupNode> leaves;
  PCCPointSet3& pointCloud;
  const int distanceSearchEncoder;
  const bool isInter;
  const bool interSkipEnabled;
  const PCCPointSet3& compensatedPointCloud;

  // for coding
  const GeometryParameterSet& gps;
  const GeometryBrickHeader& gbh;
  GeometryOctreeContexts& ctxtMemOctree;

  AdaptiveBitModel ctxFaces;
  std::vector<TrisoupNodeEdgeVertex> eVerts;
  std::vector<TrisoupCentroidVertex> cVerts;
  std::vector<TrisoupNodeFaceVertex> fVerts;
  std::vector<std::array<int, 6>> neiNodeIdxVec;
  std::vector<bool> isNonClosedForwardNode;
  std::map<int64_t, uint8_t> triNodeEnable;

  AdaptiveBitModel multiFlagctx;
  AdaptiveBitModel axiFlagctx[3][2];
  AdaptiveBitModel posFlagctx[3][2];

  // Box
  point_t BBorig;
  point_t keyshift;

  // local variables for loop on wedges
  std::vector<int8_t> TriSoupVerticesQP;
  std::vector<int8_t> TriSoupEdgeLocalQP;
  std::vector<int8_t> TriSoupVertices2bits;
  std::queue<bool> TriSoupVerticesIsSkipped;
  std::vector<uint16_t> neighbNodes;
  std::queue<std::array<int, 18>> edgePattern;
  std::queue<int8_t> TriSoupVerticesPred;
  std::vector<int> segmentUniqueIndex;

  // for slice tracking
  int uniqueIndex = 0;
  int firstVertexToCode = 0;
  int nodeIdxC = 0;
  int firstNodeToRender = 0;

  // for edge coding
  std::queue<point_t> posForEdgeOfVertex;

  // for colocated edge tracking
  std::vector<int64_t> currentFrameEdgeKeys;
  int colocatedEdgeIdx = 0;

  // for rendering
  int nRecPoints = 0;
  int idxSegment = 0;
  std::vector<int64_t> renderedBlock;

  // for local attributes
  bool useLocalAttr;
  const point_t slabBlockSize;
  // buffer of points for a slice of slabs
  // todo: better structure; a queue?
  std::vector<std::vector<point_t> > laPoints[2];
  int slabIdxMod2 = 0;
  point_t localSlabBlockIdxStart = 0;
  point_t numSlabBlocksPerDim = -1;
  point_t localSlabBlockStart = -1;
  PCCPointSet3 localPointCloud;

  bool haloFlag;
  bool vertexMergeFlag;
  bool vertexConsistencyFlag;
  int thickness ;
  bool isCentroidResActivated;
  bool isFaceVertexActivated;

  // quantization
  int midBlock = 1;
  int maxBlock = 1;
  int LUT_maxMagnitude[55] = {};
  int LUT_nBitMagnitude[55] = {};
  int haloTriangleQP[55] = {};

  int thVertexDeterminationEnc = 0;
  int thVertexDeterminationSkip = 0;

  bool enableSkip = false;

  // constructor
  Trisoup(int blockWidth, PCCPointSet3& pointCloud,
    int distanceSearchEncoder, bool isInter, const PCCPointSet3& compensatedPointCloud,
    const GeometryParameterSet& gps, const GeometryBrickHeader& gbh,
    GeometryOctreeContexts& ctxtMemOctree,
    bool useLocalAttr, point_t slabBlockSize
    )
    : leaves()
    , blockWidth(blockWidth)
    , pointCloud(pointCloud)
    , distanceSearchEncoder(distanceSearchEncoder)
    , isInter(isInter)
    , interSkipEnabled(isInter&& gps.trisoup_skip_mode_enabled_flag)
    , compensatedPointCloud(compensatedPointCloud)
    , gps(gps)
    , gbh(gbh)
    , ctxtMemOctree(ctxtMemOctree)
    , useLocalAttr(useLocalAttr)
    , slabBlockSize(slabBlockSize)
    , currWedgePos( Vec3<int32_t>{INT32_MIN, INT32_MIN, INT32_MIN})
    , edgesNeighNodes{ 0,0,0,0,0,0,0,0 }
  {
  }

  //---------------------------------------------------------------------------

  void flush2SlabBuff(
    std::vector<int64_t>::iterator itBegingBlock,
    int nPointsInBlock,
    int slabIdxMod2,
    int slabBlockIdx,
    bool isSkipped
  )
  {
    auto last = itBegingBlock + nPointsInBlock;
    std::sort(itBegingBlock, last);

    if (!isSkipped) {
      last = std::unique(itBegingBlock, itBegingBlock + nPointsInBlock);
    }

    const decltype(&laPoints[slabIdxMod2][slabBlockIdx]) laPointsSlabBlock[8] = {
      &laPoints[slabIdxMod2][slabBlockIdx],
      &laPoints[slabIdxMod2][slabBlockIdx + 1],
      &laPoints[slabIdxMod2][slabBlockIdx + numSlabBlocksPerDim[2]],
      &laPoints[slabIdxMod2][slabBlockIdx + numSlabBlocksPerDim[2] + 1],
      &laPoints[slabIdxMod2 ^ 1][slabBlockIdx],
      &laPoints[slabIdxMod2 ^ 1][slabBlockIdx + 1],
      &laPoints[slabIdxMod2 ^ 1][slabBlockIdx + numSlabBlocksPerDim[2]],
      &laPoints[slabIdxMod2 ^ 1][slabBlockIdx + numSlabBlocksPerDim[2] + 1],
    };

    const auto slabBlockEnd = localSlabBlockStart + slabBlockSize;

    for (auto it = itBegingBlock; it != last; it++) {
      point_t p = { int(*it >> 40), int(*it >> 20) & 0xFFFFF, int(*it) & 0xFFFFF };
      int idx =
        ((p[0] >= slabBlockEnd[0]) << 2)
        + ((p[1] >= slabBlockEnd[1]) << 1)
        + ((p[2] >= slabBlockEnd[2]) << 0);
      laPointsSlabBlock[idx]->push_back(p);
    }
  }

  //---------------------------------------------------------------------------

  void flush2PointCloud(
    int& nRecPoints,
    std::vector<int64_t>::iterator itBegingBlock,
    int nPointsInBlock,
    PCCPointSet3& recPointCloud,
    bool isSkipped
  )
  {
    auto last = itBegingBlock + nPointsInBlock;
    std::sort(itBegingBlock, last);

    if (!isSkipped) {
      last = std::unique(itBegingBlock, itBegingBlock + nPointsInBlock);
    }

    // Move list of points to pointCloud
    int nPointInCloud = recPointCloud.getPointCount();
    int nPointInNode = last - itBegingBlock;

    if (nPointInCloud < nRecPoints + nPointInNode)
      recPointCloud.resize(nRecPoints + nPointInNode + PC_PREALLOCATION_SIZE);

    for (auto it = itBegingBlock; it != last; it++)
      recPointCloud[nRecPoints++] = { int(*it >> 40), int(*it >> 20) & 0xFFFFF, int(*it) & 0xFFFFF };
  }

  //---------------------------------------------------------------------------

  bool eligibleNonClosed(const Vec3<int32_t>& v1, const Vec3<int32_t>& v2, int width)
  {
    width = (width << kTrisoupFpBits) - kTrisoupFpHalf;
    int cnt = (v1[0] < 0 || v1[0] >= width) && (v1[0] == v2[0]);
    cnt += (v1[1] < 0 || v1[1] >= width) && (v1[1] == v2[1]);
    cnt += (v1[2] < 0 || v1[2] >= width) && (v1[2] == v2[2]);

    return cnt == 1;
  }

  //---------------------------------------------------------------------------

  bool getNeighboringCentroid(
    Vec3<int32_t>& neiCentroid, uint8_t& dominantAxis, int nodeIdx,
    Vec3<int32_t>& v1, Vec3<int32_t>& v2, int blockWidth, bool& backwardNei
  ) {
    int width = (blockWidth << kTrisoupFpBits) - kTrisoupFpHalf;
    int neiCentroIdx = -1;

    for (int i = 0; i < 3; ++i) {
      if (v1[i] < 0 && v2[i] < 0) {
        dominantAxis = i;
        neiCentroIdx = neiNodeIdxVec[nodeIdx][i];
        backwardNei = true;
      }
      if (v1[i] >= width && v2[i] >= width)
      {
        dominantAxis = i;
      }
    }

    if (neiCentroIdx == -1 || !cVerts[neiCentroIdx].valid)
      return false;

    neiCentroid = cVerts[neiCentroIdx].pos;
    neiCentroid[dominantAxis] -= blockWidth << kTrisoupFpBits;
    return true;
  }

  // find face vertex position from connection of centroid vertices of
  // current and neighbour nodes one by one.
  // (totally this function is called up to three times per node.)
  void findTrisoupFaceVertex(
    Vec3<int32_t> c0,
    Vec3<int32_t> c1,
    const int axis, // order : -x,-y,-z
    Vec3<int32_t >* fVert)
  {
    int W = blockWidth << kTrisoupFpBits;
    int32_t c0face = W - kTrisoupFpHalf;
    c0[axis] += W;
    int32_t denom = c0[axis] - c1[axis];

    int64_t t = 0;  // denom=0 means that c0 and c1 are on the same grid(t=0).
    if (denom)
      t = divApprox(c0face - c1[axis] << kTrisoupFpBits, denom, 0);

    Vec3<int32_t > faceVertex(c1 + (t * (c0 - c1)  + kTrisoupFpHalf >> kTrisoupFpBits));
    fVert[0] = faceVertex;
    fVert[0][axis] = -kTrisoupFpHalf;
    fVert[1] = faceVertex;
    fVert[1][axis] = c0face;
  }


  void placeFaceVertex(
    const TrisoupNodeEdgeVertex& eVerts, int axis, Vec3<int32_t >& fvert, int eIdx[2])
  {
    // if there were two or three edge vertices on the face,
    // then to select the nearest bridge segment
    // between two edge vertices from tentative face vertex of current node,
    // edge vertices within current node are already sorted,
    // and two vertices are selected which make the nearest segment from temtative face vertex.
    // if the surface has a hole within the current node,
    // the sequential number of edge vertex couldn't be found,
    // false is returned without creating a face vertex.
    int evCnt = eVerts.vertices.size();
    int distMin = 0x7fffffff; // initial value must be larger than any case.
    int facePos = fvert[axis];

    for (int evIdx = 0; evIdx < (evCnt == 3 ? 1 : evCnt); evIdx++) {
      int ev0 = evIdx;
      int ev1 = evIdx + 1;
      if (ev1 >= evCnt)
        ev1 -= evCnt;

      Vec3<int32_t> V0 = eVerts.vertices[ev0];
      Vec3<int32_t> V1 = eVerts.vertices[ev1];
      if (facePos == V0[axis] && facePos == V1[axis]) {
        Vec3<int32_t> middlePoint = (V0 + V1) >> 1;
        Vec3<int32_t> distVec = (middlePoint - fvert) >> kTrisoupFpBits;
        int dist = distVec[0] * distVec[0] + distVec[1] * distVec[1] + distVec[2] * distVec[2];
        if (distMin > dist) {
          eIdx[0] = ev0;
          eIdx[1] = ev1;
          distMin = dist;
        }
      }
    }
  }



  // finding vector:
  //   1. gravityCenter to Centroid of current node
  //   2. gravityCenter to Centroid of current neighbour node
  //   3. face vertex vector on face of current node
  // and confirm directions of these three vectors are
  // not invert with inner product.
  bool validateFaceVertex(
    const TrisoupNodeEdgeVertex& eVerts,
    const TrisoupCentroidVertex& cVerts0,
    const TrisoupCentroidVertex& cVerts1,
    int e0, int e1, Vec3<int32_t>& fVert)
  {
    // unit vector between two edge vertices on boundary face
    Vec3<int64_t> vec01 = eVerts.vertices[e1] - eVerts.vertices[e0];
    int64_t norm01 = vec01[0] * vec01[0] + vec01[1] * vec01[1] + vec01[2] * vec01[2] >> kTrisoupFpBits;

    Vec3<int64_t> ef = fVert - eVerts.vertices[e0];
    int64_t en = ef * vec01 >> kTrisoupFpBits;
    Vec3<int64_t> efTangent = norm01 ? norm01 * ef - (en * vec01)  >> kTrisoupFpBits : ef;

    int64_t dp0 = (cVerts0.pos - cVerts0.gravityCenter) * efTangent;
    int64_t dp1 = (cVerts1.pos - cVerts1.gravityCenter) * efTangent;
    return dp0 > 0 && dp1 > 0;
  }

  //---------------------------------------------------------------------------
  int generateTrianglesInNodeRasterScan(
    const TrisoupNode& leaf, int i,
    std::vector<int64_t>& renderedBlock,
    int thickness,
    bool isFaceVertexActivated)
  {
    Vec3<int32_t> nodepos = leaf.pos;
    int nPointsInBlock = 0;

    const auto& eVertices = eVerts[i].vertices;
    const int eVerCnt = eVertices.size();

    // Skip leaves that have fewer than 2 vertices
    if (eVerCnt < 2 || (eVerCnt == 2 && !cVerts[i].valid)) {
      for (int j = 0; j < eVerCnt; j++) {
        Vec3<int32_t> point = eVertices[j] + kTrisoupFpHalf >> kTrisoupFpBits;
        // vertex to list of points
        if (checkIfPointIsInsideNode(point, blockWidth - 1)) {
          Vec3<int64_t> renderedPoint = nodepos + point;
          renderedBlock[nPointsInBlock++] =
            (renderedPoint[0] << 40)
            + (renderedPoint[1] << 20)
            + renderedPoint[2];
        }
      }
      return nPointsInBlock;
    }

    // case use two subcentroids
    const bool useSubPos = cVerts[i].useSubPos;
    std::vector<std::vector<Vec3<int32_t>>> multinodeVertices(useSubPos ? 2 : 1);
    decltype(fVerts[i].vertices)* fVertices = isFaceVertexActivated ? &fVerts[i].vertices : nullptr;
    if (useSubPos) {
      for (int j = 0; j < eVerCnt; j++)
        multinodeVertices[cVerts[i].subCentroMask & 1 << j ? 0 : 1].push_back(eVertices[j]);
    }

    // case use one centroid and not non-closed
    if (!useSubPos && !isNonClosedForwardNode[i]) {

      for (int j = 0; j < eVerCnt; j++) {
        multinodeVertices[0].push_back(eVertices[j]);
        if (isFaceVertexActivated)
          for (int k = 0; k < fVertices->size(); k++)
            if (j == fVerts[i].EdgeVerIdxBeforeFaceVer[k])
              multinodeVertices[0].push_back((*fVertices)[k]);
      }
    }

    // case use one centroid and non-closed
    if (!useSubPos && isNonClosedForwardNode[i]) {
      for (int j = 0; j < eVerCnt; ++j) {
        if (j != 0)
          multinodeVertices[0].push_back(eVertices[j]);

        if (isFaceVertexActivated)
          for (int k = 0; k < fVertices->size(); ++k)
            if ((*fVertices)[k][0] == eVertices[j][0] || (*fVertices)[k][1] == eVertices[j][1]
                || (*fVertices)[k][2] == eVertices[j][2]) // are vert on same face ?
              multinodeVertices[0].push_back((*fVertices)[k]);

        if (j == 0)
          multinodeVertices[0].push_back(eVertices[j]);
      }
    }

    int quIndex = leaf.quIndex;
    int qpNode =
      quIndex != -1 ? ctxtMemOctree.listOfQUs[quIndex].localQP : gbh.trisoup_QP;
    int haloTriangle = haloTriangleQP[qpNode];

    for (int n = 0; n < multinodeVertices.size(); n++) {
      std::vector<Vec3<int32_t>>& nodeVertices = multinodeVertices[n];
      // Divide vertices into triangles around centroid
      // and upsample each triangle by an upsamplingFactor.
      int vtxCount = nodeVertices.size();
      Vec3<int32_t> blockCentroid = useSubPos ? cVerts[i].subPos[n] : cVerts[i].pos;
      Vec3<int32_t> v2 = blockCentroid;
      Vec3<int32_t> v1 = nodeVertices[0];
      Vec3<int32_t> foundvoxel = blockCentroid + kTrisoupFpHalf >> kTrisoupFpBits;
      if (checkIfPointIsInsideNode(foundvoxel, blockWidth - 1)) {
        Vec3<int64_t> renderedPoint = nodepos + foundvoxel;
        renderedBlock[nPointsInBlock++] =
          (renderedPoint[0] << 40) + (renderedPoint[1] << 20) + renderedPoint[2];
      }

      for (int vtxIndex = 0; vtxIndex < vtxCount; vtxIndex++) {
        int j1 = vtxIndex + 1;
        if (vtxCount == 2 && j1 >= vtxCount) {
          break;
        }
        else if (j1 >= vtxCount) {
          j1 -= vtxCount;
        }

        Vec3<int32_t> v0 = v1;
        v1 = nodeVertices[j1];


        // choose rasterization perpendicular  direction
        Vec3<int32_t> edge1 = v1 - v0;
        Vec3<int32_t> edge2 = v2 - v0;
        Vec3<int32_t> a = crossProduct(edge2, edge1) >> kTrisoupFpBits;
        Vec3<int32_t> h = a.abs();
        int directionOk = (h[0] > h[1] && h[0] > h[2]) ? 0 : h[1] > h[2] ? 1 : 2;

        // check if rasterization is valid; if not skip triangle which is too small
        if (h[directionOk] <= kTrisoupFpOne) // < 2*kTrisoupFpOne should be ok
          continue;

        int64_t inva = divApprox(int64_t(1) << precDivA, h[directionOk], 0);
        inva = a[directionOk] > 0 ? inva : -inva;

        // range
        int minRange[3];
        int maxRange[3];
        for (int k = 0; k < 3; k++) {
          minRange[k] =
            std::max(0,
              std::min(std::min(v0[k], v1[k]), v2[k]) + kTrisoupFpHalf >> kTrisoupFpBits);
          maxRange[k] =
            std::min(blockWidth - 1,
              std::max(std::max(v0[k], v1[k]), v2[k]) + kTrisoupFpHalf >> kTrisoupFpBits);
        }
        Vec3<int32_t> s0 = {
          (minRange[0] << kTrisoupFpBits) - v0[0],
          (minRange[1] << kTrisoupFpBits) - v0[1],
          (minRange[2] << kTrisoupFpBits) - v0[2]
        };

        // ensure there is enough space in the block buffer
        if (renderedBlock.size() <= nPointsInBlock + blockWidth * blockWidth)
          renderedBlock.resize(renderedBlock.size() + blockWidth * blockWidth);

        // applying ray tracing along direction
        if (directionOk == 0)
          rayTracingAlongdirection_samp1<0>(
            renderedBlock, nPointsInBlock, blockWidth, nodepos, minRange,
            maxRange, edge1, edge2, s0, inva, haloTriangle, thickness);

        if (directionOk == 1)
          rayTracingAlongdirection_samp1<1>(
            renderedBlock, nPointsInBlock, blockWidth, nodepos, minRange,
            maxRange, edge1, edge2, s0, inva, haloTriangle, thickness);

        if (directionOk == 2)
          rayTracingAlongdirection_samp1<2>(
            renderedBlock, nPointsInBlock, blockWidth, nodepos, minRange,
            maxRange, edge1, edge2, s0, inva, haloTriangle, thickness);

      }  // end loop on triangles
    }
    return nPointsInBlock;
  }

  //---------------------------------------------------------------------------

  int renderSkip(
    std::vector<int64_t>& renderedBlock,
    int nPointsInBlock,
    TrisoupNode& node
  )
  {
    int nPointsInPred = node.predEnd - node.predStart;
    // ensure there is enough space in the block buffer
    if (renderedBlock.size() <= nPointsInBlock + nPointsInPred)
      renderedBlock.resize(renderedBlock.size() + nPointsInPred);

    for (int ptIdx = node.predStart; ptIdx < node.predEnd; ++ptIdx) {
      auto& point = compensatedPointCloud[ptIdx];
      renderedBlock[nPointsInBlock++] =
        (int64_t(point[0]) << 40) + (int64_t(point[1]) << 20) + point[2];
    }
    return nPointsInBlock;
  }

  //---------------------------------------------------------------------------

  void clearTrisoupElements(void)
  {
    eVerts.clear();
    cVerts.clear();
    fVerts.clear();
    neiNodeIdxVec.clear();
    return;
  }

  //---------------------------------------------------------------------------
  int quantizeQP(int dist, int count, int QP) {
    // division on the encoder side; 8 bit precision !!! division for inter pred
    int pos = fpReduce<8>(divApprox(dist, count, 16)); // (dist << 8) / count;
    return quantizeQP(pos, QP);
  }

  int quantizeQP(int pos, int QP) {
    pos -= midBlock;

    bool sign = pos >= 0;
    pos = std::abs(pos);
    int stepInv = LUT_QP_Trisoup_inv[QP];
    int magnitude  = midBlock * stepInv >> 21;
    // division on the encoder side; 8 bit precision !!! division for inter pred
    int posQ = std::min(magnitude, pos * stepInv >> 21);

    return sign ? posQ + 1 : -posQ - 1;
  }

  // quantize to 2 bits for contextual coding
  int quantizeQP2bits(int pos) {
    pos -= midBlock;

    int bit1 = pos >= 0;
    pos = std::abs(pos);
    int bit0 = pos >= (blockWidth << 6);

    return (bit1 << 1) + (bit1 ? bit0 : !bit0);
  }

  // dequantize  accordign to QP
  int dequantizeQP(int posQ, int QP) {
    bool sign = posQ >= 0;
    posQ = std::abs(posQ);

    int step = LUT_QP_Trisoup[QP];
    int pos = (posQ - 1) * step + (step >> 1); // dequantize at center
    if (posQ * step >= midBlock) // detect extreme interval
      if (step < 128) // if quantization interval is very small, push to bounds
        pos = ((posQ - 1) * step + 3 * midBlock) >> 2;
      else
        pos = ((posQ - 1) * step + midBlock) >> 1;

    return sign ? std::min(maxBlock, midBlock + pos) : std::max(0, midBlock - pos);
  }

  // map to 2 bits from quantized value
  int map2bits(int posQ, int QP) {
    if (posQ == 0)
      return -1;

    return quantizeQP2bits(dequantizeQP(posQ, QP));
  }

protected:
  //---------------------------------------------------------------------------
  bool nextIsAvailable() const { return edgesNeighNodes[0] < leaves.size(); }

  //---------------------------------------------------------------------------
  bool changeSlice() const {
    return  (!nextIsAvailable()) || currWedgePos[0] > lastWedgex;
  }

  //---------------------------------------------------------------------------
  void goNextWedge(const std::array<bool, 8>& isNeigbourSane) {
    if (isNeigbourSane[0])
      edgesNeighNodes[7] = edgesNeighNodes[0];

    // move ++ sane neigbours
    for (int i = 0; i < 7; i++)
      edgesNeighNodes[i] += isNeigbourSane[i];

    if (edgesNeighNodes[0] >= leaves.size())
      return;

    currWedgePos = leaves[edgesNeighNodes[0]].pos - offsets[0];
    for (int i = 1; i < 7; ++i) {
      if (edgesNeighNodes[i] >= leaves.size())
        break;
      auto wedgePos = leaves[edgesNeighNodes[i]].pos - offsets[i];
      if (currWedgePos > wedgePos) {
        currWedgePos = wedgePos;
      }
    }
  }
};


//============================================================================

}  // namespace pcc
