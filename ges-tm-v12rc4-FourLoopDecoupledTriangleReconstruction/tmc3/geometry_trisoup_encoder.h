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

#include "geometry_trisoup.h"

namespace pcc {

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
  int scaleQ,
  int &res,
  int vtxCount,
  int blockWidth,
  bool eligDoubleCentro,
  int dominantAxis,
  bool& isDoubleCentro);

//----------------------------------------------------------------------------

void encodeCentroidResidual(
  int resCentroQ,
  pcc::EntropyEncoder* arithmeticEncoder,
  GeometryOctreeContexts& ctxtMemOctree,
  int resCentroQPred,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound);

//============================================================================

struct TrisoupNodeEncoder {
  Vec3<int32_t> pos;
  uint32_t _start;
  uint32_t _end;
  uint32_t predStart;
  uint32_t predEnd;
  bool isSkiped;
  // local Quality Unit
  int quIndex;

  TrisoupNodeEncoder(const PCCOctree3Node& from)
    : pos(from.pos),
    _start(from.start),
    _end(from.end),
    predStart(from.predStart),
    predEnd(from.predEnd),
    quIndex(from.quIndex),
    isSkiped(from.isSkiped)
  {}

  int start() const { return _start; }
  int end() const { return _end; }
};

//============================================================================

struct TrisoupEncoder: public Trisoup<TrisoupNodeEncoder> {
  // for rendering
  PCCPointSet3 recPointCloud;

  // For coding
  pcc::EntropyEncoder* const arithmeticEncoder;

  // For local attributes
  struct PCCTMC3Encoder3* const encoder;
  std::vector<int32_t> numPointsSlabBlocks;

  //--------------------------------------------------------------------------

  TrisoupEncoder(int blockWidth, PCCPointSet3& pointCloud,
    int distanceSearchEncoder, bool isInter, const PCCPointSet3& compensatedPointCloud,
    const GeometryParameterSet& gps, const GeometryBrickHeader& gbh,
    pcc::EntropyEncoder* const arithmeticEncoder,
    GeometryOctreeContexts& ctxtMemOctree,
    bool useLocalAttr, point_t slabBlockSize,
    struct PCCTMC3Encoder3* const encoder)
    : Trisoup(blockWidth, pointCloud, distanceSearchEncoder,
      isInter, compensatedPointCloud, gps, gbh, ctxtMemOctree, useLocalAttr,
      slabBlockSize)
    , arithmeticEncoder(arithmeticEncoder)
    , encoder(encoder)
  {}

  //---------------------------------------------------------------------------

  void initEncoder();

  void sliceEncoder(bool isFinalPass);

  void finishSliceEncoder();

  //---------------------------------------------------------------------------

  void processLocalAttributesEncoder(PCCPointSet3& recPointCloud, bool isLast);

  void processGlobalAttributesEncoder(PCCPointSet3& recPointCloud);

  //--------------------------------------------------------------------------

private:

  //--------------------------------------------------------------------------

  void  encodeOneTriSoupVertexRasterScan(
    int8_t vertexQP,
    pcc::EntropyEncoder* arithmeticEncoder,
    GeometryOctreeContexts& ctxtMemOctree,
    std::vector<int8_t>& TriSoupVertices2bits,
    int neigh,
    std::array<int, 18>& patternIdx,
    int8_t interPredictor,
    int8_t colocatedVertex,
    int i,
    int qpEdge);

  //---------------------------------------------------------------------------

  void encodeAxiAndPosFlags(
    pcc::EntropyEncoder* const arithmeticEncoder,
    const TrisoupOriginalPCinfo oriPCinfo,
    AdaptiveBitModel (&axiFlagctx)[3][2],
    AdaptiveBitModel (&posFlagctx)[3][2],
    Vec3<int>& axiFlag_ctx,
    Vec3<int>& posFlag_ctx,
    int index);

  //---------------------------------------------------------------------------

  void prepareVertexConsistencyEncoder(
    pcc::EntropyEncoder* const arithmeticEncoder,
    TrisoupOriginalPCinfo& oriPCinfo);

  //--------------------------------------------------------------------------
  void CentroidsResidualEncoder(
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
    std::array<Vec3<int32_t>, 2>& subCentroid);

  //--------------------------------------------------------------------------

  void GenerateCentroidsInNodeRasterScanEncoder(
    const TrisoupNodeEncoder& leaf,
    const std::vector<int8_t>& TriSoupVerticesQP,
    const std::vector<int8_t>& TriSoupEdgeLocalQP,
    int& idxSegment,
    const bool isCentroidResActivated,
    std::vector<int>& segmentUniqueIndex,
    TrisoupOriginalPCinfo* poriPCinfo,
    int nodeIdx);
  //--------------------------------------------------------------------------

  void FinalizeNonCloseNodesEncoder(
    const TrisoupNodeEncoder& leaf,
    const bool isCentroidResActivated,
    int nodeIdx);

//--------------------------------------------------------------------------

  void generateFaceVerticesInNodeRasterScanEncoder(const int nodeIdx);
};

//============================================================================

}  // namespace pcc
