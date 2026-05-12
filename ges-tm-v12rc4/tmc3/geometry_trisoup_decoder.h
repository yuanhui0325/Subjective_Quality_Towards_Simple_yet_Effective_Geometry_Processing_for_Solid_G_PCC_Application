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

int decodeCentroidResidual(
  pcc::EntropyDecoder* arithmeticDecoder,
  GeometryOctreeContexts& ctxtMemOctree,
  int resCentroQPred,
  int ctxMinMax,
  int lowBoundSurface,
  int highBoundSurface,
  int lowBound,
  int highBound);

//============================================================================

struct TrisoupNodeDecoder {
  Vec3<int32_t> pos;
  uint32_t predStart;
  uint32_t predEnd;
  bool isSkiped;
  // local Quality Unit
  int quIndex;

  TrisoupNodeDecoder(const PCCOctree3Node& from)
    : pos(from.pos),
    predStart(from.predStart),
    predEnd(from.predEnd),
    quIndex(from.quIndex),
    isSkiped(from.isSkiped)
  {}

  int start() const { return 0; }
  int end() const { return 0; }
};

//============================================================================

struct TrisoupDecoder: public Trisoup<TrisoupNodeDecoder> {
  // For coding
  pcc::EntropyDecoder* const arithmeticDecoder;

  // For local attributes
  struct PCCTMC3Decoder3* const decoder;

  //--------------------------------------------------------------------------

  TrisoupDecoder(int blockWidth, PCCPointSet3& pointCloud,
    int distanceSearchEncoder, bool isInter, const PCCPointSet3& compensatedPointCloud,
    const GeometryParameterSet& gps, const GeometryBrickHeader& gbh,
    pcc::EntropyDecoder* const arithmeticDecoder,
    GeometryOctreeContexts& ctxtMemOctree,
    bool useLocalAttr, point_t slabBlockSize,
    struct PCCTMC3Decoder3* const decoder)
    : Trisoup(blockWidth, pointCloud, distanceSearchEncoder,
      isInter, compensatedPointCloud, gps, gbh, ctxtMemOctree, useLocalAttr,
      slabBlockSize)
    , arithmeticDecoder(arithmeticDecoder)
    , decoder(decoder)
  {}

  //--------------------------------------------------------------------------

  void initDecoder();

  void sliceDecoder(bool isFinalPass);

  void finishSliceDecoder();

  //--------------------------------------------------------------------------

  void processLocalAttributesDecoder(PCCPointSet3& recPointCloud, bool isLast);

  void processGlobalAttributesDecoder(PCCPointSet3& recPointCloud);

  //--------------------------------------------------------------------------

private:

  //--------------------------------------------------------------------------

  void  decodeOneTriSoupVertexRasterScan(
    pcc::EntropyDecoder* const arithmeticDecoder,
    GeometryOctreeContexts& ctxtMemOctree,
    std::vector<int8_t>& TriSoupVerticesQP,
    std::vector<int8_t>& TriSoupVertices2bits,
    int neigh,
    std::array<int, 18>& patternIdx,
    int8_t interPredictor,
    int8_t colocatedVertex,
    int i,
    int qpEdge);

  //---------------------------------------------------------------------------

  void decodeFlags(
    pcc::EntropyDecoder* const decoder,
    uint8_t& axiFlag,
    Vec3<bool>& posFlag,
    AdaptiveBitModel (&axiFlagctx)[3][2],
    AdaptiveBitModel (&posFlagctx)[3][2],
    Vec3<int>& axiFlag_ctx,
    Vec3<int>& posFlag_ctx,
    int idx);

  //---------------------------------------------------------------------------

  void prepareVertexConsistencyDecoder(
    pcc::EntropyDecoder* const arithmeticDecoder,
    TrisoupOriginalPCinfo& oriPCinfo);

  //---------------------------------------------------------------------------
  void CentroidsResidualDecoder(
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
    std::array<Vec3<int32_t>, 2>& subCentroid);

  //---------------------------------------------------------------------------
  void GenerateCentroidsInNodeRasterScanDecoder(
    const TrisoupNodeDecoder& leaf,
    const std::vector<int8_t>& TriSoupVerticesQP,
    const std::vector<int8_t>& TriSoupEdgeLocalQP,
    int& idxSegment,
    const bool isCentroidResActivated,
    std::vector<int>& segmentUniqueIndex,
    TrisoupOriginalPCinfo* poriPCinfo,
    int nodeIdx);

  //--------------------------------------------------------------------------
  void FinalizeNonCloseNodesDecoder(
    const TrisoupNodeDecoder& leaf,
    const bool isCentroidResActivated,
    int nodeIdx);
  //--------------------------------------------------------------------------

  void generateFaceVerticesInNodeRasterScanDecoder(const int nodeIdx);

  //---------------------------------------------------------------------------

};

//============================================================================

}  // namespace pcc
