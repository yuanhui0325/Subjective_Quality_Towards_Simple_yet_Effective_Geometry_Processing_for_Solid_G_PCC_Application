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

#include "AttributeDecoder.h"

#include "AttributeCommon.h"
#include "constants.h"
#include "entropy.h"
#include "hls.h"
#include "io_hls.h"
#include "RAHT.h"

namespace pcc {

using namespace RAHT;

//----------------------------------------------------------------------------

PCCResidualsDecoder::PCCResidualsDecoder(
  const AttributeBrickHeader& abh, const AttributeContexts& ctxtMem, const PredModeContexts& ctxtMemPredMode)
  : AttributeContexts(ctxtMem)
  , PredModeContexts(ctxtMemPredMode)
{}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::start(
  const SequenceParameterSet& sps, const char* buf, int buf_len)
{
  arithmeticDecoder.setBuffer(buf_len, buf);
  arithmeticDecoder.setBypassBinCodingWithoutProbUpdate(sps.bypass_bin_coding_without_prob_update);
  arithmeticDecoder.start();
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::stop()
{
  arithmeticDecoder.stop();
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeZeroBlock(int coeffCnt, int numCoeffNot0, int numCoeffTotal, bool zeroParent, bool enableAveragePrediction)
{
  int ctx = (numCoeffNot0 > 0) + (numCoeffNot0 > 1) + (numCoeffNot0 > 3);
  ctx = 2 * ctx + zeroParent;
  ctx = 3 * ctx + (4 * numCoeffNot0 >= numCoeffTotal) + (2 * numCoeffNot0 >= numCoeffTotal);
  ctx = 3 * ctx + (coeffCnt > 2) + (coeffCnt > 5);

  return arithmeticDecoder.decode(ctxZeroBlock[ctx][enableAveragePrediction]);
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeZeroCoeffs(int numCoeffNot0, int numCoeffTotal, int blockIndex, int existsNoZeroInBlock, bool enableAveragePrediction)
{
  static constexpr std::array<int, 8> blockIndexMap {0, 2, 1, 3, 0, 3, 3, 3};

  int ctx = (numCoeffNot0 > 0) + (numCoeffNot0 > 1) + (numCoeffNot0 > 3);
  ctx = 3 * ctx + (4 * numCoeffNot0 >= numCoeffTotal) + (2 * numCoeffNot0 >= numCoeffTotal);

  ctx = 4 * ctx + blockIndexMap[blockIndex];
  ctx = 3 * ctx + existsNoZeroInBlock;

  return arithmeticDecoder.decode(ctxZeroCoeffs[ctx][enableAveragePrediction]);
}

//----------------------------------------------------------------------------

int
PCCResidualsDecoder::decodeCoeffMag(int k1, int k2, int k3, bool enableAveragePrediction, int ctx0)
{
  int ctx = enableAveragePrediction + 2 * (k3 == 0) * ctx0;
  if (!arithmeticDecoder.decode(ctxCoeffGtN[0][k1][ctx]))
    return 0;

  if (!arithmeticDecoder.decode(ctxCoeffGtN[1][k2][enableAveragePrediction]))
    return 1;

  int coeff_abs_minus2 = arithmeticDecoder.decodeExpGolomb(
    1, ctxCoeffRemPrefix[k3], ctxCoeffRemSuffix[k3]);

  return coeff_abs_minus2 + 2;
}

//----------------------------------------------------------------------------

void
PCCResidualsDecoder::decode(int32_t value[3], bool enableAveragePrediction, int numCoeffNot0, int numCoeffTotal, bool is420)
{
  int b0 = 1, b1 = 1, b2 = 1, b3 = 1;
  int ctx0 = (4 * numCoeffNot0 >= numCoeffTotal) + (2 * numCoeffNot0 >= numCoeffTotal);
  if (!is420) {
    value[1] = decodeCoeffMag(0, 0, 1, enableAveragePrediction, ctx0);
    b0 = value[1] == 0;
    b1 = value[1] <= 1;
    value[2] = decodeCoeffMag(1 + b0, 1 + b1, 1, enableAveragePrediction, ctx0);
    b2 = value[2] == 0;
    b3 = value[2] <= 1;
  }
  else {
    value[1] = 0;
    value[2] = 0;
  }
  value[0] = decodeCoeffMag(3 + (b0 << 1) + b2, 3 + (b1 << 1) + b3, 0, enableAveragePrediction, ctx0);

  if (b0 && b2)
    value[0] += 1;

  if (value[0] && arithmeticDecoder.decode())
    value[0] = -value[0];
  if (value[1] && arithmeticDecoder.decode())
    value[1] = -value[1];
  if (value[2] && arithmeticDecoder.decode())
    value[2] = -value[2];
}

//----------------------------------------------------------------------------

int32_t
PCCResidualsDecoder::decode()
{
  auto mag = decodeCoeffMag(0, 0, 0,false, 0) + 1;
  bool sign = arithmeticDecoder.decode();
  return sign ? -mag : mag;
}

//----------------------------------------------------------------------------

bool
PCCResidualsDecoder::decodeSkip(int64_t averageDGeom)
{
  return arithmeticDecoder.decode(ctxSkip[
    (averageDGeom > 0)
    + (averageDGeom > 2)
    + (averageDGeom > 8)
    ]);
}

//============================================================================
// AttributeDecoderIntf

AttributeDecoderIntf::~AttributeDecoderIntf() = default;

//============================================================================
// AttributeDecoder factory

std::unique_ptr<AttributeDecoderIntf>
makeAttributeDecoder()
{
  return std::unique_ptr<AttributeDecoder>(new AttributeDecoder());
}

//============================================================================
// AttributeDecoder Members

void
AttributeDecoder::decodeSlabBlock(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  const GeometryBrickHeader& gbh,
  AttributeBrickHeader& abh,
  int geom_num_points_minus1,
  int minGeomNodeSizeLog2,
  const char* payload,
  size_t payloadLen,
  PCCPointSet3& slabBlockPointCloud,
  AttributeInterPredParams& attrInterPredParams)
{
  PCCPointSet3 tmp;
  int64_t sumDGeom = 0;
  if (attrInterPredParams.enableAttrInterPred
      && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
    // compensatedPointCloud is needed by geometry
    tmp.swap(attrInterPredParams.compensatedPointCloud);
    attrInterPredParams.attributes_mc.clear();
    // copy geometry but not attributes
    //attrInterPredParams.compensatedPointCloud.addRemoveAttributes(false, false);
    attrInterPredParams.compensatedPointCloud.clear();
    attrInterPredParams.compensatedPointCloud.appendPartition(
      slabBlockPointCloud, 0, slabBlockPointCloud.size(), false);
    // allocates attributes
    attrInterPredParams.compensatedPointCloud.addRemoveAttributes(slabBlockPointCloud);
    if (attr_aps.dual_motion_field_flag)
      sumDGeom = attrInterPredParams.decodeMotionAndBuildCompensated(attr_aps, *_pMotionDecoder);
    else
      sumDGeom = attrInterPredParams.buildCompensatedSlabBlock(attr_aps, gps.motion);
  }

  BlockRefBoundaries blockRefBoundaries;
  BlockBoundaries* boundariesInfo = gbh.uniqueSlabBlock ? nullptr
  : &prepareBoundariesInfo(
    attrInterPredParams.getSlabBlockStart(),
    attrInterPredParams.getSlabBlockSize(),
    attr_desc.attr_num_dimensions_minus1 + 1,
    boundaries, blockRefBoundaries);

  // ---- test skip ----
  const bool skipEnabledFlag = attr_aps.slab_block_skip_enabled_flag;

  if (attrInterPredParams.hasLocalMotion() && skipEnabledFlag)
  {
    const int numPoints = slabBlockPointCloud.size();
    if(_pDecoder->decodeSkip(numPoints ? divApprox(sumDGeom, numPoints, 0) : 0)) {
      if (attr_desc.attr_num_dimensions_minus1 == 0) {
        auto attr = &slabBlockPointCloud.getReflectance(0);
        auto attrPredictor = &attrInterPredParams.compensatedPointCloud.getReflectance(0);
        for (int i = 0; i < numPoints; i++) {
          *attr++ = *attrPredictor++;
        }
      } else if (attr_desc.attr_num_dimensions_minus1 == 2) {
        auto attr = &slabBlockPointCloud.getColor(0);
        auto attrPredictor = &attrInterPredParams.compensatedPointCloud.getColor(0);
        for (int i = 0; i < numPoints; i++) {
          *attr++ = *attrPredictor++;
        }
      } else {
        assert(
          attr_desc.attr_num_dimensions_minus1 == 0
          || attr_desc.attr_num_dimensions_minus1 == 2);
      }
      if (attrInterPredParams.enableAttrInterPred
          && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
        tmp.swap(attrInterPredParams.compensatedPointCloud);
      }
      return;
    }
  }

  if (attrInterPredParams.enableAttrInterPred && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
    for (int i = 0; i < attrInterPredParams.compensatedPointCloud.size(); ++i)
      attrInterPredParams.compensatedPointCloud[i] -= attrInterPredParams.getSlabBlockStart();
  }
  for (int i = 0; i < slabBlockPointCloud.size(); ++i)
    slabBlockPointCloud[i] -= attrInterPredParams.getSlabBlockStart();

  if (attr_desc.attr_num_dimensions_minus1 == 0) {
      decodeReflectancesRaht(
        attr_desc, attr_aps, abh, _qpSet, *_pDecoder, slabBlockPointCloud,
        attrInterPredParams,
        point_t({}), // Now, everything is aligned to zero
        attrInterPredParams.getSlabBlockSize()-1,
        boundariesInfo,
        &blockRefBoundaries);
  } else if (attr_desc.attr_num_dimensions_minus1 == 2) {
      decodeColorsRaht(
        attr_desc, attr_aps, abh, _qpSet, *_pDecoder, slabBlockPointCloud,
        attrInterPredParams,
        point_t({}), // Now, everything is aligned to zero
        attrInterPredParams.getSlabBlockSize()-1,
        boundariesInfo,
        &blockRefBoundaries);
  } else {
    assert(
      attr_desc.attr_num_dimensions_minus1 == 0
      || attr_desc.attr_num_dimensions_minus1 == 2);
  }
  if (attrInterPredParams.enableAttrInterPred
      && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
    tmp.swap(attrInterPredParams.compensatedPointCloud);
  }

  for (int i = 0; i < slabBlockPointCloud.size(); ++i)
    slabBlockPointCloud[i] += attrInterPredParams.getSlabBlockStart();
}

//----------------------------------------------------------------------------

void
AttributeDecoder::startDecode(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  const char* payload,
  size_t payloadLen,
  const AttributeContexts& ctxtMem,
  const PredModeContexts& ctxtMemPredMode,
  const MotionEntropy& ctxtMemDualMotion
  )
{
  _qpSet = deriveQpSet(attr_desc, attr_aps, abh);

  _pDecoder.reset(new PCCResidualsDecoder(abh, ctxtMem, ctxtMemPredMode));
  _pDecoder->start(sps, payload, payloadLen);

  _pMotionDecoder.reset(new MotionEntropyDecoder(
    ctxtMemDualMotion, &_pDecoder->arithmeticDecoder));
}

//----------------------------------------------------------------------------

void
AttributeDecoder::finishDecode(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& attr_desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  AttributeContexts& ctxtMem,
  PredModeContexts& ctxtMemPredMode,
  MotionEntropy& ctxtMemDualMotion)
{
  _pDecoder->stop();

  // save the context state for re-use by a future slice if required
  ctxtMem = _pDecoder->getAttrCtx();
  ctxtMemPredMode = _pDecoder->getModeCtx();
  ctxtMemDualMotion = _pMotionDecoder->getCtx();
}

//----------------------------------------------------------------------------

template<const int attribCount>
inline void
decodeRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud,
  const AttributeInterPredParams& attrInterPredParams,
  point_t blockStart,
  point_t blockSizeMinus1,
  RAHT::BlockBoundaries* blockBoundaries,
  RAHT::BlockRefBoundaries* blockRefBoundaries)
{
  const int voxelCount = pointCloud.getPointCount();

  // Morton codes
  std::vector<int64_t> mortonCode;
  std::vector<attr_t> attributes;
  auto indexOrd = sortedPointCloud(attribCount, pointCloud, mortonCode, attributes);
  attributes.resize(voxelCount * attribCount);

  // Entropy decode
  std::vector<int> coefficients(attribCount * voxelCount, 0);
  std::vector<Qps> pointQpOffsets;
  pointQpOffsets.reserve(voxelCount);

  for (auto index : indexOrd) {
    pointQpOffsets.push_back(qpSet.regionQpOffset(pointCloud[index]));
  }

  if (attrInterPredParams.hasLocalMotion()) {
    auto& attributes_mc = attrInterPredParams.attributes_mc;
    sortedPointCloud(
      attribCount, attrInterPredParams.compensatedPointCloud, indexOrd,
      attributes_mc);

    regionAdaptiveHierarchicalInverseTransform(
      aps, desc, abh, qpSet, pointQpOffsets.data(), attribCount,
      voxelCount, mortonCode.data(), attributes.data(), attributes_mc.data(),
      decoder, blockStart, blockSizeMinus1, blockBoundaries, blockRefBoundaries);
  } else {
    regionAdaptiveHierarchicalInverseTransform(
      aps, desc, abh, qpSet, pointQpOffsets.data(), attribCount,
      voxelCount, mortonCode.data(), attributes.data(), nullptr,
      decoder, blockStart, blockSizeMinus1, blockBoundaries, blockRefBoundaries);
  }

  auto attribute = attributes.begin();
  if (attribCount == 3) {
    for (auto index : indexOrd) {
      auto& color = pointCloud.getColor(index);
      color[0] = *attribute++;
      color[1] = *attribute++;
      color[2] = *attribute++;
    }
  } else if (attribCount == 1) {
    for (auto index : indexOrd) {
      auto& refl = pointCloud.getReflectance(index);
      refl = *attribute++;
    }
  }
}

//----------------------------------------------------------------------------

void
AttributeDecoder::decodeReflectancesRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud,
  const AttributeInterPredParams& attrInterPredParams,
  point_t blockStart,
  point_t blockSizeMinus1,
  RAHT::BlockBoundaries* blockBoundaries,
  RAHT::BlockRefBoundaries* blockRefBoundaries)
{
  decodeRaht<1>(
    desc, aps, abh, qpSet, decoder, pointCloud,
    attrInterPredParams, blockStart, blockSizeMinus1, blockBoundaries,
    blockRefBoundaries);
}

void
AttributeDecoder::decodeColorsRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCResidualsDecoder& decoder,
  PCCPointSet3& pointCloud,
  const AttributeInterPredParams& attrInterPredParams,
  point_t blockStart,
  point_t blockSizeMinus1,
  RAHT::BlockBoundaries* blockBoundaries,
  RAHT::BlockRefBoundaries* blockRefBoundaries)
{
  decodeRaht<3>(
    desc, aps, abh,qpSet, decoder, pointCloud,
    attrInterPredParams, blockStart, blockSizeMinus1, blockBoundaries,
    blockRefBoundaries);
}

//============================================================================

} /* namespace pcc */
