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

#include "AttributeEncoder.h"

#include "constants.h"
#include "entropy.h"
#include "io_hls.h"
#include "quantization.h"
#include "RAHT.h"
#include "PCCTMC3Encoder.h"

#include <algorithm>

// todo(df): promote to per-attribute encoder parameter
//static const double kAttrPredLambdaR = 0.01;
//static const double kAttrPredLambdaC = 0.14;

namespace pcc {

using namespace RAHT;

//----------------------------------------------------------------------------

PCCResidualsEncoder::PCCResidualsEncoder(
  const AttributeParameterSet& aps,
  const AttributeBrickHeader& abh,
  const AttributeContexts& ctxtMem,
  const PredModeContexts& ctxtMemPredMode)
  : AttributeContexts(ctxtMem)
  , PredModeContexts(ctxtMemPredMode)
{
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::start(const SequenceParameterSet& sps, int pointCount)
{
  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = pointCount * 3 * 2 + 1024;
  arithmeticEncoder.setBuffer(maxAcBufLen, nullptr);
  arithmeticEncoder.setBypassBinCodingWithoutProbUpdate(
    sps.bypass_bin_coding_without_prob_update);
  arithmeticEncoder.start();
}

//----------------------------------------------------------------------------

int
PCCResidualsEncoder::stop()
{
  return arithmeticEncoder.stop();
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeZeroBlock(bool flag, int coeffCnt, int numCoeffNot0, int numCoeffTotal, bool zeroParent, bool enableAveragePrediction)
{
  int ctx = (numCoeffNot0 > 0) + (numCoeffNot0 > 1) + (numCoeffNot0 > 3);
  ctx = 2 * ctx + zeroParent;
  ctx = 3 * ctx + (4 * numCoeffNot0 >= numCoeffTotal) + (2 * numCoeffNot0 >= numCoeffTotal);
  ctx = 3 * ctx + (coeffCnt > 2) + (coeffCnt > 5);

  arithmeticEncoder.encode(flag, ctxZeroBlock[ctx][enableAveragePrediction]);
}



//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeZeroCoeffs(bool flag, int numCoeffNot0, int numCoeffTotal, int blockIndex, int existsNoZeroInBlock, bool enableAveragePrediction)
{
  static constexpr std::array<int, 8> blockIndexMap {0, 2, 1, 3, 0, 3, 3, 3};

  int ctx = (numCoeffNot0 > 0) + (numCoeffNot0 > 1) + (numCoeffNot0 > 3);
  ctx = 3 * ctx + (4 * numCoeffNot0 >= numCoeffTotal) + (2 * numCoeffNot0 >= numCoeffTotal);

  ctx = 4 * ctx + blockIndexMap[blockIndex];
  ctx = 3 * ctx + existsNoZeroInBlock;

  arithmeticEncoder.encode(flag, ctxZeroCoeffs[ctx][enableAveragePrediction]);
}


//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeCoeffMag(uint32_t value, int k1, int k2, int k3, bool enableAveragePrediction, int ctx0)
{
  int ctx = enableAveragePrediction + 2 * (k3 == 0) * ctx0;
  arithmeticEncoder.encode(value > 0, ctxCoeffGtN[0][k1][ctx]);
  if (!value)
    return;

  arithmeticEncoder.encode(--value > 0, ctxCoeffGtN[1][k2][enableAveragePrediction]);
  if (!value)
    return;

  arithmeticEncoder.encodeExpGolomb(
    --value, 1, ctxCoeffRemPrefix[k3], ctxCoeffRemSuffix[k3]);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(int32_t value0, int32_t value1, int32_t value2, bool enableAveragePrediction, int numCoeffNot0, int numCoeffTotal, bool is420)
{
  int mag0 = abs(value0);
  int mag1 = abs(value1);
  int mag2 = abs(value2);

  int b0 = (mag1 == 0);
  int b1 = (mag1 <= 1);
  int b2 = (mag2 == 0);
  int b3 = (mag2 <= 1);
  int ctx0 = (4 * numCoeffNot0 >= numCoeffTotal) + (2 * numCoeffNot0 >= numCoeffTotal);
  if (!is420) {
    encodeCoeffMag(mag1, 0, 0, 1, enableAveragePrediction, ctx0);
    encodeCoeffMag(mag2, 1 + b0, 1 + b1, 1, enableAveragePrediction, ctx0);
  }
  auto mag0minusX = b0 && b2 ? mag0 - 1 : mag0;
  assert(mag0minusX >= 0);
  encodeCoeffMag(mag0minusX, 3 + (b0 << 1) + b2, 3 + (b1 << 1) + b3, 0, enableAveragePrediction, ctx0);

  if (mag0)
    arithmeticEncoder.encode(value0 < 0);
  if (mag1)
    arithmeticEncoder.encode(value1 < 0);
  if (mag2)
    arithmeticEncoder.encode(value2 < 0);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encodeSkip(bool flag, int64_t averageDGeom)
{
  arithmeticEncoder.encode(flag, ctxSkip[
    (averageDGeom > 0)
    + (averageDGeom > 2)
    + (averageDGeom > 8)
  ]);
}

//----------------------------------------------------------------------------

void
PCCResidualsEncoder::encode(int32_t value)
{
  int mag = abs(value) - 1;
  encodeCoeffMag(mag, 0, 0, 0, false, 0);
  arithmeticEncoder.encode(value < 0);
}

//============================================================================
// AttributeEncoderIntf

AttributeEncoderIntf::~AttributeEncoderIntf() = default;

//============================================================================
// AttributeEncoder factory

std::unique_ptr<AttributeEncoderIntf>
makeAttributeEncoder()
{
  return std::unique_ptr<AttributeEncoder>(new AttributeEncoder());
}

//============================================================================
// AttributeEncoder Members

void AttributeEncoder::encodeSlabBlock(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& slabBlockPointCloud,
  PayloadBuffer* payload,
  const EncoderAttributeParams& attrEncParams,
  AttributeInterPredParams& attrInterPredParams)
{
  PCCPointSet3 tmp;
  PCCPointSet3 tmpOrig;
  int64_t sumDGeom = 0;
  // local motion encoding and compensation performed by slab
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
      sumDGeom = attrInterPredParams.encodeMotionAndBuildCompensated(
        attr_aps, *_pMotionEncoder, _fp16SizeMotionBits);
    else
      sumDGeom = attrInterPredParams.buildCompensatedSlabBlock(attr_aps, gps.motion);
  }

  BlockRefBoundaries blockRefBoundaries;
  BlockBoundaries* boundariesInfo = gbh.uniqueSlabBlock ? nullptr
  : &prepareBoundariesInfo(
    attrInterPredParams.getSlabBlockStart(),
    attrInterPredParams.getSlabBlockSize(),
    desc.attr_num_dimensions_minus1 + 1,
    boundaries, blockRefBoundaries);

  auto quantizers = _qpSet.quantizers(0);
  // TODO: shouldn't stepQ be scaled according to fractional bits or internal bitdepth ?
  double stepQ = quantizers[0].stepSize() / 65536.;
  double lambda = stepQ * stepQ * attrEncParams.slabBlockSkipRDO_lambdaFactor;

  BlockBoundaries _tmp0;
  auto RDO = _pEncoder->arithmeticEncoder.makeRDO(
    lambda, _pEncoder->getAttrCtx(), _pEncoder->getModeCtx(),
    boundariesInfo ? *boundariesInfo : _tmp0);

  bool skip_taken = false;
  // ---- test skip ----
  const bool skipEnabledFlag = attr_aps.slab_block_skip_enabled_flag;
  const bool fastSkip = attrEncParams.slabBlockSkipFastRDO;
  int64_t errorBlock_skip = 0;
  int64_t averageDGeom = 0;
  if (attrInterPredParams.hasLocalMotion() && skipEnabledFlag)
  {
    const int numPoints = slabBlockPointCloud.size();
    averageDGeom = numPoints ? divApprox(sumDGeom, numPoints, 0) : 0;
    tmpOrig = slabBlockPointCloud; // Save for distortion computation in RDO

    if (desc.attr_num_dimensions_minus1 == 0) {
      auto attr = &slabBlockPointCloud.getReflectance(0);
      auto attrPredictor = &attrInterPredParams.compensatedPointCloud.getReflectance(0);
      for (int i = 0; i < numPoints; i++) {
        for (int k = 0; k < 1; k++) {
          int64_t temp = ((*attr++) - (*attrPredictor++));
          errorBlock_skip += temp * temp;
        }
      }
    } else if (desc.attr_num_dimensions_minus1 == 2) {
      auto attr = &slabBlockPointCloud.getColor(0);
      auto attrPredictor = &attrInterPredParams.compensatedPointCloud.getColor(0);
      for (int i = 0; i < numPoints; i++) {
        for (int k = 0; k < 3; k++) {
          int64_t temp = ((*attr)[k] - (*attrPredictor)[k]);
          errorBlock_skip += temp * temp;
        }
        attr++;
        attrPredictor++;
      }
    } else {
      assert(
        desc.attr_num_dimensions_minus1 == 0
        || desc.attr_num_dimensions_minus1 == 2);
    }

    if (fastSkip) {
      double meanError = errorBlock_skip / double(numPoints * (desc.attr_num_dimensions_minus1 + 1));
      double mseFactor = attrEncParams.slabBlockSkipFastRDO_mseFactor;

      skip_taken = meanError < mseFactor * stepQ * stepQ;

      _pEncoder->encodeSkip(skip_taken, averageDGeom);
    }
    else {
      RDO.start();
      // --- test skip ---
      RDO.startAlternative();
      _pEncoder->encodeSkip(true, averageDGeom);
      RDO.finishAlternative(errorBlock_skip).first;

      // --- start test normal ---
      RDO.startAlternative();
      _pEncoder->encodeSkip(false, averageDGeom);
    }
  }

  if (!skip_taken) {
    if (attrInterPredParams.enableAttrInterPred && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
      for (int i = 0; i < attrInterPredParams.compensatedPointCloud.size(); ++i)
        attrInterPredParams.compensatedPointCloud[i] -= attrInterPredParams.getSlabBlockStart();
    }
    for (int i = 0; i < slabBlockPointCloud.size(); ++i)
      slabBlockPointCloud[i] -= attrInterPredParams.getSlabBlockStart();

    if (desc.attr_num_dimensions_minus1 == 0) {
        encodeReflectancesTransformRaht(
          desc, attr_aps, _abh, _qpSet, slabBlockPointCloud, *_pEncoder,
          attrEncParams, attrInterPredParams,
          point_t({}), // Now, everything is aligned to zero
          attrInterPredParams.getSlabBlockSize()-1,
          boundariesInfo,
          &blockRefBoundaries);
    } else if (desc.attr_num_dimensions_minus1 == 2) {
        encodeColorsTransformRaht(
          desc, attr_aps, _abh, _qpSet, slabBlockPointCloud, *_pEncoder,
          attrEncParams, attrInterPredParams,
          point_t({}), // Now, everything is aligned to zero
          attrInterPredParams.getSlabBlockSize()-1,
          boundariesInfo,
          &blockRefBoundaries);
    } else {
      assert(
        desc.attr_num_dimensions_minus1 == 0
        || desc.attr_num_dimensions_minus1 == 2);
    }

    if (!fastSkip) {
      int64_t errorBlock_noSkip = 0;
      if (attrInterPredParams.hasLocalMotion() && skipEnabledFlag)
      {
        const int numPoints = slabBlockPointCloud.size();
        if (desc.attr_num_dimensions_minus1 == 0) {
          auto attr = &slabBlockPointCloud.getReflectance(0);
          auto attrOrig = &tmpOrig.getReflectance(0);
          for (int i = 0; i < numPoints; i++) {
            for (int k = 0; k < 1; k++) {
              int64_t temp = ((*attr++) - (*attrOrig++));
              errorBlock_noSkip += temp * temp;
            }
          }
        } else if (desc.attr_num_dimensions_minus1 == 2) {
          auto attr = &slabBlockPointCloud.getColor(0);
          auto attrOrig = &tmpOrig.getColor(0);
          for (int i = 0; i < numPoints; i++) {
            for (int k = 0; k < 3; k++) {
              int64_t temp = ((*attr)[k] - (*attrOrig)[k]);
              errorBlock_noSkip += temp * temp;
            }
            attr++;
            attrOrig++;
          }
        }

        auto rdo_res = RDO.finishAlternative(errorBlock_noSkip);
        skip_taken = !rdo_res.second;
        RDO.finish();
      }
    }
  }

  if (skip_taken) {
    const int numPoints = slabBlockPointCloud.size();

    if (desc.attr_num_dimensions_minus1 == 0) {
      auto attr = &slabBlockPointCloud.getReflectance(0);
      auto attrPredictor = &attrInterPredParams.compensatedPointCloud.getReflectance(0);
      for (int i = 0; i < numPoints; i++) {
        *attr++ = *attrPredictor++;
      }
    } else if (desc.attr_num_dimensions_minus1 == 2) {
      auto attr = &slabBlockPointCloud.getColor(0);
      auto attrPredictor = &attrInterPredParams.compensatedPointCloud.getColor(0);
      for (int i = 0; i < numPoints; i++) {
        *attr++ = *attrPredictor++;
      }
    }
  }

  if (attrInterPredParams.enableAttrInterPred
      && !attrInterPredParams.mSOctreeRef.nodes.empty()) {
    tmp.swap(attrInterPredParams.compensatedPointCloud);
  }

  if (!fastSkip || !skip_taken) {
    for (int i = 0; i < slabBlockPointCloud.size(); ++i)
      slabBlockPointCloud[i] += attrInterPredParams.getSlabBlockStart();
  }
}

//----------------------------------------------------------------------------

void AttributeEncoder::startEncode(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh,
  const AttributeContexts& ctxtMem,
  const PredModeContexts& ctxtMemPredMode,
  const MotionEntropy& ctxtMemDualMotion,
  uint32_t pointCountInPointCloud)
{
  _abh = abh;

  _pEncoder.reset(new PCCResidualsEncoder(attr_aps, _abh, ctxtMem, ctxtMemPredMode));
  _pEncoder->start(sps, int(pointCountInPointCloud));

  _pMotionEncoder.reset(new MotionEntropyEncoder(
    ctxtMemDualMotion, &_pEncoder->arithmeticEncoder));

  _qpSet = deriveQpSet(desc, attr_aps, _abh);

  _fp16SizeMotionBits = 0;
}

//----------------------------------------------------------------------------

void AttributeEncoder::finishEncode(
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const AttributeDescription& desc,
  const AttributeParameterSet& attr_aps,
  AttributeContexts& ctxtMem,
  PredModeContexts& ctxtMemPredMode,
  MotionEntropy& ctxtMemDualMotion,
  PayloadBuffer* payload)
{
  uint32_t acDataLen = _pEncoder->stop();

  // write abh
  write(sps, attr_aps, _abh, payload);
  _abh = AttributeBrickHeader();

  std::copy_n(
    _pEncoder->arithmeticEncoder.buffer(), acDataLen,
    std::back_inserter(*payload));

  // save the context state for re-use by a future slice if required
  ctxtMem = _pEncoder->getAttrCtx();
  ctxtMemPredMode = _pEncoder->getModeCtx();
  ctxtMemDualMotion = _pMotionEncoder->getCtx();
}

//----------------------------------------------------------------------------
// may be put directly in RAHT transform ?
template<const int attribCount>
void
mapAttributesToPointCloud(
  const int voxelCount,
  const int* indexOrd,
  const attr_t* attributes,
  PCCPointSet3& pointCloud)
{
  auto attribute = attributes;
  if (attribCount == 3) {
    for (int n = 0; n < voxelCount; ++n) {
      auto index = indexOrd[n];
      auto& color = pointCloud.getColor(index);
      color[0] = *attribute++;
      color[1] = *attribute++;
      color[2] = *attribute++;
    }
  } else if (attribCount == 1) {
    for (int n = 0; n < voxelCount; ++n) {
      auto index = indexOrd[n];
      auto& refl = pointCloud.getReflectance(index);
      refl = *attribute++;
    }
  }
}

//----------------------------------------------------------------------------

template<const int attribCount>
inline void
encodeRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder,
  const EncoderAttributeParams& attrEncParams,
  const AttributeInterPredParams& attrInterPredParams,
  point_t blockStart,
  point_t blockSizeMinus1,
  RAHT::BlockBoundaries* blockBoundaries,
  RAHT::BlockRefBoundaries* blockRefBoundaries)
{
  const int voxelCount = pointCloud.getPointCount();

  // Allocate arrays.
  std::vector<int64_t> mortonCode;
  std::vector<attr_t> attributes;
  std::vector<Qps> pointQpOffsets;

  // Populate input arrays.
  auto indexOrd =
    sortedPointCloud(attribCount, pointCloud, mortonCode, attributes, true);
  pointQpOffsets.reserve(voxelCount);
  for (auto index : indexOrd) {
    pointQpOffsets.push_back(qpSet.regionQpOffset(pointCloud[index]));
  }
  if (attrInterPredParams.hasLocalMotion()) {
    auto& attributes_mc = attrInterPredParams.attributes_mc;
    // Allocate arrays.
    sortedPointCloud(
      attribCount, attrInterPredParams.compensatedPointCloud, indexOrd,
      attributes_mc);

    // Transform.
    regionAdaptiveHierarchicalTransformRSO(
      aps, desc, attrEncParams, abh, qpSet, pointQpOffsets.data(), attribCount,
      voxelCount, mortonCode.data(), attributes.data(), attributes_mc.data(),
      encoder, blockStart, blockSizeMinus1, blockBoundaries, blockRefBoundaries);
  } else {
    // Transform.
    regionAdaptiveHierarchicalTransformRSO(
      aps, desc, attrEncParams, abh, qpSet, pointQpOffsets.data(), attribCount,
      voxelCount, mortonCode.data(), attributes.data(), nullptr,
      encoder, blockStart, blockSizeMinus1, blockBoundaries, blockRefBoundaries);
  }

  mapAttributesToPointCloud<attribCount>(
    voxelCount, indexOrd.data(), attributes.data(), pointCloud);
}

//----------------------------------------------------------------------------

void
AttributeEncoder::encodeReflectancesTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder,
  const EncoderAttributeParams& attrEncParams,
  const AttributeInterPredParams& attrInterPredParams,
  point_t blockStart,
  point_t blockSizeMinus1,
  RAHT::BlockBoundaries* blockBoundaries,
  RAHT::BlockRefBoundaries* blockRefBoundaries)
{
  encodeRaht<1>(
    desc, aps,abh, qpSet, pointCloud, encoder, attrEncParams,
    attrInterPredParams, blockStart, blockSizeMinus1, blockBoundaries,
    blockRefBoundaries);
}

void
AttributeEncoder::encodeColorsTransformRaht(
  const AttributeDescription& desc,
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpSet,
  PCCPointSet3& pointCloud,
  PCCResidualsEncoder& encoder,
  const EncoderAttributeParams& attrEncParams,
  const AttributeInterPredParams& attrInterPredParams,
  point_t blockStart,
  point_t blockSizeMinus1,
  RAHT::BlockBoundaries* blockBoundaries,
  RAHT::BlockRefBoundaries* blockRefBoundaries)
{
  encodeRaht<3>(
    desc, aps, abh, qpSet, pointCloud, encoder, attrEncParams,
    attrInterPredParams, blockStart, blockSizeMinus1, blockBoundaries,
    blockRefBoundaries);
}

//============================================================================

} /* namespace pcc */
