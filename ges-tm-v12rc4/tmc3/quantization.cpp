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

#include "quantization.h"

#include "constants.h"
#include "hls.h"
#include "tables.h"

namespace pcc {

//============================================================================

 // 8 bit precision dead zone half width
const int64_t Quantizer::LUTDZ[2][4] = {
  { 171 , 171 , 171, 171 },  // intra luma, chroma; inter luma, chroma
  { 130 , 150 , 171, 171 }};

// 8 bit precision deq position
const int64_t Quantizer::LUTdeQinf[4] =
  { 128 , 128 , 102, 64 }; // intra luma, chroma; inter luma, chroma

//============================================================================

Quantizer::Quantizer(int qp, int internalBitPrecision)
{
  qp = std::max(qp, 4);
  int qpShift = qp / 6 + internalBitPrecision;
  _stepSize = kQpStep[qp % 6] << qpShift;
  _stepSizeRecip = kQpStepRecip[qp % 6] >> qpShift;
}

//============================================================================

Qps
deriveQps(
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh)
{
  int sliceQpLuma = attr_aps.init_qp_minus4 + 4;
  int sliceQpChroma = attr_aps.aps_chroma_qp_offset;

  if (attr_aps.aps_slice_qp_deltas_present_flag) {
    sliceQpLuma += abh.attr_qp_delta_luma;
    sliceQpChroma += abh.attr_qp_delta_chroma;
  }

  return {sliceQpLuma, sliceQpChroma};
}

//============================================================================

QpRegionList
deriveQpRegions(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh)
{
  QpRegionList regions;
  regions.reserve(abh.qpRegions.size());

  for (int i = 0; i < abh.qpRegions.size(); i++) {
    regions.emplace_back();
    auto& region = regions.back();
    const auto& src = abh.qpRegions[i];

    region.qpOffset = src.attr_region_qp_offset;
    region.region.min = src.regionOrigin;
    region.region.max = src.regionOrigin + src.regionSize;
  }

  return regions;
}

//============================================================================

QpSet
deriveQpSet(
  const AttributeDescription& attrDesc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh)
{
  QpSet qpset;
  qpset.baseQp = deriveQps(attr_aps, abh);
  qpset.regions = deriveQpRegions(attr_aps, abh);

  // The mimimum Qp = 4 is always lossless; the maximum varies according to
  // bitdepth.
  qpset.maxQp = 51 + 6 * (attrDesc.bitdepth - 8);
  qpset.attr_frac_bits = attrDesc.internalBitdepth - attrDesc.bitdepth;

  return qpset;
}

//============================================================================
// Determines the quantizers at a given layer
Quantizers
QpSet::quantizers(Qps qpOffset) const
{
  int qp0 = PCCClip(baseQp[0] + qpOffset[0], 4, maxQp);
  int qp1 = PCCClip(baseQp[1] + qpOffset[1] + qp0, 4, maxQp);

  return {Quantizer(qp0, attr_frac_bits), Quantizer(qp1, attr_frac_bits)};
}

//============================================================================
// Determines the quantizers for a point at a given layer
Quantizers
QpSet::quantizers(const Vec3<int32_t>& point) const
{
  for (const auto& region : regions) {
    if (region.region.contains(point))
      return quantizers(region.qpOffset);
  }

  return quantizers({0, 0});
}

//============================================================================
//for RAHT region QP Offset
Qps
QpSet::regionQpOffset(const Vec3<int32_t>& point) const
{
  for (const auto& region : regions) {
    if (region.region.contains(point))
      return region.qpOffset;
  }

  return {0, 0};
}

//============================================================================

}  // namespace pcc
