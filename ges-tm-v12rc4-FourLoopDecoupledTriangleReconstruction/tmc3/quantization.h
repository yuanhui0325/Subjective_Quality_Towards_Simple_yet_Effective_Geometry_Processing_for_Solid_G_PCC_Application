/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2019, ISO/IEC
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

#include <array>
#include <vector>

#include "constants.h"
#include "PCCMath.h"

namespace pcc {

struct AttributeDescription;
struct AttributeParameterSet;
struct AttributeBrickHeader;

//============================================================================
// Quantisation methods

class Quantizer {
public:
  Quantizer() = default;
  // Derives step sizes from qp
  Quantizer(int qp, int internalBitsPrecision = 0);
  Quantizer(const Quantizer&) = default;
  Quantizer& operator=(const Quantizer&) = default;

  // The quantizer's step size
  int stepSize() const { return _stepSize; }

  // Quantise a value
  int64_t quantize(int64_t x, int typeNode, int W) const;

  // Scale (inverse quantise) a quantised value
  int64_t scale(int64_t x, int typeNode, int W) const;
  int64_t getStepSize() const;

private:
  // Quantisation step size
  int _stepSize;

  // Reciprocal stepsize for forward quantisation optimisation
  int _stepSizeRecip;
  // adpative quantization
  static const int64_t LUTDZ[2][4];
  static const int64_t LUTdeQinf[4];
};

//---------------------------------------------------------------------------

inline int64_t
Quantizer::quantize(int64_t x, int typeNode, int W) const
{
  // Forward quantisation avoids division by using the multiplicative inverse
  // with 22 fractional bits.
  const int64_t fracBits = 22 + 6 + kFixedPointAttributeShift;

  // NB, the folowing offsets quantizes with a different deadzone to the
  // reconstruction function.
  const int64_t DZ = LUTDZ[W > 64][typeNode];
  const int64_t offset = (1ll << fracBits) * (256 - DZ) >> 8;

  return (x >= 0) * ((x * _stepSizeRecip + offset) >> fracBits)
    - (x < 0) * ((offset - x * _stepSizeRecip) >> fracBits);
}

//---------------------------------------------------------------------------

inline int64_t
Quantizer::scale(int64_t x, int typeNode, int W) const
{
  const int64_t DZ = LUTDZ[W > 64][typeNode];
  int64_t deQ = LUTdeQinf[typeNode];
  return fpReduce<kFixedPointAttributeShift + 8>((256 * x + ((x >= 1) - (x <= -1)) * (DZ - 256 + deQ)) * _stepSize);
}

//---------------------------------------------------------------------------
inline int64_t
Quantizer::getStepSize() const
{
  return _stepSize;
}

//============================================================================
// Encapslation of multi-component attribute quantizer values.

typedef std::array<Quantizer, 2> Quantizers;
typedef std::array<int, 2> Qps;

//============================================================================

struct QpRegionOffset {
  Qps qpOffset;
  Box3<int32_t> region;
};

typedef std::vector<QpRegionOffset> QpRegionList;

//============================================================================

struct QpSet {
  Qps baseQp;
  QpRegionList regions;
  int maxQp;
  int attr_frac_bits;

  // Derive the quantizers at a given layer after applying qpOffset
  Quantizers quantizers(Qps qpOffset) const;

  // Derive the quantizer for a point at a particular layer
  Quantizers quantizers(const Vec3<int32_t>& point) const;

  Qps regionQpOffset(const Vec3<int32_t>& point) const;
};

//============================================================================
// Determine the Qps for a particular layer in an attribute slice
Qps deriveQps(
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh);

// Determine a list of Qp offsets per region
QpRegionList deriveQpRegions(
  const AttributeParameterSet& attr_aps, const AttributeBrickHeader& abh);

// Determine the Qp configuration for an attribute slice
QpSet deriveQpSet(
  const AttributeDescription& attrDesc,
  const AttributeParameterSet& attr_aps,
  const AttributeBrickHeader& abh);


//============================================================================

}  // namespace pcc
