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

#include <stdint.h>

#include "Attribute.h"
#include "AttributeCommon.h"
#include "PayloadBuffer.h"
#include "PCCTMC3Common.h"
#include "quantization.h"

namespace pcc {

//============================================================================
// Opaque definitions (Internal detail)
//============================================================================
// An encapsulation of the entropy decoding methods used in attribute coding

  class PCCResidualsDecoder
    : protected AttributeContexts
    , protected PredModeContexts
  {
  public:
    PCCResidualsDecoder(
      const AttributeBrickHeader& abh, const AttributeContexts& ctxtMem, const PredModeContexts& ctxtMemPredMode);

    EntropyDecoder arithmeticDecoder;

    const AttributeContexts& getAttrCtx() const { return *this; }
    const PredModeContexts& getModeCtx() const { return *this; }

    void start(const SequenceParameterSet& sps, const char* buf, int buf_len);
    void stop();

    int decodeZeroBlock(int coeffCnt, int numCoeffNot0, int numCoeffTotal, bool zeroParent, bool enableAveragePrediction);
    int decodeZeroCoeffs(int numCoeffNot0, int numCoeffTotal, int c, int existsNoZeroInBlock, bool enableAveragePrediction);
    int decodeCoeffMag(int k1, int k2, int k3, bool enableAveragePrediction, int ctx0);
    void decode(int32_t values[3], bool enableAveragePrediction, int numCoeffNot0, int numCoeffTotal, bool is420);
    int32_t decode();
    bool decodeSkip(int64_t averageDGeom);

    attr::PredMode decodePredMode(int ctxMode)
    {
      return arithmeticDecoder.decode(modeIsIntra[ctxMode])
        ? attr::PredMode::Intra : attr::PredMode::Inter;
    }
  };

//============================================================================

class AttributeDecoder : public AttributeDecoderIntf {
public:
  void decodeSlabBlock(
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
    AttributeInterPredParams& attrInterPredParams
  ) override;

  void startDecode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh,
    const char* payload,
    size_t payloadLen,
    const AttributeContexts& ctxtMem,
    const PredModeContexts& ctxtMemPredMode,
    const MotionEntropy& ctxtMemDualMotion
  ) override;

  void finishDecode(
    const SequenceParameterSet& sps,
    const GeometryParameterSet& gps,
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    const AttributeBrickHeader& abh,
    AttributeContexts& ctxtMem,
    PredModeContexts& ctxtMemPredMode,
    MotionEntropy& ctxtMemDualMotion
  ) override;

  MotionEntropyDecoder* getMotionDecoder() override { return _pMotionDecoder.get(); }
protected:
  // todo(df): consider alternative encapsulation

  void decodeReflectancesRaht(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    AttributeBrickHeader& abh,
    const QpSet& qpSet,
    PCCResidualsDecoder& decoder,
    PCCPointSet3& pointCloud,
    const AttributeInterPredParams& attrInterPredParams,
    point_t blockStart = 0,
    point_t blockSizeMinus1 = 0,
    RAHT::BlockBoundaries* blockBoundaries = nullptr,
    RAHT::BlockRefBoundaries* blockRefBoundaries = nullptr);

  void decodeColorsRaht(
    const AttributeDescription& desc,
    const AttributeParameterSet& aps,
    AttributeBrickHeader& abh,
    const QpSet& qpSet,
    PCCResidualsDecoder& decoder,
    PCCPointSet3& pointCloud,
    const AttributeInterPredParams& attrInterPredParams,
    point_t blockStart = 0,
    point_t blockSizeMinus1 = 0,
    RAHT::BlockBoundaries* blockBoundaries = nullptr,
    RAHT::BlockRefBoundaries* blockRefBoundaries = nullptr);

private:
  // for local attributes
  std::unique_ptr<PCCResidualsDecoder> _pDecoder;
  // for dual Motion
  std::unique_ptr<MotionEntropyDecoder> _pMotionDecoder;

  QpSet _qpSet;

  ringbuf<std::pair<point_t/*slab block start*/,RAHT::BlockBoundaries>> boundaries;
};

//============================================================================

} /* namespace pcc */
