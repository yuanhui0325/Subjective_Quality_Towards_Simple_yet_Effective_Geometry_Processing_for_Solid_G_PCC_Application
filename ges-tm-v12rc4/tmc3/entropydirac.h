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

#include "dependencies/schroedinger/schroarith.h"

#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <memory>
#include <vector>
#include <cmath>
#include <cstring>
#include <stdexcept>

#include "PCCMisc.h"

namespace pcc {
namespace dirac {

  //==========================================================================

  const uint16_t diraclut[256] = {
    //LUT corresponds to window = 16 @ p0=0.5 & 256 @ p=1.0
    0,    2,    5,    8,    11,   15,   20,   24,   29,   35,   41,   47,
    53,   60,   67,   74,   82,   89,   97,   106,  114,  123,  132,  141,
    150,  160,  170,  180,  190,  201,  211,  222,  233,  244,  256,  267,
    279,  291,  303,  315,  327,  340,  353,  366,  379,  392,  405,  419,
    433,  447,  461,  475,  489,  504,  518,  533,  548,  563,  578,  593,
    609,  624,  640,  656,  672,  688,  705,  721,  738,  754,  771,  788,
    805,  822,  840,  857,  875,  892,  910,  928,  946,  964,  983,  1001,
    1020, 1038, 1057, 1076, 1095, 1114, 1133, 1153, 1172, 1192, 1211, 1231,
    1251, 1271, 1291, 1311, 1332, 1352, 1373, 1393, 1414, 1435, 1456, 1477,
    1498, 1520, 1541, 1562, 1584, 1606, 1628, 1649, 1671, 1694, 1716, 1738,
    1760, 1783, 1806, 1828, 1851, 1874, 1897, 1920, 1935, 1942, 1949, 1955,
    1961, 1968, 1974, 1980, 1985, 1991, 1996, 2001, 2006, 2011, 2016, 2021,
    2025, 2029, 2033, 2037, 2040, 2044, 2047, 2050, 2053, 2056, 2058, 2061,
    2063, 2065, 2066, 2068, 2069, 2070, 2071, 2072, 2072, 2072, 2072, 2072,
    2072, 2071, 2070, 2069, 2068, 2066, 2065, 2063, 2060, 2058, 2055, 2052,
    2049, 2045, 2042, 2038, 2033, 2029, 2024, 2019, 2013, 2008, 2002, 1996,
    1989, 1982, 1975, 1968, 1960, 1952, 1943, 1934, 1925, 1916, 1906, 1896,
    1885, 1874, 1863, 1851, 1839, 1827, 1814, 1800, 1786, 1772, 1757, 1742,
    1727, 1710, 1694, 1676, 1659, 1640, 1622, 1602, 1582, 1561, 1540, 1518,
    1495, 1471, 1447, 1422, 1396, 1369, 1341, 1312, 1282, 1251, 1219, 1186,
    1151, 1114, 1077, 1037, 995,  952,  906,  857,  805,  750,  690,  625,
    553,  471,  376,  255};

  //==========================================================================

  struct SchroContext {
    uint16_t probability = 0x8000;  // p=0.5

    template<class... Args>
    void reset(Args...)
    {
      probability = 0x8000;
    }

    template <int FPP, unsigned K=2>
    void getEntropy(int64_t h[2]) const
    {
      h[0] = fpEntropyProbaLUT<16, 16>(probability);
      h[1] = fpEntropyProbaLUT<16, 16>((1ULL << 16) - probability);
    }

  };

  //--------------------------------------------------------------------------
  // The approximate (7 bit) probability of a symbol being 1 or 0 according
  // to a context model.

  inline int approxSymbolProbability(int bit, SchroContext& model)
  {
    int p = std::max(1, model.probability >> 9);
    return bit ? 128 - p : p;
  }

  //==========================================================================

  struct SchroContextFixed {
    //uint16_t probability = 0x8000; // p=0.5
  };

  //==========================================================================
  // context definition that automatically binarises an m-ary symbol

  struct SchroMAryContext {
    SchroMAryContext() = default;

    SchroMAryContext(int numsyms) { set_alphabet(numsyms); }

    void set_alphabet(int numsyms);

    int numsyms;
    std::vector<uint16_t> probabilities;
  };

  //==========================================================================

  struct InvalidContext {
    void set_alphabet(int) {}
  };

  //==========================================================================

  class ArithmeticEncoder {
  public:
    ArithmeticEncoder() = default;
    ArithmeticEncoder(size_t bufferSize, std::nullptr_t)
    {
      setBuffer(bufferSize, nullptr);
    }

    //------------------------------------------------------------------------

    size_t getNumBytes() const
    {
      return _bufWr - _buf;
    }

    //------------------------------------------------------------------------
    void getProbabilityLUT(uint16_t* out) const
    {
      for (int i = 0; i < 512; i++)
        out[i] = impl.lut[i];
    }

    //------------------------------------------------------------------------

    void setBuffer(size_t size, uint8_t* buffer)
    {
      _bufSize = size;
      if (buffer) {
        _buf = _bufWr = buffer;
        allocatedBuffer.reset(nullptr);
      } else {
        allocatedBuffer.reset(new uint8_t[size]);
        _buf = _bufWr = allocatedBuffer.get();
      }
    }

    //------------------------------------------------------------------------

    void setBypassBinCodingWithoutProbUpdate(bool bypass_bin_coding_without_prob_update)
    {
      _bypass_bin_coding_without_prob_update = bypass_bin_coding_without_prob_update;
    }

    //------------------------------------------------------------------------

    void start()
    {
      schro_arith_encode_init(&impl, &writeByteCallback, this);
    }

    //------------------------------------------------------------------------

    size_t stop()
    {
      schro_arith_flush(&impl);

      return _bufWr - _buf;
    }

    //------------------------------------------------------------------------
    // fixed point with 16 bits decimal precision

    uint64_t getNumBitsEstimate() const
    {
      return _getNumBitsEstimate(
        _bufWr - _buf, impl.carry, impl.first_byte, impl.cntr, impl.range[1]);
    }

    //------------------------------------------------------------------------

    const char* buffer() { return reinterpret_cast<char*>(_buf); }

    //------------------------------------------------------------------------

    void encode(int bit, SchroContextFixed&) { encode(bit); }

    //------------------------------------------------------------------------

    void encode(int bit)
    {
      if (_bypass_bin_coding_without_prob_update)
        schro_arith_encode_bypass_bit(&impl, bit);
      else {
        uint16_t probability = 0x8000;  // p=0.5
        schro_arith_encode_bit(&impl, &probability, bit);
      }
    }

    //------------------------------------------------------------------------

    void encode(int data, SchroMAryContext& model);

    //------------------------------------------------------------------------

    void encode(int data, InvalidContext& model) { assert(0); }

    void encode(int bit, SchroContext& model)
    {
      schro_arith_encode_bit(&impl, &model.probability, bit);
    }

    //------------------------------------------------------------------------

    void encode(
      int bit,
      const int& offset,
      SchroContext& model,
      uint16_t* obufSingleBound)
    {
      uint16_t* probability = &model.probability;
      uint16_t& lowTh = obufSingleBound[offset + 1];
      uint16_t& HighTh = obufSingleBound[offset];
      if (*probability > HighTh) {
        *probability = HighTh;
        HighTh += (diraclut[255 - ((HighTh) >> 8)] >> 2);
        if (offset > 0 && HighTh > obufSingleBound[offset - 1]) {
          HighTh = obufSingleBound[offset - 1];
        }
      } else if (*probability < lowTh) {
        *probability = lowTh;
        lowTh -= (diraclut[(lowTh) >> 8] >> 2);
        if (offset < 31 && lowTh < obufSingleBound[offset + 2]) {
          lowTh = obufSingleBound[offset + 2];
        }
      }
      schro_arith_encode_bit(&impl, probability, bit);
      return;
    }

    //------------------------------------------------------------------------

    template<class... Types>
    struct _RDO{
      // n.b. the constructor of dataStart and dataEnd is now deferred to avoid
      //      systematic call to the constructors at the creation of an _RDO
      //      object even if the dataStart and dataEnd are not used.
      //      Thus, for an _RDO object created at a given scope, many memory
      //      operations can be avoided when the _RDO object is declared but
      //      finally not used.
      //
      //      If constructor operations are still too expensive, the
      //      declaration of the _RDO object will have to be made at a higher
      //      scope level (outside of a loop for instance),
      //      even if less readable.
      //
      //      The several contexts may also be expensive as they are all
      //      systematically initialized. At some point context initialization
      //      could also be deferred.
      //
      //      It may also be planed to not provide context tables as part of
      //      `args` parameter, but to iteratively save the state of the
      //      each context beeing modified, and only these contexts.

      _RDO(ArithmeticEncoder& encoder, double lambda, Types&... args)
      : encoder(encoder)
      , lambda(lambda)
      , dataIO(args...)
      {}


      // must be called before any arithmetic encoding operation for the first
      // RDO alternative, and before any modification on supervised variables
      void start()
      {
        bestCost = std::numeric_limits<double>::infinity();
        currentAlternativeIdx = bestAlternativeIdx = -1;
        bufferBytesStart = encoder.getNumBytes();
        schro_arith_get_encode_state(&encoder.impl, &arithStateStart);
        if (!dataStart.isConstructed) {
          dataStart.construct(dataIO.read());
        }
        else {
          dataIO.read(dataStart.val);
        }
      }

      // must be called before any arithmetic encoding, and before any
      // modification on supervised variables, for a given RDO alternative
      void startAlternative()
      {
        ++currentAlternativeIdx;
        bufferBytesEnd = encoder.getNumBytes();
        if (currentAlternativeIdx > 0) {
          if (!dataEnd.isConstructed) {
            dataEnd.construct(dataIO.read());
          }
          else {
            dataIO.read(dataEnd.val);
          }
          schro_arith_set_encode_state(&encoder.impl, &arithStateStart);
          dataIO.write(dataStart.val);
        }
      }

      // must be called after all arithmetic encoding, and after all required
      // modifications on supervised variables, for a given RDO alternative
      std::pair<double/*cost*/, bool/*taken*/> finishAlternative(double D=0.)
      {
        bool taken = false;

        SchroArithEncState currentArithStateEnd;
        auto currentBufferBytesEnd = encoder.getNumBytes();
        schro_arith_get_encode_state(&encoder.impl, &currentArithStateEnd);

        double cost = D + lambda * getNumBitsEstimateCurrent(
          currentBufferBytesEnd, currentArithStateEnd) / 65536.0;

        if (cost < bestCost) {
          taken = true;
          auto numBytes = currentBufferBytesEnd - bufferBytesEnd;
          if (currentAlternativeIdx > 0) {
            std::memmove(
              encoder._buf + bufferBytesStart,
              encoder._buf + bufferBytesEnd,
              numBytes
              );
          }
          bestCost = cost;
          bestAlternativeIdx = currentAlternativeIdx;
          bufferBytesEnd = bufferBytesStart + numBytes;
          arithStateEnd = currentArithStateEnd;
        }
        encoder._bufWr = encoder._buf + bufferBytesEnd;

        return std::make_pair(cost, taken);
      }

      // must be called after all alternatives are finished and before any
      // other arithmetic encoding, or modification on supervised variables
      void finish()
      {
        if (bestAlternativeIdx != currentAlternativeIdx) {
          schro_arith_set_encode_state(&encoder.impl, &arithStateEnd);
          dataIO.write(dataEnd.val);
        }
      }

      uint64_t getNumBitsEstimateCurrent(
        size_t currentBufferBytesEnd,
        const SchroArithEncState& currentArithStateEnd
      ) const
      {
        return ArithmeticEncoder::_getNumBitsEstimate(
          currentBufferBytesEnd,
          currentArithStateEnd.carry,
          currentArithStateEnd.first_byte,
          currentArithStateEnd.cntr,
          currentArithStateEnd.range[1])
        - ArithmeticEncoder::_getNumBitsEstimate(
          bufferBytesEnd,
          arithStateStart.carry,
          arithStateStart.first_byte,
          arithStateStart.cntr,
          arithStateStart.range[1]);
      }

      ArithmeticEncoder& encoder;
      _dataIOInterface<Types...> dataIO;

      template <class T> struct DeferredConstructor {
        // union is used to defer construction of val
        union {
          uint8_t _val_constructor_disabler;
          T val;
        };

        bool isConstructed;

        DeferredConstructor()
          : isConstructed(false)
        {}

        DeferredConstructor(const DeferredConstructor& from)
          : isConstructed(from.isConstructed)
        {
          if (isConstructed) {
            new (&val) T(from.val);
          }
        }

        DeferredConstructor(DeferredConstructor&& from)
          : isConstructed(std::move(from.isConstructed))
        {
          if (isConstructed) {
            new (&val) T(std::move(from.val));
          }
        }

        DeferredConstructor& operator =(const DeferredConstructor& from)
        {
          if (isConstructed) {
            if (from.isConstructed)
              val = from.val;
            else {
              // TODO: can we avoid this else section ?
              isConstructed = destruct();
            }
          } else if(from.isConstructed) {
            new (&val) T(from.val);
          }

          return *this;
        }

        DeferredConstructor& operator =(const DeferredConstructor&& from)
        {
          if (isConstructed) {
            if (from.isConstructed)
              val = std::move(from.val);
            else {
              // TODO: can we avoid this else section ?
              isConstructed = destruct();
            }
          } else if(from.isConstructed) {
            new (&val) T(std::move(from.val));
          }

          return *this;
        }

        ~DeferredConstructor() { if (isConstructed) destruct(); }

        void construct() { isConstructed = new (&val) T; }
        void construct(T&& from) { isConstructed = new (&val) T (from); }
      private:
        bool destruct() { val.~T(); return false; }
      };

      DeferredConstructor<std::tuple<Types...>> dataStart;
      DeferredConstructor<std::tuple<Types...>> dataEnd;

      double lambda;
      double bestCost;
      int bestAlternativeIdx;
      int currentAlternativeIdx;
      size_t bufferBytesStart;
      size_t bufferBytesEnd;
      SchroArithEncState arithStateStart;
      SchroArithEncState arithStateEnd;
    };

    template<class... Types>
    inline _RDO<Types...> makeRDO(double lambda, Types&... args)
    { return _RDO<Types...>(*this, lambda, args...); }

    template<class... Types>
    struct _RateEstimator{
      _RateEstimator(ArithmeticEncoder& encoder, Types&... args)
      : encoder(encoder)
      , dataIO(args...)
      {}

      // TODO: use a single function with test performed within a lambda
      // function passed as parameter (instead of requiring start/finish calls)

      // must be called before any arithmetic encoding operation specific to
      // the rate estimation, and before any specific modification on
      // supervised variables
      void start()
      {
        schro_arith_get_encode_state(&encoder.impl, &arithStateStart);
        dataIO.read(dataStart);
        bufferBytesStart = encoder.getNumBytes();
      }

      // must be after all arithmetic encoding operation specific to
      // the rate estimation, and before any non related modification on
      // supervised variables, or arithmetic encoding operation.
      uint64_t finish()
      {
        uint64_t res = encoder.getNumBitsEstimate()
        - ArithmeticEncoder::_getNumBitsEstimate(
          bufferBytesStart,
          arithStateStart.carry,
          arithStateStart.first_byte,
          arithStateStart.cntr,
          arithStateStart.range[1]);
        schro_arith_set_encode_state(&encoder.impl, &arithStateStart);
        dataIO.write(dataStart);
        encoder._bufWr = encoder._buf + bufferBytesStart;
        return res;
      }

      ArithmeticEncoder& encoder;
      SchroArithEncState arithStateStart;
      _dataIOInterface<Types...> dataIO;
      std::tuple<Types...> dataStart;
      size_t bufferBytesStart;
    };

    template<class... Types>
    inline _RateEstimator<Types...> makeRateEstimator(Types&... args)
    { return _RateEstimator<Types...>(*this, args...); }

    template<class Fct, class... Types>
    inline int64_t getRateEstimate(Fct fct, Types&... args)
    {
      SchroArithEncState arithStateStart;
      schro_arith_get_encode_state(&impl, &arithStateStart);
      std::tuple<Types...> dataStart(args...);
      std::tuple<Types&...> dataRef(args...);
      auto bufferBytesStart = getNumBytes();
      int64_t rateBegin = getNumBitsEstimate();
      fct();
      auto bufferBytesEnd = getNumBytes();
      int64_t rateEnd = getNumBitsEstimate();
      schro_arith_set_encode_state(&impl, &arithStateStart);
      dataRef = dataStart;
      _bufSize = _bufSize + bufferBytesEnd - bufferBytesStart;
      _bufWr = _buf + bufferBytesStart;
      return rateEnd - rateBegin;
    }

    //------------------------------------------------------------------------
  private:
    static uint64_t _getNumBitsEstimate(
      uint64_t buffSize, int carry, bool first_byte, int cntr, uint32_t range)
    {
      int decimalPrec = range & 0x0F;
      int rangeIdx = range - 0x4000 >> 4;
      return ((buffSize + carry + !first_byte << 3) + cntr << 16)
        + _rangeBits[rangeIdx]
        + ((_rangeBits[rangeIdx+1] - _rangeBits[rangeIdx]) * decimalPrec + 8 >> 4);
    }

    //------------------------------------------------------------------------

    static void writeByteCallback(uint8_t byte, void* thisptr)
    {
      auto _this = reinterpret_cast<ArithmeticEncoder*>(thisptr);
      if (_this->_bufSize == 0) {
        if (!_this->allocatedBuffer)
          throw std::runtime_error("Aec stream overflow");
        _this->realloc();
      }
      _this->_bufSize--;
      *_this->_bufWr++ = byte;
    }

    //------------------------------------------------------------------------

    void realloc(uint8_t** new_buffer = nullptr, size_t* new_size = nullptr)
    {
      size_t currentSize = _bufSize + (_bufWr - _buf);
      if (currentSize > 0x10000000)
        // Never use more than 512MB for a bitstream chunk
        throw std::runtime_error("Aec stream overflow");
      size_t newSize = currentSize << 1;
      uint8_t* newBuf = new uint8_t[newSize];
      std::copy(_buf, _buf + currentSize, newBuf);
      allocatedBuffer.reset(newBuf);
      _bufSize += currentSize;
      _bufWr = allocatedBuffer.get() + (_bufWr - _buf);
      _buf = allocatedBuffer.get();
      if (new_buffer) *new_buffer = newBuf;
      if (new_size) *new_size = newSize;
    }

    static void reallocCallback(void* thisptr, uint8_t** new_buffer, size_t* new_size)
    {
      auto _this = reinterpret_cast<ArithmeticEncoder*>(thisptr);
      _this->realloc(new_buffer, new_size);
    }

  private:
    ::SchroArith impl;
    uint8_t* _buf;
    uint8_t* _bufWr;
    size_t _bufSize;
    std::unique_ptr<uint8_t[]> allocatedBuffer;

    // Controls separate coding of bypass bins
    bool _bypass_bin_coding_without_prob_update = false;

    static const uint16_t _rangeBits[1026];
  };

  //==========================================================================

  class ArithmeticDecoder {
  public:
    void setBuffer(size_t size, const char* buffer)
    {
      _buffer = reinterpret_cast<const uint8_t*>(buffer);
      _bufferLen = size;
    }

    //------------------------------------------------------------------------

    void setBypassBinCodingWithoutProbUpdate(bool bypass_bin_coding_without_prob_update)
    {
      _bypass_bin_coding_without_prob_update = bypass_bin_coding_without_prob_update;
    }

    //------------------------------------------------------------------------

    void start()
    {
      schro_arith_decode_init(&impl, &readByteCallback, this);
    }

    //------------------------------------------------------------------------

    void stop() { schro_arith_decode_flush(&impl); }

    //------------------------------------------------------------------------
    // Terminate the arithmetic decoder, and reinitialise to start decoding
    // the next entropy stream.

    void flushAndRestart()
    {
      stop();
      schro_arith_decode_init(&impl, &readByteCallback, this);
    }

    //------------------------------------------------------------------------

    int decode(SchroContextFixed&) { return decode(); }

    //------------------------------------------------------------------------

    int decode()
    {
      if (_bypass_bin_coding_without_prob_update)
        return schro_arith_decode_bypass_bit(&impl);
      else {
        uint16_t probability = 0x8000;  // p=0.5
        return schro_arith_decode_bit(&impl, &probability);
      }
    }

    //------------------------------------------------------------------------

    int decode(SchroMAryContext& model);

    //------------------------------------------------------------------------

    int decode(InvalidContext& model)
    {
      assert(0);
      return 0;
    }

    //------------------------------------------------------------------------

    int decode(SchroContext& model)
    {
      return schro_arith_decode_bit(&impl, &model.probability);
    }

    //------------------------------------------------------------------------

    int
    decode(const int& offset, SchroContext& model, uint16_t* obufSingleBound)
    {
      uint16_t* probability = &model.probability;
      uint16_t& lowTh = obufSingleBound[offset + 1];
      uint16_t& HighTh = obufSingleBound[offset];
      if (*probability > HighTh) {
        *probability = HighTh;
        HighTh += (diraclut[255 - ((HighTh) >> 8)] >> 2);
        if (offset > 0 && HighTh > obufSingleBound[offset - 1]) {
          HighTh = obufSingleBound[offset - 1];
        }
      } else if (*probability < lowTh) {
        *probability = lowTh;
        lowTh -= (diraclut[(lowTh) >> 8] >> 2);
        if (offset < 31 && lowTh < obufSingleBound[offset + 2]) {
          lowTh = obufSingleBound[offset + 2];
        }
      }

      return schro_arith_decode_bit(&impl, probability);
    }

    //------------------------------------------------------------------------
  private:
    static uint8_t readByteCallback(void* thisptr)
    {
      auto _this = reinterpret_cast<ArithmeticDecoder*>(thisptr);
      if (!_this->_bufferLen)
        return 0xff;
      _this->_bufferLen--;
      return *_this->_buffer++;
    }

    //------------------------------------------------------------------------

  private:
    ::SchroArith impl;

    // the user supplied buffer.
    const uint8_t* _buffer;

    // the length of the user supplied buffer
    size_t _bufferLen;

    // Controls separate coding of bypass bins
    bool _bypass_bin_coding_without_prob_update = false;
  };

  //==========================================================================

}  // namespace dirac

using StaticBitModel = dirac::SchroContextFixed;
using StaticMAryModel = dirac::InvalidContext;
using AdaptiveBitModel = dirac::SchroContext;
using AdaptiveBitModelFast = dirac::SchroContext;
using AdaptiveMAryModel = dirac::SchroMAryContext;

//============================================================================

}  // namespace pcc
