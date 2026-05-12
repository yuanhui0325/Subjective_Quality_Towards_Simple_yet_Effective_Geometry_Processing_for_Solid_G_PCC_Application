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

#include "geometry_octree.h"

#include <algorithm>
#include <iterator>
#include <climits>

#include "PCCMisc.h"
#include "geometry_params.h"
#include "quantization.h"
#include "tables.h"

namespace pcc {

//============================================================================

const int CtxModelDynamicOBUF::kContextsInitProbability[] = {
  65461, 65160, 64551, 63637, 62426, 60929, 59163, 57141, 54884, 52413, 49753,
  46929, 43969, 40899, 37750, 34553, 31338, 28135, 24977, 21893, 18914, 16067,
  13382, 10883, 8596,  6542,  4740,  3210,  1967,  1023,  388,   75
};


//============================================================================

void
GeometryOctreeContexts::resetMap(bool forTrisoup)
{
  for (int i = 0; i < 2; i++) {

    const int n2 = 6;
    const int nT2 = forTrisoup ? 2 : 0;
    _MapOccupancy[i][0].reset(6 + n2 + 1 - nT2, 18 - 6 - n2);
    _MapOccupancy[i][1].reset(6 + n2 + 1 - nT2, 18 - 6 - n2);
    _MapOccupancy[i][2].reset(6 + n2 + 1 - nT2, 18 - 6 - n2);
    _MapOccupancy[i][3].reset(4 + n2 + 1 - nT2, 18 - 6 - n2);
    _MapOccupancy[i][4].reset(6 + n2 + 1 - nT2, 18 - 6 - n2);
    _MapOccupancy[i][5].reset(6 + n2 + 1 - nT2, 18 - 6 - n2);
    _MapOccupancy[i][6].reset(6 + n2 + 1 - nT2, 18 - 6 - n2);
    _MapOccupancy[i][7].reset(4 + n2 + 1 - nT2, 18 - 6 - n2);

    const int n3 = 5;
    const int nT3 = forTrisoup ? 4 : 0;
    _MapOccupancySparse[i][0].reset(6 + n3 + 1 - nT3, 9  - n3);
    _MapOccupancySparse[i][1].reset(6 + n3 + 1 - nT3, 12 - n3);
    _MapOccupancySparse[i][2].reset(6 + n3 + 1 - nT3, 12 - n3);
    _MapOccupancySparse[i][3].reset(6 + n3 + 1 - nT3, 11 - n3);
    _MapOccupancySparse[i][4].reset(6 + n3 + 1 - nT3, 9  - n3);
    _MapOccupancySparse[i][5].reset(6 + n3 + 1 - nT3, 12 - n3);
    _MapOccupancySparse[i][6].reset(6 + n3 + 1 - nT3, 12 - n3);
    _MapOccupancySparse[i][7].reset(6 + n3 + 1 - nT3, 11 - n3);
  }

  _OBUFleafNumber = 0;
  std::fill_n(_BufferOBUFleaves, kLeavesBufferSize, 0);

  // octree intra
  static const uint8_t initValueOcc0[64] = {  127, 17, 82, 38, 127, 105, 141, 81, 127, 15, 45, 43, 116, 105, 152, 115, 127, 53, 21, 20, 127, 127, 127, 37, 127, 127, 127, 127, 127, 127, 127, 127, 171, 186, 170, 240, 182, 209, 223, 240, 44, 101, 101, 74, 65, 66, 134, 199, 47, 27, 141, 113, 126, 61, 240, 151, 45, 68, 113, 101, 47, 84, 153, 234, };
  static const uint8_t initValueOcc1[64] = {  240, 240, 222, 240, 175, 181, 127, 127, 120, 152, 132, 116, 57, 127, 127, 127, 105, 185, 127, 87, 105, 116, 65, 69, 66, 105, 58, 43, 44, 49, 18, 15, 228, 240, 138, 240, 178, 198, 114, 152, 173, 240, 204, 127, 70, 141, 127, 127, 184, 192, 105, 116, 121, 181, 35, 46, 58, 87, 114, 73, 51, 15, 101, 40, };
  static const uint8_t initValueOcc2[64] = {  194, 240, 173, 190, 115, 129, 87, 87, 168, 161, 116, 92, 127, 127, 26, 96, 160, 106, 96, 127, 86, 109, 105, 127, 116, 68, 80, 27, 116, 116, 46, 19, 240, 240, 205, 114, 215, 194, 134, 78, 225, 182, 191, 141, 122, 127, 58, 127, 200, 214, 124, 89, 188, 161, 91, 59, 126, 126, 74, 152, 80, 96, 59, 127, };
  static const uint8_t initValueOcc3[64] = {  59, 121, 160, 210, 171, 211, 240, 231, 127, 56, 149, 125, 127, 115, 230, 204, 55, 127, 78, 192, 127, 182, 197, 218, 35, 39, 15, 72, 96, 87, 151, 139, 46, 141, 152, 240, 114, 162, 240, 240, 87, 69, 127, 96, 44, 67, 129, 155, 53, 105, 141, 73, 96, 105, 198, 128, 15, 35, 96, 57, 127, 96, 127, 96, };
  static const uint8_t initValueOcc4[64] = {  23, 30, 130, 66, 139, 127, 30, 105, 113, 127, 87, 127, 127, 127, 127, 127, 166, 146, 70, 15, 209, 116, 141, 90, 114, 138, 71, 15, 127, 127, 127, 127, 204, 240, 198, 219, 232, 240, 142, 240, 151, 139, 87, 127, 209, 190, 43, 141, 141, 181, 116, 127, 240, 210, 88, 127, 73, 170, 65, 61, 140, 194, 48, 65, };
  static const uint8_t initValueOcc5[64] = {  240, 99, 240, 69, 189, 96, 105, 80, 154, 233, 152, 141, 127, 152, 127, 127, 166, 48, 57, 15, 97, 41, 43, 15, 127, 116, 127, 127, 127, 85, 127, 127, 235, 214, 177, 154, 240, 240, 161, 61, 219, 185, 152, 208, 157, 90, 127, 127, 117, 138, 69, 30, 154, 80, 62, 15, 141, 121, 127, 127, 127, 41, 127, 105, };
  static const uint8_t initValueOcc6[64] = {  227, 199, 188, 103, 212, 141, 205, 55, 240, 240, 210, 141, 178, 70, 127, 127, 240, 84, 139, 73, 139, 60, 127, 59, 161, 127, 127, 127, 80, 65, 127, 127, 201, 195, 127, 69, 175, 80, 87, 39, 115, 240, 127, 175, 116, 168, 127, 127, 115, 96, 42, 23, 65, 65, 49, 15, 96, 141, 127, 127, 105, 127, 127, 127, };
  static const uint8_t initValueOcc7[64] = {  141, 141, 139, 146, 127, 144, 177, 218, 127, 63, 127, 115, 127, 164, 240, 194, 127, 127, 73, 97, 127, 190, 186, 128, 73, 16, 15, 88, 116, 127, 80, 161, 127, 116, 116, 240, 42, 166, 161, 230, 96, 47, 127, 127, 58, 88, 116, 109, 105, 116, 15, 61, 15, 80, 73, 155, 15, 15, 15, 45, 36, 73, 57, 121,};
  _MapOccupancy[0][0].init(initValueOcc0);
  _MapOccupancy[0][1].init(initValueOcc1);
  _MapOccupancy[0][2].init(initValueOcc2);
  _MapOccupancy[0][3].init(initValueOcc3);
  _MapOccupancy[0][4].init(initValueOcc4);
  _MapOccupancy[0][5].init(initValueOcc5);
  _MapOccupancy[0][6].init(initValueOcc6);
  _MapOccupancy[0][7].init(initValueOcc7);

  if (!forTrisoup)
    return;

  MapOBUFTriSoup[0].reset(15 - 1, 7);      // flag
  MapOBUFTriSoup[1].reset(15 - 1, 6);      // first bit position
  MapOBUFTriSoup[2].reset(15 - 1, 6 + 1);  // second bit position
  MapOBUFTriSoup[3].reset(15 - 2, 6 + 2);  // third bit position
  MapOBUFTriSoup[4].reset(13 - 1, 6 + 3);  // fourth bit position

  _OBUFleafNumberTrisoup = 0;
  std::fill_n(_BufferOBUFleavesTrisoup, kLeavesBufferSize, 0);

  // Intra
  static const uint8_t initValue0[128] = { 15,15,15,15,15,15,15,15,15,15,42,96,71,37,15,15,22,51,15,15,30,27,15,15,64,15,48,15,224,171,127,24,127,34,80,46,141,44,66,49,127,116,140,116,105,39,127,116,114,46,172,109,60,73,181,161,112,65,240,159,127,127,127,87,183,127,116,116,195,88,152,141,228,141,127,80,127,127,160,92,224,167,129,135,240,183,240,184,240,240,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127 };
  static const uint8_t initValue1[64] = { 116,127,118,15,104,56,97,15,96,15,29,15,95,15,46,15,196,116,182,53,210,104,163,69,169,15,114,15,121,15,167,63,240,127,184,92,240,163,197,77,239,73,179,59,213,48,185,108,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127, };
  static const uint8_t initValue2[128] = { 141,127,127,127,189,81,36,127,143,105,103,116,201,60,38,116,116,127,15,127,153,59,15,116,69,105,15,127,158,93,36,79,141,161,116,127,197,102,53,127,177,125,88,79,209,75,102,28,95,74,72,56,189,62,78,18,88,116,28,45,237,100,152,35,141,240,127,127,208,133,101,141,186,210,168,98,201,124,138,15,195,194,103,94,229,82,167,23,92,197,112,59,185,87,156,79,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127,127 };
  MapOBUFTriSoup[0].init(initValue0);
  MapOBUFTriSoup[1].init(initValue1);
  MapOBUFTriSoup[2].init(initValue2);
}

//============================================================================

void
GeometryOctreeContexts::clearMap()
{
  for (int j = 0; j < 2; j++)
    for (int i = 0; i < 8; i++) {
      _MapOccupancy[j][i].clear();
      _MapOccupancySparse[j][i].clear();
    }

  std::cout << "Size used buffer OBUF LEAF = " << _OBUFleafNumber << "\n";

  MapOBUFTriSoup[0].clear();
  MapOBUFTriSoup[1].clear();
  MapOBUFTriSoup[2].clear();
  MapOBUFTriSoup[3].clear();
  MapOBUFTriSoup[4].clear();

  // colocated edges
  refFrameEdgeKeys.clear();
  refFrameEdgeValue.clear();
  refFrameEdgeQP.clear();
}

//============================================================================

}  // namespace pcc
