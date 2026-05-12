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

#include "quantization.h"
#include "hls.h"
#include "ply.h"
#include <vector>
#include "pointset_processing.h"
#include "AttributeEncoder.h"
#include "AttributeDecoder.h"

namespace pcc {

//============================================================================

struct EncoderAttributeParams;

//============================================================================

void regionAdaptiveHierarchicalTransformRSO(
  const AttributeParameterSet& aps,
  const AttributeDescription& desc,
  const EncoderAttributeParams& attrEncParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int attribCount,
  const int voxelCount,
  int64_t* positions,
  attr_t* attributes,
  const attr_t* attributes_mc,
  PCCResidualsEncoder& encoder,
  point_t blockStart,
  point_t blockSizeMinusOne,
  RAHT::BlockBoundaries* blockBoundaries,
  RAHT::BlockRefBoundaries* blockRefBoundaries);

void regionAdaptiveHierarchicalInverseTransform(
  const AttributeParameterSet& aps,
  const AttributeDescription& desc,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int attribCount,
  const int voxelCount,
  int64_t* positions,
  attr_t* attributes,
  const attr_t* attributes_mc,
  PCCResidualsDecoder& decoder,
  point_t blockStart,
  point_t blockSizeMinusOne,
  RAHT::BlockBoundaries* blockBoundaries,
  RAHT::BlockRefBoundaries* blockRefBoundaries);

//============================================================================

namespace RAHT {

enum class NeighbourAvailable: uint8_t {
  kNeighbourBlock_No = 0,
  kNeighbourBlock_Yes = 1,
  kNeighbourBlock_Ignore = 2,
  kNeighbourBound_No = 4,
  kNeighbourBound_Yes = 5,
  kNeighbourBound_Ignore = 6
};

template <typename It1, typename It2>
struct RSO_OneLevelNeighboursTraversal {
  using It3 = std::vector<RahtBoundaryNode>::const_iterator;

  static constexpr int numParentNeigh = 27;
  static constexpr int numChildNeigh = 28;
  static constexpr int numBoundaryParentNeigh = 9;
  static constexpr int numBoundaryChildNeigh = 16;

  const It1 parentBegin;
  const It1 parentEnd;
  It1 parent[numParentNeigh];
  It1 curr;
  const It2 childBegin;
  const It2 childEnd;
  It2 child[numChildNeigh];
  const It3 boundaryParentBegin[3];
  const It3 boundaryParentEnd[3];
  It3 boundaryParent[3][numBoundaryParentNeigh];
  const It3 boundaryChildBegin[3];
  const It3 boundaryChildEnd[3];
  It3 boundaryChild[3][numBoundaryChildNeigh];
  const BlockRefBoundaries* refBounds;
  const int layer;
  const bool useChild;

  const int parentLevel = layer + 1;
  const int childLevel = layer;
  const int64_t maskXparentLevel = RSO_RAHT::maskX & ~RSO_RAHT::packX((1 << parentLevel) - 1);
  const int64_t maskYparentLevel = RSO_RAHT::maskY & ~RSO_RAHT::packY((1 << parentLevel) - 1);
  const int64_t maskZparentLevel = RSO_RAHT::maskZ & ~RSO_RAHT::packZ((1 << parentLevel) - 1);

  const int64_t maskXYparentLevel = maskXparentLevel | maskYparentLevel;
  const int64_t maskXYZparentLevel = maskXYparentLevel | maskZparentLevel;
  const int64_t offsetParentTube = RSO_RAHT::packZ(1 << parentLevel);
  const int64_t offsetBoundaryParentTube[3] = {
    RSO_RAHT::packY(1 << parentLevel),
    RSO_RAHT::packZ(1 << parentLevel),
    RSO_RAHT::packZ(1 << parentLevel)
  };
  const int64_t offsetChildTube = RSO_RAHT::packZ(1 << childLevel);
  const int64_t offsetBoundaryChildTube[3] = {
    RSO_RAHT::packY(1 << childLevel),
    RSO_RAHT::packZ(1 << childLevel),
    RSO_RAHT::packZ(1 << childLevel)
  };

  // TODO see how to avoid computing all the offsets at each object creations

  const std::array<int64_t, numParentNeigh> parentPosOffset = {
    RSO_RAHT::pack( -1, -1, -1 ) << parentLevel,
    RSO_RAHT::pack( -1, -1,  0 ) << parentLevel, // parentNeigh 7
    RSO_RAHT::pack( -1, -1, +1 ) << parentLevel,
    RSO_RAHT::pack( -1,  0, -1 ) << parentLevel, // parentNeigh 8
    RSO_RAHT::pack( -1,  0,  0 ) << parentLevel, // parentNeigh 9
    RSO_RAHT::pack( -1,  0, +1 ) << parentLevel, // parentNeigh 13
    RSO_RAHT::pack( -1, +1, -1 ) << parentLevel,
    RSO_RAHT::pack( -1, +1,  0 ) << parentLevel, // parentNeigh 15
    RSO_RAHT::pack( -1, +1, +1 ) << parentLevel,
    RSO_RAHT::pack(  0, -1, -1 ) << parentLevel, // parentNeigh 10
    RSO_RAHT::pack(  0, -1,  0 ) << parentLevel, // parentNeigh 11
    RSO_RAHT::pack(  0, -1, +1 ) << parentLevel, // parentNeigh 14
    RSO_RAHT::pack(  0,  0, -1 ) << parentLevel, // parentNeigh 12
    RSO_RAHT::pack(  0,  0,  0 ) << parentLevel, // parentNeigh 0
    RSO_RAHT::pack(  0,  0, +1 ) << parentLevel, // parentNeigh 3
    RSO_RAHT::pack(  0, +1, -1 ) << parentLevel, // parentNeigh 16
    RSO_RAHT::pack(  0, +1,  0 ) << parentLevel, // parentNeigh 2
    RSO_RAHT::pack(  0, +1, +1 ) << parentLevel, // parentNeigh 6
    RSO_RAHT::pack( +1, -1, -1 ) << parentLevel,
    RSO_RAHT::pack( +1, -1,  0 ) << parentLevel, // parentNeigh 17
    RSO_RAHT::pack( +1, -1, +1 ) << parentLevel,
    RSO_RAHT::pack( +1,  0, -1 ) << parentLevel, // parentNeigh 18
    RSO_RAHT::pack( +1,  0,  0 ) << parentLevel, // parentNeigh 1
    RSO_RAHT::pack( +1,  0, +1 ) << parentLevel, // parentNeigh 5
    RSO_RAHT::pack( +1, +1, -1 ) << parentLevel,
    RSO_RAHT::pack( +1, +1,  0 ) << parentLevel, // parentNeigh 4
    RSO_RAHT::pack( +1, +1, +1 ) << parentLevel,
  };

  static constexpr int currParentIdx = 13;
  static constexpr int numParentsPerTube = 3;
  static constexpr int numChildsPerTube = 4;

  int32_t parentAvailableMask;

  const std::array<int64_t, numChildNeigh> childPosOffset = {
    RSO_RAHT::pack( -1, -1, -1 ) << childLevel,
    RSO_RAHT::pack( -1, -1,  0 ) << childLevel, // parentNeigh 7 child 6
    RSO_RAHT::pack( -1, -1,  1 ) << childLevel, // parentNeigh 7 child 7
    RSO_RAHT::pack( -1, -1, +2 ) << childLevel,

    RSO_RAHT::pack( -1,  0, -1 ) << childLevel, // parentNeigh 8 child 5
    RSO_RAHT::pack( -1,  0,  0 ) << childLevel, // parentNeigh 9 child 4
    RSO_RAHT::pack( -1,  0,  1 ) << childLevel, // parentNeigh 9 child 5
    RSO_RAHT::pack( -1,  0, +2 ) << childLevel, // parentNeigh 13 child 4

    RSO_RAHT::pack( -1,  1, -1 ) << childLevel, // parentNeigh 8 child 7
    RSO_RAHT::pack( -1,  1,  0 ) << childLevel, // parentNeigh 9 child 6
    RSO_RAHT::pack( -1,  1,  1 ) << childLevel, // parentNeigh 9 child 7
    RSO_RAHT::pack( -1,  1, +2 ) << childLevel, // parentNeigh 13 child 6

    RSO_RAHT::pack( -1, +2, -1 ) << childLevel,
    RSO_RAHT::pack( -1, +2,  0 ) << childLevel, // parentNeigh 15 child 4
    RSO_RAHT::pack( -1, +2,  1 ) << childLevel, // parentNeigh 15 child 5
    RSO_RAHT::pack( -1, +2, +2 ) << childLevel,

    RSO_RAHT::pack(  0, -1, -1 ) << childLevel, // parentNeigh 10 child 3
    RSO_RAHT::pack(  0, -1,  0 ) << childLevel, // parentNeigh 11 child 2
    RSO_RAHT::pack(  0, -1,  1 ) << childLevel, // parentNeigh 11 child 3
    RSO_RAHT::pack(  0, -1, +2 ) << childLevel, // parentNeigh 14 child 2

    RSO_RAHT::pack(  0,  0, -1 ) << childLevel, // parentNeigh 12 child 1

    RSO_RAHT::pack(  0,  1, -1 ) << childLevel, // parentNeigh 12 child 3

    RSO_RAHT::pack(  1, -1, -1 ) << childLevel, // parentNeigh 10 child 7
    RSO_RAHT::pack(  1, -1,  0 ) << childLevel, // parentNeigh 11 child 6
    RSO_RAHT::pack(  1, -1,  1 ) << childLevel, // parentNeigh 11 child 7
    RSO_RAHT::pack(  1, -1, +2 ) << childLevel, // parentNeigh 14 child 6

    RSO_RAHT::pack(  1,  0, -1 ) << childLevel, // parentNeigh 12 child 5

    RSO_RAHT::pack(  1,  1, -1 ) << childLevel, // parentNeigh 12 child 7

    /* never available in lexicographic order

    RSO_RAHT::pack(  0, +2, -1 ) << childLevel, // parentNeigh 16 child 1
    RSO_RAHT::pack(  1, +2, -1 ) << childLevel, // parentNeigh 16 child 5

    RSO_RAHT::pack( +2, -1,  0 ) << childLevel, // parentNeigh 17 child 2
    RSO_RAHT::pack( +2, -1,  1 ) << childLevel, // parentNeigh 17 child 3

    RSO_RAHT::pack( +2,  0, -1 ) << childLevel, // parentNeigh 18 child 1

    RSO_RAHT::pack( +2,  1, -1 ) << childLevel  // parentNeigh 18 child 3
    */
  };

  int32_t childAvailableMask;

  const std::array<int64_t, numBoundaryParentNeigh> boundaryParentPosOffset[3] = {
    {
    RSO_RAHT::pack( -1, -1, /*-1*/0 ) << parentLevel,
    RSO_RAHT::pack( -1,  0, /*-1*/0 ) << parentLevel, // parentNeigh 8
    RSO_RAHT::pack( -1, +1, /*-1*/0 ) << parentLevel,
    RSO_RAHT::pack(  0, -1, /*-1*/0 ) << parentLevel, // parentNeigh 10
    RSO_RAHT::pack(  0,  0, /*-1*/0 ) << parentLevel, // parentNeigh 12
    RSO_RAHT::pack(  0, +1, /*-1*/0 ) << parentLevel, // parentNeigh 16
    RSO_RAHT::pack( +1, -1, /*-1*/0 ) << parentLevel,
    RSO_RAHT::pack( +1,  0, /*-1*/0 ) << parentLevel, // parentNeigh 18
    RSO_RAHT::pack( +1, -1, /*-1*/0 ) << parentLevel,
    },
    {
    RSO_RAHT::pack( -1, /*-1*/0, -1 ) << parentLevel,
    RSO_RAHT::pack( -1, /*-1*/0,  0 ) << parentLevel, // parentNeigh 7
    RSO_RAHT::pack( -1, /*-1*/0, +1 ) << parentLevel,
    RSO_RAHT::pack(  0, /*-1*/0, -1 ) << parentLevel, // parentNeigh 10
    RSO_RAHT::pack(  0, /*-1*/0,  0 ) << parentLevel, // parentNeigh 11
    RSO_RAHT::pack(  0, /*-1*/0, +1 ) << parentLevel, // parentNeigh 14
    RSO_RAHT::pack( +1, /*-1*/0, -1 ) << parentLevel,
    RSO_RAHT::pack( +1, /*-1*/0,  0 ) << parentLevel, // parentNeigh 17
    RSO_RAHT::pack( +1, /*-1*/0, +1 ) << parentLevel,
    },
    {
    RSO_RAHT::pack( /*-1*/0, -1, -1 ) << parentLevel,
    RSO_RAHT::pack( /*-1*/0, -1,  0 ) << parentLevel, // parentNeigh 7
    RSO_RAHT::pack( /*-1*/0, -1, +1 ) << parentLevel,
    RSO_RAHT::pack( /*-1*/0,  0, -1 ) << parentLevel, // parentNeigh 8
    RSO_RAHT::pack( /*-1*/0,  0,  0 ) << parentLevel, // parentNeigh 9
    RSO_RAHT::pack( /*-1*/0,  0, +1 ) << parentLevel, // parentNeigh 13
    RSO_RAHT::pack( /*-1*/0, +1, -1 ) << parentLevel,
    RSO_RAHT::pack( /*-1*/0, +1,  0 ) << parentLevel, // parentNeigh 15
    RSO_RAHT::pack( /*-1*/0, +1, +1 ) << parentLevel,
    }
  };

  int32_t boundaryParentAvailableMask[3];

  const std::array<int64_t, numBoundaryChildNeigh> boundaryChildPosOffset[3] = {
    {
    RSO_RAHT::pack( -1, -1, /*-1*/0 ) << childLevel,
    RSO_RAHT::pack( -1,  0, /*-1*/0 ) << childLevel, // parentNeigh 8 child 5
    RSO_RAHT::pack( -1,  1, /*-1*/0 ) << childLevel, // parentNeigh 8 child 7
    RSO_RAHT::pack( -1, +2, /*-1*/0 ) << childLevel,
    RSO_RAHT::pack(  0, -1, /*-1*/0 ) << childLevel, // parentNeigh 10 child 3
    RSO_RAHT::pack(  0,  0, /*-1*/0 ) << childLevel, // parentNeigh 12 child 1
    RSO_RAHT::pack(  0,  1, /*-1*/0 ) << childLevel, // parentNeigh 12 child 3
    RSO_RAHT::pack(  0, +2, /*-1*/0 ) << childLevel, // parentNeigh 16 child 1
    RSO_RAHT::pack(  1, -1, /*-1*/0 ) << childLevel, // parentNeigh 10 child 7
    RSO_RAHT::pack(  1,  0, /*-1*/0 ) << childLevel, // parentNeigh 12 child 5
    RSO_RAHT::pack(  1,  1, /*-1*/0 ) << childLevel, // parentNeigh 12 child 7
    RSO_RAHT::pack(  1, +2, /*-1*/0 ) << childLevel, // parentNeigh 16 child 5
    RSO_RAHT::pack( +2, -1, /*-1*/0 ) << childLevel,
    RSO_RAHT::pack( +2,  0, /*-1*/0 ) << childLevel, // parentNeigh 18 child 1
    RSO_RAHT::pack( +2,  1, /*-1*/0 ) << childLevel, // parentNeigh 18 child 3
    RSO_RAHT::pack( +2, +2, /*-1*/0 ) << childLevel,
    },
    {
    RSO_RAHT::pack( -1, /*-1*/0, -1 ) << childLevel,
    RSO_RAHT::pack( -1, /*-1*/0,  0 ) << childLevel, // parentNeigh 7 child 6
    RSO_RAHT::pack( -1, /*-1*/0,  1 ) << childLevel, // parentNeigh 7 child 7
    RSO_RAHT::pack( -1, /*-1*/0, +2 ) << childLevel,
    RSO_RAHT::pack(  0, /*-1*/0, -1 ) << childLevel, // parentNeigh 10 child 3
    RSO_RAHT::pack(  0, /*-1*/0,  0 ) << childLevel, // parentNeigh 11 child 2
    RSO_RAHT::pack(  0, /*-1*/0,  1 ) << childLevel, // parentNeigh 11 child 3
    RSO_RAHT::pack(  0, /*-1*/0, +2 ) << childLevel, // parentNeigh 14 child 2
    RSO_RAHT::pack(  1, /*-1*/0, -1 ) << childLevel, // parentNeigh 10 child 7
    RSO_RAHT::pack(  1, /*-1*/0,  0 ) << childLevel, // parentNeigh 11 child 6
    RSO_RAHT::pack(  1, /*-1*/0,  1 ) << childLevel, // parentNeigh 11 child 7
    RSO_RAHT::pack(  1, /*-1*/0, +2 ) << childLevel, // parentNeigh 14 child 6
    RSO_RAHT::pack( +2, /*-1*/0, -1 ) << childLevel,
    RSO_RAHT::pack( +2, /*-1*/0,  0 ) << childLevel, // parentNeigh 17 child 2
    RSO_RAHT::pack( +2, /*-1*/0,  1 ) << childLevel, // parentNeigh 17 child 3
    RSO_RAHT::pack( +2, /*-1*/0, +2 ) << childLevel,
    },
    {
    RSO_RAHT::pack( /*-1*/0, -1, -1 ) << childLevel,
    RSO_RAHT::pack( /*-1*/0, -1,  0 ) << childLevel, // parentNeigh 7 child 6
    RSO_RAHT::pack( /*-1*/0, -1,  1 ) << childLevel, // parentNeigh 7 child 7
    RSO_RAHT::pack( /*-1*/0, -1, +2 ) << childLevel,
    RSO_RAHT::pack( /*-1*/0,  0, -1 ) << childLevel, // parentNeigh 8 child 5
    RSO_RAHT::pack( /*-1*/0,  0,  0 ) << childLevel, // parentNeigh 9 child 4
    RSO_RAHT::pack( /*-1*/0,  0,  1 ) << childLevel, // parentNeigh 9 child 5
    RSO_RAHT::pack( /*-1*/0,  0, +2 ) << childLevel, // parentNeigh 13 child 4
    RSO_RAHT::pack( /*-1*/0,  1, -1 ) << childLevel, // parentNeigh 8 child 7
    RSO_RAHT::pack( /*-1*/0,  1,  0 ) << childLevel, // parentNeigh 9 child 6
    RSO_RAHT::pack( /*-1*/0,  1,  1 ) << childLevel, // parentNeigh 9 child 7
    RSO_RAHT::pack( /*-1*/0,  1, +2 ) << childLevel, // parentNeigh 13 child 6
    RSO_RAHT::pack( /*-1*/0, +2, -1 ) << childLevel,
    RSO_RAHT::pack( /*-1*/0, +2,  0 ) << childLevel, // parentNeigh 15 child 4
    RSO_RAHT::pack( /*-1*/0, +2,  1 ) << childLevel, // parentNeigh 15 child 5
    RSO_RAHT::pack( /*-1*/0, +2, +2 ) << childLevel,
    }
  };

  int32_t boundaryChildAvailableMask[3];

  RSO_OneLevelNeighboursTraversal(
    const BlockRefBoundaries* refBounds,
    It1 parentBegin,
    It1 parentEnd,
    bool useChild,
    It2 childBegin,
    It2 childEnd,
    int layer)
    : refBounds(refBounds)
    , parentBegin(parentBegin)
    , parentEnd(parentEnd)
    , useChild(useChild)
    , childBegin(childBegin)
    , childEnd(childEnd)
    , boundaryParentBegin {
        refBounds ? refBounds->layer_cbegin(0, layer + 1) : It3(),
        refBounds ? refBounds->layer_cbegin(1, layer + 1) : It3(),
        refBounds ? refBounds->layer_cbegin(2, layer + 1) : It3()
      }
    , boundaryParentEnd {
        refBounds ? refBounds->layer_cend(0, layer + 1) : It3(),
        refBounds ? refBounds->layer_cend(1, layer + 1) : It3(),
        refBounds ? refBounds->layer_cend(2, layer + 1) : It3()
      }
    , boundaryChildBegin {
        refBounds ? refBounds->layer_cbegin(0, layer) : It3(),
        refBounds ? refBounds->layer_cbegin(1, layer) : It3(),
        refBounds ? refBounds->layer_cbegin(2, layer) : It3()
    }
    , boundaryChildEnd{
        refBounds ? refBounds->layer_cend(0, layer) : It3(),
        refBounds ? refBounds->layer_cend(1, layer) : It3(),
        refBounds ? refBounds->layer_cend(2, layer) : It3()
    }
    , layer(layer)
  {
    curr = parentBegin;
    std::fill_n(parent, numParentNeigh, parentBegin);
    std::fill_n(&child[0], numChildNeigh, childBegin);
    // ** Boundaries **
    std::fill_n(boundaryParentAvailableMask, 3, 0);
    std::fill_n(boundaryChildAvailableMask, 3, 0);
    if (refBounds) {
      // parents
      for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
        if (refBounds->nodes[faceIdx]) {
          std::fill_n(boundaryParent[faceIdx], numBoundaryParentNeigh,
            boundaryParentBegin[faceIdx]);
        }
      }
      // Childs
      for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
        if (refBounds->nodes[faceIdx]) {
          std::fill_n(boundaryChild[faceIdx], numBoundaryChildNeigh,
            boundaryChildBegin[faceIdx]);
        }
      }
    }
    updateIterators();
  }

  void next() {
    assert(curr < parentEnd);
    ++curr;
    updateIterators();
  }

  template <int firstNodeToUpdate, int numNodesToUpdate, int numNodesPerTube,
    class It, typename T0, typename T1, typename T2, typename T3>
  void updateIteratorsTubes(T0& nodes, const It& itEnd, T1& mask, const T2& curPos, const T3& offsets, const T2& offsetInTube)
  {
    for (int i = firstNodeToUpdate; i < firstNodeToUpdate + numNodesToUpdate; ++i) {
      auto neighpos = curPos + offsets[i];
      if (numNodesPerTube > 1)
        nodes[i] = nodes[i + 1];
      while (nodes[i] != itEnd && nodes[i]->pos < neighpos) ++nodes[i];
      mask |= (nodes[i] != itEnd && nodes[i]->pos == neighpos) << i;
      const int tubeStart = i;
      for (int j = tubeStart + 1; j < tubeStart + numNodesPerTube; ++j, ++i) {
        neighpos += offsetInTube;
        nodes[j] = nodes[i];
        nodes[j] += nodes[j] != itEnd && nodes[j]->pos < neighpos;
        mask |= (nodes[j] != itEnd && nodes[j]->pos == neighpos) << j;
      }
    }
  }

  void updateIterators()
  {
    if (curr != parentEnd) {
      const auto curPos = curr->pos;

      parentAvailableMask = 0;
      updateIteratorsTubes<0, numParentNeigh, numParentsPerTube>(
        parent, parentEnd, parentAvailableMask, curPos, parentPosOffset, offsetParentTube);

      if (useChild) {
        childAvailableMask = 0;
        updateIteratorsTubes<0, 20, numChildsPerTube>(
          child, childEnd, childAvailableMask, curPos, childPosOffset, offsetChildTube);
        updateIteratorsTubes<20, 1, 1>(
          child, childEnd, childAvailableMask, curPos, childPosOffset, offsetChildTube);
        updateIteratorsTubes<21, 1, 1>(
          child, childEnd, childAvailableMask, curPos, childPosOffset, offsetChildTube);
        updateIteratorsTubes<22, 4, numChildsPerTube>(
          child, childEnd, childAvailableMask, curPos, childPosOffset, offsetChildTube);
        updateIteratorsTubes<26, 1, 1>(
          child, childEnd, childAvailableMask, curPos, childPosOffset, offsetChildTube);
        updateIteratorsTubes<27, 1, 1>(
          child, childEnd, childAvailableMask, curPos, childPosOffset, offsetChildTube);
      }

      if (refBounds) {
        static constexpr int64_t maskFaceUpdate[3] =
          { RSO_RAHT::maskZ, RSO_RAHT::maskY, RSO_RAHT::maskX };

        for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
          boundaryParentAvailableMask[faceIdx] = 0;
          // only update if coordinate is 0 (slab blocks origin are offset to 0)
          if (!(curPos & maskFaceUpdate[faceIdx])) {
            updateIteratorsTubes<0, numBoundaryParentNeigh, numParentsPerTube>(
              boundaryParent[faceIdx], boundaryParentEnd[faceIdx],
              boundaryParentAvailableMask[faceIdx], curPos,
              boundaryParentPosOffset[faceIdx], offsetBoundaryParentTube[faceIdx]);
          }
        }
        if (useChild) {
          for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
            boundaryChildAvailableMask[faceIdx] = 0;
            // only update if coordinate is 0 (slab blocks origin are offset to 0)
            if (!(curPos & maskFaceUpdate[faceIdx])) {
              updateIteratorsTubes<0, numBoundaryChildNeigh, numChildsPerTube>(
                boundaryChild[faceIdx], boundaryChildEnd[faceIdx],
                boundaryChildAvailableMask[faceIdx], curPos,
                boundaryChildPosOffset[faceIdx], offsetBoundaryChildTube[faceIdx]);
            }
          }
        }
      }
    }
  }

  void getNeighborsMode(
    uint8_t occupancy,
    int& voteInterWeight,
    int& voteIntraWeight,
    int* voteMode,
    int& numCoeffNot0,
    int& numCoeffTotal,
    int boundariesPosMask)
  {
    static const uint8_t neighMasks[numParentNeigh] = {
        1,   3,   2,   5,  15,  10,   4,  12,   8,
       17,  51,  34,  85, 255, 170,  68, 204, 136,
       16,  48,  32,  80, 240, 160,  64, 192, 128 };

    constexpr int numPN = 18;
    static constexpr int pIdx[numPN] =
    { 1, 3, 4, 5, 7, 9, 10, 11, 12, 14, 15, 16,  17, 19, 21, 22, 23, 25 };

    constexpr int numBPN = 5;
    static constexpr int bPIdx[numBPN] =
    { 1, 3, 4, 5, 7 };

    static constexpr int bPIdxPIdx[3][numBPN] =
    {{ 3, 9, 12, 15, 21 }, { 1, 9, 10, 11, 19 }, { 1, 3, 4, 5, 7 }};

    int voteAverageMode[4] = { 0,0,0,0 }; // Null, intra, inter, size;

    for (int j = 0; j < numPN; ++j) {
      const int i = pIdx[j];
      const int32_t bitMask = 1 << i;

      if (occupancy & neighMasks[i] && parentAvailableMask & bitMask) {
        auto& uncle = parent[i];
        voteMode[uncle->mode]++;

        if (uncle->decoded) {
          auto cousin_mode = uncle->child_mode;
          auto cousin_Wintra = uncle->childsWintra;
          auto num_cousin_from_uncle = popcnt(uncle->occupancy);
          voteMode[cousin_mode] += 3 * num_cousin_from_uncle;
          voteAverageMode[1] += 3 * num_cousin_from_uncle * cousin_Wintra;
          voteAverageMode[2] += 3 * num_cousin_from_uncle * (128 - cousin_Wintra);

          numCoeffNot0 += uncle->CntNonZero;
          numCoeffTotal += num_cousin_from_uncle;
        }
        else {
          voteAverageMode[1] += uncle->Wintra;
          voteAverageMode[2] += 128 - uncle->Wintra;
        }
      }
    }

    if (refBounds && boundariesPosMask & 7) {
      for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
        if(boundaryParentAvailableMask[faceIdx]) {
          for (int j = 0; j < numBPN; ++j) {
            const int i = bPIdx[j];
            const int idx = bPIdxPIdx[faceIdx][j];
            const int32_t bitMask = 1 << i;

            if (occupancy & neighMasks[idx] && boundaryParentAvailableMask[faceIdx] & bitMask) {
              const auto& uncle = *(boundaryParent[faceIdx][i]);
              voteMode[uncle.mode]++;

              auto cousin_mode = uncle.child_mode;
              auto cousin_Wintra = uncle.childsWintra;
              auto num_cousin_from_uncle = popcnt(uncle.occupancy);
              voteMode[cousin_mode] += 3 * num_cousin_from_uncle;
              voteAverageMode[1] += 3 * num_cousin_from_uncle * cousin_Wintra;
              voteAverageMode[2] += 3 * num_cousin_from_uncle * (128 - cousin_Wintra);
              numCoeffNot0 += uncle.CntNonZero;
              numCoeffTotal += num_cousin_from_uncle;
            }
          }
        }
      }
    }

    auto& _parent = parent[currParentIdx];

    voteIntraWeight = voteAverageMode[1] + 4 * _parent->Wintra;
    voteInterWeight = voteAverageMode[2] + 4 * (128 - _parent->Wintra);
  }

  template<bool haarFlag, int numAttrs>
  inline void
  intraDcPred(
    const int occupancy,
    It1 first,
    It2 firstChild,
    int64_t* predBuf,
    const RahtPredictionParams& rahtPredParams,
    int boundariesPosMask,
    BlockRefBoundaries* refBounds,
    bool inheritDc) const;
};

//============================================================================
// expand a set of eight weights into three levels

template <bool haarFlag>
struct mkWeightTree {
  template<bool skipkernel, class Kernel = RahtKernel>
  static void
  apply(int64_t weights[8 + 8 + 8 + 8 + 24])
  {
    auto in = &weights[0];
    auto out = &weights[8];

    for (int i = 0; i < 4; i++) {
      out[0] = in[0] + in[1];
      out[4] = (in[0] && in[1]) * out[0];  // single node, no high frequencies
      in += 2;
      out++;
    }
    out += 4;
    for (int i = 0; i < 4; i++) {
      out[0] = in[0] + in[1];
      out[4] = (in[0] && in[1]) * out[0];  // single node, no high frequencies
      in += 2;
      out++;
    }
    out += 4;
    for (int i = 0; i < 4; i++) {
      out[0] = in[0] + in[1];
      out[4] = (in[0] && in[1]) * out[0];  // single node, no high frequencies
      in += 2;
      out++;
    }

    for (int i = 0; i < 24; i += 2) {
      weights[i + 32] = (!!weights[i]) << kFPFracBits;
      weights[i + 33] = (!!weights[i + 1]) << kFPFracBits;

      if (weights[i] && weights[i + 1]) {
        if (skipkernel) {
          weights[i + 32] = weights[i];
          weights[i + 33] = weights[i + 1];
        }
        else {
          Kernel w(weights[i], weights[i + 1]);
          weights[i + 32] = w.getW0();
          weights[i + 33] = w.getW1();
        }
      }
    }
  }
};

//----------------------------------------------------------------------------

template <>
struct mkWeightTree<true> {
  template<bool skipkernel, class Kernel = HaarKernel>
  static void
  apply(bool weights[8 + 8 + 8 + 8 + 24])
  {
    auto in = &weights[0];
    auto out = &weights[8];

    for (int i = 0; i < 4; i++) {
      out[0] = in[0] || in[1];
      out[4] = in[0] && in[1];  // single node, no high frequencies
      in += 2;
      out++;
    }
    out += 4;
    for (int i = 0; i < 4; i++) {
      out[0] = in[0] || in[1];
      out[4] = in[0] && in[1];  // single node, no high frequencies
      in += 2;
      out++;
    }
    out += 4;
    for (int i = 0; i < 4; i++) {
      out[0] = in[0] || in[1];
      out[4] = in[0] && in[1];  // single node, no high frequencies
      in += 2;
      out++;
    }

    for (int i = 0; i < 24; i++) {
      weights[i + 32] = weights[i];
    }
  }
};

//============================================================================
// Generate the spatial prediction of a block.

template<typename It1, typename It2>
template<bool haarFlag, int numAttrs>
void
RSO_OneLevelNeighboursTraversal<It1, It2>::
intraDcPred(
  const int occupancy,
  It1 first,
  It2 firstChild,
  int64_t* predBuf,
  const RahtPredictionParams& rahtPredParams,
  int boundariesPosMask,
  BlockRefBoundaries* refBounds,
  bool inheritDc) const
{
  constexpr int currParentIdxPred = 9;
  constexpr int numParentNeighPred = 19;
  constexpr int numChildNeighPred = 30; // 24 + up to 6 more on boundaries
  // neighbors values
  int32_t nv[numParentNeighPred+numChildNeighPred][numAttrs];
  // mask for neighbors availability
  int64_t neigh_mask = 0;

  // parent neighbors prediction indices = neighbors indices
  static const int npIdx[numParentNeigh] = {
    -1,  0, -1,  1,  2,  3, -1,  4, -1,
     5,  6,  7,  8,  9, 10, 11, 12, 13,
    -1, 14, -1, 15, 16, 17, -1, 18, -1,
  };
  static constexpr int npMask = 0x2EBFEBA;

  // boundary parent neighbors prediction indices
  static const int bnpIdx[3][numBoundaryParentNeigh] = {
    {
      -1, /* 0, -1, /**/  1, /* 2,  3, /**/ -1, /* 4, -1, */
       5, /* 6,  7, /**/  8, /* 9, 10, /**/ 11, /*12, 13, */
      -1, /*14, -1, /**/ 15, /*16, 17, /**/ -1, /*18, -1, */
    },
    {
      -1,  0, -1, /* /  1,  2,  3, /* / -1,  4, -1, */
       5,  6,  7, /* /  8,  9, 10, /* / 11, 12, 13, */
      -1, 14, -1, /* / 15, 16, 17, /* / -1, 18, -1, */
    },
    {
      -1,  0, -1, /**/ 1,  2,  3,  /**/ -1,  4, -1,
    /* 5,  6,  7, /* /  8,  9, 10, /* / 11, 12, 13, */
    /*-1, 14, -1, /* / 15, 16, 17, /* / -1, 18, -1, */
    },
  };
  static constexpr int bnpMask = 0xBA;

  // child neighbors prediction indices
  static const int cnpIdx[numChildNeigh] = {
    -1, 19, 20, -1, /**/ 21, 22, 23, 24, /**/ 25, 26, 27, 28, /**/ -1, 29, 30, -1,
    31, 32, 33, 34, /**/ 35,             /**/ 36,             /* / 37,*/
    38, 39, 40, 41, /**/ 42,             /**/ 43,             /* / 44,*/
  /*-1, 45, 46, -1, /* / 47,             /* / 48,             /* / -1 */
  };
  static constexpr int cnpMask = 0xFFF6FF6;

  // boundary child neighbors prediction indices
  static const int bcnpIdx[3][numBoundaryChildNeigh] = {
    {
      -1, /*19, 20, -1, */ 21, /*22, 23, 24, */ 25, /*26, 27, 28, */ -1, /*29, 30, -1,*/
      31, /*32, 33, 34, */ 35,             /**/ 36,             /**/ 37,
      38, /*39, 40, 41, */ 42,             /**/ 43,             /**/ 44,
      -1, /*45, 46, -1, */ 47,             /**/ 48,             /**/ -1,
    },
    {
      -1, 19, 20, -1, /* / 21, 22, 23, 24, /* / 25, 26, 27, 28, /* / -1, 29, 30, -1,*/
      31, 32, 33, 34, /* / 35,             /* / 36,             /* / 37,*/
      38, 39, 40, 41, /* / 42,             /* / 43,             /* / 44,*/
      -1, 45, 46, -1, /* / 47,             /* / 48,             /* / -1,*/
    },
    {
      -1, 19, 20, -1, /**/ 21, 22, 23, 24, /**/ 25, 26, 27, 28, /**/ -1, 29, 30, -1,/*
      31, 32, 33, 34, /* / 35,             /* / 36,             /* / 37,
      38, 39, 40, 41, /* / 42,             /* / 43,             /* / 44,
      -1, 45, 46, -1, /* / 47,             /* / 48,             /* / -1,*/
    },
  };
  static constexpr int bcnpMask = 0x6FF6;

  if (inheritDc) {
    int32_t tobeMasked = parentAvailableMask & npMask;
    for (int i = 0; i < numParentNeigh; i++, tobeMasked >>= 1) {
      if (tobeMasked & 1) {
        const int neighIdx = npIdx[i];
        const auto& node = *parent[i];
        neigh_mask |= 1LL << neighIdx;
        for (int k = 0; k < numAttrs; ++k)
          nv[neighIdx][k] = node.reconstructedAttr[k];
      }
    }
  }

  if (refBounds && boundariesPosMask & 7) {
    for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
      if(boundaryParentAvailableMask[faceIdx]) {
        int32_t tobeMasked = boundaryParentAvailableMask[faceIdx] & bnpMask;
        for (int i = 0; i < numBoundaryParentNeigh; i++, tobeMasked >>= 1) {
          if (tobeMasked & 1) {
            const int neighIdx = bnpIdx[faceIdx][i];
            const auto& node = *boundaryParent[faceIdx][i];
            neigh_mask |= 1LL << neighIdx;
            for (int k = 0; k < numAttrs; ++k)
              nv[neighIdx][k] = node.reconstructedAttr[k];
          }
        }
      }
    }
  }

  if (rahtPredParams.subnode_prediction_enabled_flag) {
    int32_t tobeMasked = childAvailableMask & cnpMask;
    for (int i = 0; i < numChildNeigh; i++, tobeMasked >>= 1) {
      if (tobeMasked & 1) {
        const int neighIdx = cnpIdx[i];
        const auto& node = *child[i];
        neigh_mask |= 1LL << neighIdx;
        for (int k = 0; k < numAttrs; ++k)
          nv[neighIdx][k] = node.reconstructedAttr[k];
      }
    }

    if (refBounds && boundariesPosMask & 7) {
      for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
        if(boundaryChildAvailableMask[faceIdx]) {
          int32_t tobeMasked = boundaryChildAvailableMask[faceIdx] & bcnpMask;
          for (int i = 0; i < numBoundaryChildNeigh; i++, tobeMasked >>= 1) {
            if (tobeMasked & 1) {
              const int neighIdx = bcnpIdx[faceIdx][i];
              const auto& node = *boundaryChild[faceIdx][i];
              neigh_mask |= 1LL << neighIdx;
              for (int k = 0; k < numAttrs; ++k)
                nv[neighIdx][k] = node.reconstructedAttr[k];
            }
          }
        }
      }
    }
  }

  int64_t limitLow = 0;
  int64_t limitHigh = INT64_MAX;

  if (inheritDc) {
    constexpr int64_t ratioThreshold1 = 2;
    constexpr int64_t ratioThreshold2 = 25;
    limitLow = ratioThreshold1 * nv[currParentIdxPred][0];
    limitHigh = ratioThreshold2 * nv[currParentIdxPred][0];
  }

  static const int uncleSelect[8][7] = {
    // child 0, parents for pred:
    { 0, 1, 2, 5, 6, 8, 9 },
    // child 1, parents for pred:
    { 0, 2, 3, 6, 7, 9, 10 },
    // child 2, parents for pred:
    { 1, 2, 4, 8, 9, 11, 12 },
    // child 3, parents for pred:
    { 2, 3, 4, 9, 10, 12, 13 },
    // child 4, parents for pred:
    { 5, 6, 8, 9, 14, 15, 16 },
    // child 5, parents for pred:
    { 6, 7, 9, 10, 14, 16, 17 },
    // child 6, parents for pred:
    { 8, 9, 11, 12, 15, 16, 18 },
    // child 7, parents for pred:
    { 9, 10, 12, 13, 16, 17, 18 },
  };

  const auto& pW = rahtPredParams.predWeightParent;

  // weight for replacement childs
  const auto& cW = rahtPredParams.predWeightChild;

  // childIdx, predIdx -> childNeighborIdx
  static const int cousinSelect[8][7] =
  {
    // child 0
    {    19, 21, 22, 31, 32, 35, -1 },
    // child 1
    { 20,    23, 24, 33, 34, -1, -1 },
    // child 2
    { 25, 26,    29, 36, -1, 37, -1 },
    // child 3
    { 27, 28, 30,    -1, -1, -1, -1 },
    // child 4
    { 38, 39, 42, -1,    45, 47, -1 },
    // child 5
    { 40, 41, -1, -1, 46,    -1, -1 },
    // child 6
    { 43, -1, 44, -1, 48, -1,    -1 },
    // child 7
    { -1, -1, -1, -1, -1, -1, -1,   },
  };

  int weightSum[8];

  std::fill_n(predBuf, numAttrs * 8, 0);
  std::fill_n(weightSum, 8, 0);

  int tobeMasked = occupancy;
  int weight;
  for (uint8_t childIdx = 0; childIdx < 8; childIdx++, tobeMasked >>= 1) {
    if (tobeMasked & 1) {
      for (int pIdx = 0; pIdx < 7; ++pIdx) {
        int uncleIdx = uncleSelect[childIdx][pIdx];
        int64_t rangeCheck = 10 * int64_t(nv[uncleIdx][0]);
        if (neigh_mask & (1 << uncleIdx) && rangeCheck > limitLow && rangeCheck < limitHigh) {

          int cousinIdx = cousinSelect[childIdx][pIdx];
          int idx = uncleIdx;
          if (rahtPredParams.subnode_prediction_enabled_flag && cousinIdx >= 0 && (neigh_mask & (1LL << cousinIdx))) {
            weight = cW[cousinIdx -numParentNeighPred];
            idx = cousinIdx;
          }
          else {
            weight = pW[uncleIdx];
          }
          weightSum[childIdx] += weight;
          for (int k = 0; k < numAttrs; ++k)
            predBuf[8 * k + childIdx] += int64_t(nv[idx][k]) * weight;

        }
      }
    }
  }

  // finalize prediction for root node
  if (!inheritDc) {
    // compute DC prediction
    int wSum = 0;
    int64_t pred[numAttrs] = {};
    for (int i = 0, occ = occupancy; i < 8; i++, occ >>= 1) {
      if (occ & 1 && weightSum[i]) {
        wSum += weightSum[i];
        for (int k = 0; k < numAttrs; k++)
          pred[k] += predBuf[8 * k + i];
      }
    }
    // set DC prediction
    if (wSum > 0) {
      for (int i = 0; i < 8; i++) {
        weightSum[i] = wSum;
        for (int k = 0; k < numAttrs; k++)
          predBuf[8 * k + i] = pred[k];
      }
    }
  }

  // normalise
  for (int i = 0, occ = occupancy; i < 8; i++, occ >>= 1) {
    if (occ & 1) {
      int w = weightSum[i];
      if (w > 1) {
        int shift = 0;
        uint64_t divisor = divInvDivisorApprox(w, shift);
        for (int k = 0; k < numAttrs; k++)
          predBuf[8 * k + i] = (predBuf[8 * k + i] * divisor) >> shift;
      }
      if (haarFlag) {
        for (int k = 0; k < numAttrs; k++) {
          predBuf[8 * k + i] = fpExpand<kFPFracBits>(fpReduce<kFPFracBits>(predBuf[8 * k + i]));
        }
      }
    }
  }
}

//============================================================================
// Cross Chroma Component Prediction

static constexpr int kCCCPFiltPrecisionbits = 4;
static constexpr int kCCCPtemplatesize = 128;

struct PCCRAHTComputeCCCP {
  int computeCrossChromaComponentPredictionCoeff(Vec3<int64_t>& CCCPcorr)
  {
    if (window.size() == kCCCPtemplatesize) {
      sum -= window.front();
      window.pop_front();
    }

    CCCPcorr[0] = fpReduce<10>(CCCPcorr[0]);
    CCCPcorr[1] >>= 10;
    sum += CCCPcorr;
    window.push(CCCPcorr);

    constexpr int64_t threshold = 1 << kCCCPFiltPrecisionbits; // one

    int64_t scale = 0;
    if (sum[0] && sum[1]) {
      scale = divApprox(sum[0], sum[1], kCCCPFiltPrecisionbits);
      scale = PCCClip(scale, -threshold, threshold);
    }
    return scale;
  }

private:
  Vec3<int64_t> sum = 0;
  ringbuf<Vec3<int64_t>> window = ringbuf<Vec3<int64_t>>(kCCCPtemplatesize);
};

//============================================================================
// Compute CCRP filter by floor division

static constexpr int kCCRPFiltPrecisionbits = 4;
static constexpr int kCCRPtemplatesize = 64;

struct CCRPFilter {
  void update(const Vec3<int64_t>& CCRPcorr) {
    if (window.size() == kCCRPtemplatesize) {
      sum -= window.front();
      window.pop();
    }

    window.push(CCRPcorr);
    sum += CCRPcorr;

    constexpr int64_t threshold = 1 << kCCRPFiltPrecisionbits - 1; // one half

    constexpr int prec = kCCRPFiltPrecisionbits + 1;
    Vec3<int64_t> lsum = sum;
    // dirty hack to not overflow,
    // to be removed when overflow is better avoided
    int maxPrec = ilog2(uint64_t(std::max<int64_t>({1LL, std::abs(sum[1]), std::abs(sum[2])}))) + 1;
    if (maxPrec > 45 - prec) { // max precision for divApprox
      int shift = maxPrec - 45 + prec;
      lsum[0] >>= shift;
      lsum[1] = fpReduce(shift, lsum[1]);
      lsum[2] = fpReduce(shift, lsum[2]);
    }

    int64_t temp = divApprox(fpExpand<prec>(lsum[1]), lsum[0] + 1, kCCRPFiltPrecisionbits);
    _YCbFilt = temp >= 0 ? temp >> prec : -(-temp >> prec);
    _YCbFilt = PCCClip(_YCbFilt, -threshold, threshold);

    temp = divApprox(fpExpand<prec>(lsum[2]), lsum[0] + 1, kCCRPFiltPrecisionbits);
    _YCrFilt = temp >= 0 ? temp >> prec : -(-temp >> prec);
    _YCrFilt = PCCClip(_YCrFilt, -threshold, threshold);
  }

  int64_t getYCbFilt() const { return _YCbFilt; }
  int64_t getYCrFilt() const { return _YCrFilt; }

private:
  Vec3<int64_t> sum {0, 0, 0};

  int64_t _YCbFilt = 0;
  int64_t _YCrFilt = 0;

  ringbuf<Vec3<int64_t>> window = ringbuf<Vec3<int64_t>>(kCCRPtemplatesize);
};

//============================================================================

} /* namespace RAHT */

//============================================================================

} /* namespace pcc */
