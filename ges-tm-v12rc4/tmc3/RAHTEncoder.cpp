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

#include "RAHT.h"

#include <cassert>
#include <cfloat>
#include <cinttypes>
#include <climits>
#include <cstddef>
#include <utility>
#include <vector>
#include <stdio.h>

#include "PCCTMC3Common.h"
#include "PCCTMC3Encoder.h"
#include "PCCMisc.h"

using pcc::attr::PredMode;
using namespace pcc::RAHT;

namespace pcc {

//============================================================================
// remove any non-unique leaves from a level in the uraht tree

template<bool haarFlag, int numAttrs, typename UrahtNodeEncoder>
int
reduceUniqueEncoder(
  int numNodes,
  std::vector<UrahtNodeEncoder>* weightsIn,
  std::vector<UrahtNodeEncoder>* weightsOut)
{
  // process a single level of the tree
  int64_t posPrev = -1;
  auto weightsInWrIt = weightsIn->begin();
  auto weightsInRdIt = weightsIn->cbegin();
  for (int i = 0; i < numNodes; i++) {
    const auto& node = *weightsInRdIt++;

    // copy across unique nodes
    if (node.pos != posPrev) {
      posPrev = node.pos;
      *weightsInWrIt++ = node;
      continue;
    }

    // duplicate node
    (weightsInWrIt - 1)->weight += node.weight;
    weightsOut->push_back(node);
    if (haarFlag) {
      for (int k = 0; k < numAttrs; k++) {
        auto temp = node.sumAttr[k] - (weightsInWrIt - 1)->sumAttr[k];
        (weightsInWrIt - 1)->sumAttr[k] += temp >> 1;
        weightsOut->back().sumAttr[k] = temp;
      }
    }
    else {
      for (int k = 0; k < numAttrs; k++)
        (weightsInWrIt - 1)->sumAttr[k] += node.sumAttr[k];
    }
    for (int k = 0; k < numAttrs; k++)
      (weightsInWrIt - 1)->sumAttrInter[k] += node.sumAttrInter[k];
  }

  auto res = std::distance(weightsIn->begin(), weightsInWrIt);
  weightsIn->resize(res);
  // number of nodes in next level
  return res;
}

//============================================================================
template<bool haarFlag, int numAttrs>
int
reduceDepthRasterEncoder(
  int level,
  const std::vector<UrahtNodeRSOEncoder>& weightsIn,
  std::vector<UrahtNodeRSOEncoder>& weightsOut)
{
  // level of child nodes
  const int levelIn = level - 1;

  OneLevelRAHTNodesTraversal<decltype(weightsIn.cbegin())>
    traversal(levelIn, weightsIn.cbegin(), weightsIn.cend());

  while (!traversal.finished()) {
    traversal.determineOccupancyAndPosition();

    decltype(weightsOut.begin())::value_type outNode;

    outNode.pos = traversal.nodePos;
    outNode.occupancy = traversal.nodeOccupancy;

    outNode.weight = 0;
    for (int k = 0; k < numAttrs; k++) {
      outNode.sumAttr[k] = 0;
      outNode.sumAttrInter[k] = 0;
    }

    for (int nodeIdx = 0; nodeIdx < 8; ++nodeIdx)
      if (outNode.occupancy & (1 << nodeIdx)) {
        const auto& node = *traversal.nodeIt[nodeIdx];
        outNode.weight += node.weight;
        if (!haarFlag) {
          for (int k = 0; k < numAttrs; k++) {
            outNode.sumAttr[k] += node.sumAttr[k];
            outNode.sumAttrInter[k] += node.sumAttrInter[k];
          }
        }
      }

    //attribute processign for Haar per direction in the interval [i, i2[
    if (haarFlag) {
      // TODO: simplify

      struct HaarNode {
        int32_t attr[3];
        int32_t attrInter[3];
      } haarNode[4];

      int occupancy = outNode.occupancy;
      int occupancyNext = 0;
      // z transform
      for (int nodePairIdx = 0; nodePairIdx < 4; ++nodePairIdx) {
        auto pairOccupancy = 3 & occupancy >> nodePairIdx * 2;
        if (3 == pairOccupancy) {
          // pair of node
          const auto& nodeA = *traversal.nodeIt[2 * nodePairIdx];
          const auto& nodeB = *traversal.nodeIt[2 * nodePairIdx + 1];
          auto& hN = haarNode[nodePairIdx];

          for (int k = 0; k < numAttrs; k++) {
            auto temp = nodeB.sumAttr[k] - nodeA.sumAttr[k];
            hN.attr[k] = nodeA.sumAttr[k] + (temp >> 1);
            temp = nodeB.sumAttrInter[k] - nodeA.sumAttrInter[k];
            hN.attrInter[k] = nodeA.sumAttrInter[k] + (temp >> 1);
          }
        } else if (pairOccupancy) {
          // single node
          const auto& node = *traversal.nodeIt[2 * nodePairIdx + (pairOccupancy >> 1)];
          auto& hN = haarNode[nodePairIdx];

          for (int k = 0; k < numAttrs; k++) {
            hN.attr[k] = node.sumAttr[k];
            hN.attrInter[k] = node.sumAttrInter[k];
          }
        }
        occupancyNext |= (pairOccupancy != 0) << nodePairIdx;
      }
      std::swap(occupancy, occupancyNext);
      occupancyNext = 0;
      // y transform
      for (int nodePairIdx = 0; nodePairIdx < 2; ++nodePairIdx) {
        auto pairOccupancy = 3 & occupancy >> nodePairIdx * 2;
        if (3 == pairOccupancy) {
          // pair of node
          const auto& nodeA = haarNode[2 * nodePairIdx];
          const auto& nodeB = haarNode[2 * nodePairIdx + 1];
          auto& hN = haarNode[nodePairIdx];

          for (int k = 0; k < numAttrs; k++) {
            auto temp = nodeB.attr[k] - nodeA.attr[k];
            hN.attr[k] = nodeA.attr[k] + (temp >> 1);
            temp = nodeB.attrInter[k] - nodeA.attrInter[k];
            hN.attrInter[k] = nodeA.attrInter[k] + (temp >> 1);
          }
        } else if (pairOccupancy) {
          // single node
          const auto& node = haarNode[2 * nodePairIdx + (pairOccupancy >> 1)];
          auto& hN = haarNode[nodePairIdx];

          for (int k = 0; k < numAttrs; k++) {
            hN.attr[k] = node.attr[k];
            hN.attrInter[k] = node.attrInter[k];
          }
        }
        occupancyNext |= (pairOccupancy != 0) << nodePairIdx;
      }
      std::swap(occupancy, occupancyNext);
      // x transform
      for (int nodePairIdx = 0; nodePairIdx < 1; ++nodePairIdx) {
        auto pairOccupancy = 3 & occupancy >> nodePairIdx * 2;
        if (3 == pairOccupancy) {
          // pair of node
          const auto& nodeA = haarNode[2 * nodePairIdx];
          const auto& nodeB = haarNode[2 * nodePairIdx + 1];
          auto& hN = outNode;

          for (int k = 0; k < numAttrs; k++) {
            auto temp = nodeB.attr[k] - nodeA.attr[k];
            hN.sumAttr[k] = nodeA.attr[k] + (temp >> 1);
            temp = nodeB.attrInter[k] - nodeA.attrInter[k];
            hN.sumAttrInter[k] = nodeA.attrInter[k] + (temp >> 1);
          }
        } else {
          // single node
          const auto& node = haarNode[2 * nodePairIdx + (pairOccupancy >> 1)];
          auto& hN = outNode;

          for (int k = 0; k < numAttrs; k++) {
            hN.sumAttr[k] = node.attr[k];
            hN.sumAttrInter[k] = node.attrInter[k];
          }
        }
      }
    } // end Haar attributes

    // TODO: local QP
    outNode.qp = traversal.nodeIt[0]->qp;

    weightsOut.push_back(outNode);

    traversal.next();
  }

  // number of nodes in next level
  return weightsOut.size();
}

//============================================================================
// Core transform process (for encoder)
template<bool haarFlag, int numAttrs>
inline void
uraht_rso_process_encoder(
  const AttributeParameterSet& aps,
  const EncoderAttributeParams& attrEncParams,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  int numPoints,
  int64_t* positions,
  attr_t* attributes,
  const attr_t* attributes_mc,
  PCCResidualsEncoder& encoder,
  int fracBits,
  int clipMax,
  point_t blockStart,
  point_t blockSizeMinusOne,
  BlockBoundaries* blockBoundaries,
  BlockRefBoundaries* blockRefBoundaries)
{
  const RahtPredictionParams& rahtPredParams = aps.rahtPredParams;
  // coeff entropy coder

  // --------- ascend tree per depth  -----------------
  // create leaf nodes
  int regionQpShift = 4;
  std::vector<UrahtNodeRSOEncoder>  nodesHf;
  std::vector <std::vector<UrahtNodeRSOEncoder>> nodesLfStack;

  nodesLfStack.emplace_back();
  nodesLfStack.back().reserve(numPoints);
  auto nodesLfRef = &nodesLfStack.back();
  auto attr = attributes;
  auto attrPredictor = attributes_mc;

  for (int i = 0; i < numPoints; i++) {
    UrahtNodeRSOEncoder node;
    node.pos = positions[i];
    node.weight = 1;
    node.qp = {
      int16_t(pointQpOffsets[i][0] << regionQpShift),
      int16_t(pointQpOffsets[i][1] << regionQpShift) };
    for (int k = 0; k < numAttrs; k++) {
      node.sumAttr[k] = (*attr++);
      node.sumAttrInter[k] = attributes_mc ? (*attrPredictor++) : 0;
    }
    nodesLfRef->emplace_back(node);
  }

  // -----------  bottom up per depth  --------------------
  int numNodes = nodesLfRef->size();
  // for duplicates, skipable if it is known there is no duplicate
  numNodes = reduceUniqueEncoder<haarFlag, numAttrs>(numNodes, nodesLfRef, &nodesHf);
  const bool flagNoDuplicate = nodesHf.size() == 0;
  const bool singleNode = numNodes == 1;
  int numDepth = 0;
  for (int level = 1; numDepth == 0 || numNodes > 1; ++level) {
    // one depth reduction
    nodesLfStack.emplace_back();
    nodesLfStack.back().reserve(numNodes / 3);
    nodesLfRef = &nodesLfStack.back();

    auto& nodesLfRefold = nodesLfStack[nodesLfStack.size() - 2];
    numNodes = reduceDepthRasterEncoder<haarFlag, numAttrs>(level, nodesLfRefold, *nodesLfRef);
    numDepth++;
  }

  // --------- initialize stuff ----------------
  // root node
  auto& rootNode = nodesLfStack.back()[0];
  assert(rootNode.weight == numPoints);

  bool enableACInterPred = aps.inter_prediction_enabled_flag && attributes_mc;
  const bool enableAveragePredictionBlock = rahtPredParams.enable_average_prediction && enableACInterPred;
  const bool CCRPEnabled = rahtPredParams.cross_component_residual_prediction_flag;
  const int maxlevelCCRPenabled = rahtPredParams.chroma_pred_mode_layer_threshold - 1;

  // reconstruction buffers
  int64_t clipMaxAttrRec = int64_t(clipMax) << kFPFracBits;

  // Prediction buffers
  std::array<int64_t, 3 * numAttrs * 8> SampleDomainBuff;
  std::array<int64_t, 3 * numAttrs * 8> transformBuf;
  std::array<int64_t, numAttrs * 8> predIntraCopyBuff;

  const int numBuffers = (2 + enableACInterPred) * numAttrs;
  const size_t nEltsBufs = 8 * numBuffers;

  int64_t* attrPred = &SampleDomainBuff[8 * numAttrs];
  int64_t* attrPredTransform = &transformBuf[8 * numAttrs];
  int64_t* attrReal = &SampleDomainBuff[0];
  int64_t* attrPredIntra;
  int64_t* attrPredInter;
  const int64_t* attrBestPredIt;
  int64_t* attrPredIntraTransformIt;
  int64_t* attrPredInterTransformIt;
  const int64_t* attrBestPredTransformIt;

  if (enableACInterPred) {
    attrPredInter = attrPred;
    attrPredInterTransformIt = attrPredTransform;
    attrPred += 8 * numAttrs;
    attrPredTransform += 8 * numAttrs;
  }
  if (rahtPredParams.intra_prediction_enabled_flag) {
    attrPredIntra = attrPred;
    attrPredIntraTransformIt = attrPredTransform;
  }

  const bool enableRDOQ = !haarFlag && attrEncParams.useRahtRDOQ;
  const bool enableskipCoeffRDOQ = !haarFlag && attrEncParams.useRahtskipCoeffRDOQ;

  // -------------- descend tree, loop on depth --------------
  int boundaryCurrParentNodeIdx[3] = {-1, -1, -1};
  for (int levelD = numDepth, isAttRootNode = 1; levelD > 0; /*nop*/) {
    // references
    std::vector<UrahtNodeRSOEncoder>& nodesParent = nodesLfStack[levelD];
    std::vector<UrahtNodeRSOEncoder>& nodesLf = nodesLfStack[levelD - 1];
    levelD--;

    // iterators for the neighbouring parents
    RSO_OneLevelNeighboursTraversal<decltype(nodesParent.cbegin()), decltype(nodesLf.cbegin())>
    neighbours
      (blockRefBoundaries, nodesParent.cbegin(), nodesParent.cend(),
      rahtPredParams.subnode_prediction_enabled_flag, nodesLf.cbegin(), nodesLf.cend(), levelD);

    BlockBoundariesCheckerRSO bBoundsParentChecker(blockStart, blockSizeMinusOne, levelD + 1);

    bool inheritDc = !isAttRootNode;
    isAttRootNode = 0;

    // 420
    const bool is420 = abh.is420 && levelD == 0 && !singleNode && !haarFlag;

    //CCCP Parameters
    int CccpCoeff = 0;
    PCCRAHTComputeCCCP curlevelCccp;

    //CCRP Parameters
    const bool CCRPFlag = !is420 && !haarFlag && CCRPEnabled
      && (!enableACInterPred || levelD <= maxlevelCCRPenabled);
    CCRPFilter ccrpFilter;

    int64_t boundaryPrevRowIdx[3] = {0, 0, 0};

    OneLevelRAHTNodesTraversal<decltype(nodesLf.begin())>
      traversal(levelD, nodesLf.begin(), nodesLf.end());
    for (auto parentIt = nodesParent.begin();
      parentIt < nodesParent.end();
      parentIt++) {

      traversal.determineOccupancyAndPosition();
      assert(traversal.nodeOccupancy == parentIt->occupancy);
      assert(traversal.nodePos == parentIt->pos);

      std::fill_n(SampleDomainBuff.begin(), nEltsBufs, 0);
      std::fill_n(transformBuf.begin(), nEltsBufs, 0);
      using WeightsType =
        typename std::conditional<haarFlag, bool, int64_t>::type;
      WeightsType weights[8 + 8 + 8 + 8 + 24] = {};

      int64_t interPredictor[8 * 3];
      int childTable[8] = { };

      int64_t sqrtweightsbuf[8] = { 1 << kFISqrtFracBits, 1 << kFISqrtFracBits, 1 << kFISqrtFracBits, 1 << kFISqrtFracBits, 1 << kFISqrtFracBits, 1 << kFISqrtFracBits, 1 << kFISqrtFracBits, 1 << kFISqrtFracBits }; //value for weigth is 1
      int64_t normsqrtsbuf[8] = {};
      bool skipInverseTransform = !haarFlag;

      Qps nodeQp[8] = {};
      uint8_t occupancy = 0;

      // generate weights, occupancy mask, and fwd transform buffers
      // for all siblings of the current node.
      const int nodeCnt = popcnt(parentIt->occupancy);

      for (int t = 0, j0 = 0; t < nodeCnt; t++, j0++) {
        while (!(traversal.nodeOccupancy >> j0 & 1)) ++j0;
        int nodeIdx = j0;
        auto& node = *traversal.nodeIt[j0];
        childTable[t] = nodeIdx;
        weights[nodeIdx] = node.weight;
        nodeQp[nodeIdx][0] = node.qp[0] >> regionQpShift;
        nodeQp[nodeIdx][1] = node.qp[1] >> regionQpShift;

        // inter predictor
        if (attributes_mc) {
          auto pred = &interPredictor[nodeIdx];
          if (haarFlag) {
            for (int k = 0; k < numAttrs; k++)
              pred[8 * k] = fpExpand<kFPFracBits>(int64_t(node.sumAttrInter[k]));
          }
          else {
            int64_t w = weights[nodeIdx];
            int shift = (15 - kFPFracBits + 1 >> 1) + 5 * (1 + (w > 1024) + (w > 1048576));
            int64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            int64_t divisor = fpReduce<kFPFracBits>(((int64_t(1) << 30) >> shift) * rsqrtWeight);
            divisor = fpReduce<kFPFracBits>((divisor >> shift) * rsqrtWeight);

            for (int k = 0; k < numAttrs; k++)
              pred[8 * k] = (int64_t(node.sumAttrInter[k]) * divisor) >> (30 - kFPFracBits);
          }
        }

        occupancy |= 1 << nodeIdx;

        for (int k = 0; k < numAttrs; k++) {
          attrReal[8 * k + nodeIdx] = fpExpand<kFPFracBits>(int64_t(node.sumAttr[k]));
        }
      }

      // already set during reduction
      assert(occupancy == parentIt->occupancy);
      //parentIt->occupancy = occupancy;
      int64_t sumweights = haarFlag ? 0 : parentIt->weight;
      int Wsum = sumweights;
      parentIt->decoded = true;

      using Kernel =
        typename std::conditional<haarFlag, HaarKernel, RahtKernel>::type;
      mkWeightTree<haarFlag>::template apply<false, Kernel>(weights);

      // Inter-level prediction:
      //  - Find the parent neighbours of the current node
      //  - Generate prediction for all attributes into transformIntraBuf
      //  - Subtract transformed coefficients from forward transform
      //  - The transformIntraBuf is then used for reconstruction

      const bool enableIntraPred = rahtPredParams.intra_prediction_enabled_flag && (nodeCnt > 1 || singleNode);
      const bool enableInterPred = enableACInterPred && (nodeCnt > 1 || singleNode);

      // inter prediction
      int voteMode[3] = { 0,0,0 }; // Null, intra, inter
      int voteInterWeight = 1, voteIntraWeight = 1;
      int numCoeffNot0 = 0;
      int numCoeffTotal = 0;

      int boundariesPosMask = bBoundsParentChecker.computeBoundariesPosMask(parentIt->pos);

      if (enableInterPred)
        memcpy(attrPredInter, interPredictor, 8 * numAttrs * sizeof(int64_t));

      if (enableIntraPred) {
        neighbours.getNeighborsMode(occupancy, voteInterWeight,
          voteIntraWeight, voteMode, numCoeffNot0, numCoeffTotal, boundariesPosMask);

        neighbours.template intraDcPred<haarFlag, numAttrs>(
          occupancy, nodesParent.cbegin(), nodesLf.cbegin(),
          attrPredIntra, rahtPredParams, boundariesPosMask, blockRefBoundaries, inheritDc);
      }

      if (haarFlag) {
        std::copy_n(SampleDomainBuff.begin(), nEltsBufs, transformBuf.begin());
        FwdTransformBlock222<haarFlag>
        ::template apply(numBuffers, &transformBuf[0], weights);
      }
      else {
        // normalise coefficients
        for (int n = 0; n < nodeCnt; n++) {
          const int childIdx = childTable[n];
          int64_t w = weights[childIdx];
          if (w > 1) {
            // Summed attribute values
            int shift = (15 - kFPFracBits + 1 >> 1) + 5 * ((w > 1024) + (w > 1048576));
            int64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            for (int k = 0; k < numAttrs; k++) {
              SampleDomainBuff[8 * k + childIdx] = fpReduce<kFPFracBits>((SampleDomainBuff[8 * k + childIdx] >> shift) * rsqrtWeight);
            }

            // Predicted attribute values
            int64_t sqrtWeight = fastIsqrt(weights[childIdx]);
            sqrtweightsbuf[childIdx] = sqrtWeight;

            for (auto buf = &SampleDomainBuff[8 * numAttrs];
              buf < &SampleDomainBuff[8 * numBuffers];
              buf += 8 * numAttrs) {
              for (int k = 0; k < numAttrs; k++)
                buf[8 * k + childIdx] = fpReduce<kFISqrtFracBits>(
                  buf[8 * k + childIdx] * sqrtWeight);
            }
          }
        }
        std::copy_n(SampleDomainBuff.begin(), nEltsBufs, transformBuf.begin());
        FwdTransformBlock222<haarFlag>
        ::template apply(numBuffers, &transformBuf[0], weights);
      } //else normalize

      bool noChooseMode = (haarFlag ? (1 << 2 * levelD) : sumweights) < rahtPredParams.min_weight_for_mode_selection;
      const bool enableAveragePrediction = enableAveragePredictionBlock && noChooseMode && enableIntraPred && enableInterPred;

      if (enableAveragePrediction) {
        int64_t weightIntra = divApprox(voteIntraWeight << kFPFracBits, voteInterWeight + voteIntraWeight, 0);
        int64_t weightInter = (1 << kFPFracBits) - weightIntra;

        for (int t = 0; t < 8 * numAttrs; t++) {
          predIntraCopyBuff[t] = attrPredIntra[t];
          attrPredIntra[t] = fpReduce<kFPFracBits>(attrPredIntra[t] * weightIntra + attrPredInter[t] * weightInter);
          attrPredIntraTransformIt[t] = fpReduce<kFPFracBits>(attrPredIntraTransformIt[t] * weightIntra + attrPredInterTransformIt[t] * weightInter);
          if (haarFlag) {
            attrPredIntra[t] &= kFPIntMask;
            attrPredIntraTransformIt[t] &= kFPIntMask;
          }
        }
      }

      // ---------- determine best prediction mode -------------
      PredMode predMode = PredMode::Null;
      for (int n = 0; n < nodeCnt; n++) {
        auto& childNode = traversal.nodeIt[childTable[n]];
        childNode->Wintra = parentIt->Wintra;
        childNode->mode = parentIt->mode;
      }
      parentIt->childsWintra = parentIt->Wintra;

      // TODO: try avoiding things that could be avoided when RDO will not be used
      //  - lambda computation, ...

      // TODO: shouldn't stepQ be scaled according to fractional bits or internal bitdepth ?
      double stepQ = qpset.quantizers(0)[0].getStepSize() / 65536.;
      double dLambda = stepQ * stepQ * attrEncParams.rahtRDO_lambdaFactor;

      bool skipCoeff = false;
      bool skipForced = false; // set by RDO
      int zeroForcedMask = 0; // set by RDO

      // per-coefficient operations:
      //  - subtract transform domain prediction (encoder)
      //  - write out/read in quantised coefficients
      //  - inverse quantise + add transform domain prediction
      std::array<int64_t, 8 * numAttrs> RecBuf;
      std::array<int64_t, 8 * numAttrs> CoeffRecBuf;
      std::fill_n(RecBuf.begin(), 8 * numAttrs, 0);
      if (numAttrs == 3)
        std::fill_n(CoeffRecBuf.begin(), 8 * numAttrs, 0);
      int nodelvlSum = 0;

      Vec3<int64_t> CCRPcorr = 0;
      Vec3<int64_t> CCCPcorr = 0;

      bool useRDO = false;
      if (enableIntraPred || enableInterPred) {
        if (noChooseMode)
          predMode = enableIntraPred ? PredMode::Intra : PredMode::Null;
        else if (!enableIntraPred)
          predMode = PredMode::Inter;
        else if (!enableInterPred)
          predMode = PredMode::Intra;
        else if (((parentIt->infoParent >> 6) & 1) == 1)
          predMode = parentIt->mode;
        else if (voteMode[2] == 0 && parentIt->mode == PredMode::Intra && nodeCnt <= 3)
          predMode = PredMode::Intra;
        else if (voteMode[1] == 0 && parentIt->mode == PredMode::Inter && nodeCnt <= 2)
          predMode = PredMode::Inter;
        else {
          // RDO should start here
          predMode = PredMode::Null;
          useRDO = true;
        }
      }

      int predCtxMode = 0;
      if (useRDO) {
        int cousinPred = 0;
        int cousinPredW = voteMode[1] + voteMode[2];
        if (cousinPredW)
          cousinPred = (divApprox(6 * voteMode[2], cousinPredW, 5) + 16) >> 5;
        const int extreme = (voteMode[1] == 0) + 2 * (voteMode[2] == 0);

        predCtxMode = attr::getPredCtxMode(nodeCnt, cousinPred, cousinPredW, extreme);
      }

      if (haarFlag && useRDO) {
        int64_t sacIntra = 0;
        int64_t sacInter = 0;
        for (int idxB = 0; idxB < 8; idxB++) {
          if ((idxB == 0 || weights[24 + idxB]) // there is always the DC coefficient (empty blocks are not transformed)
            && !(inheritDc && !idxB)) {  // skip the DC coefficient unless at the root of the tree
            for (int k = 0; k < numAttrs; k++) {
              const int k8idx = 8 * k + idxB;
              sacIntra += abs(transformBuf[k8idx] - attrPredIntraTransformIt[k8idx]);
              sacInter += abs(transformBuf[k8idx] - attrPredInterTransformIt[k8idx]);
            }
          }
        }
        if (14 * sacIntra < 13 * sacInter) {
          predMode = PredMode::Intra;
          useRDO = false;
        }
        else if (14 * sacInter < 13 * sacIntra) {
          predMode = PredMode::Inter;
          useRDO = false;
        }
        if (!useRDO) {
          encoder.encodePredMode(predCtxMode, predMode);
        }
      }

      auto RDO = encoder.arithmeticEncoder.makeRDO(dLambda,
        encoder.getAttrCtx(), encoder.getModeCtx().modeIsIntra[predCtxMode], skipForced, zeroForcedMask,
        /*SampleDomainBuff,*/ transformBuf, skipCoeff, RecBuf, CoeffRecBuf,
        nodelvlSum, numCoeffNot0, numCoeffTotal,
        attrBestPredIt, attrBestPredTransformIt, parentIt->CntNonZero, skipInverseTransform
        );

      if (useRDO) {
        RDO.start();
      }

      bool enableCrCCRP = CCRPFlag && (levelD <= maxlevelCCRPenabled  || CccpCoeff == 0);
      bool enableCCCP = !haarFlag && !is420 && rahtPredParams.cross_chroma_component_prediction_flag && !enableACInterPred;

      const int predModeStart = useRDO ? PredMode::Intra : predMode;
      const int predModeEnd = (useRDO ? PredMode::Inter : predMode) + 1;
      int modeTakenRDO;
      int blockIndex[8];
      for (int predMode = predModeStart; predMode < predModeEnd; ++predMode) {
        if (useRDO) {
          RDO.startAlternative();

          encoder.encodePredMode(predCtxMode, PredMode(predMode));
        }

        skipInverseTransform = skipInverseTransform && inheritDc;
        if (attr::isNull(predMode)) {
          std::fill_n(&SampleDomainBuff[8 * numAttrs], 8 * numAttrs, 0);
          std::fill_n(&transformBuf[8 * numAttrs], 8 * numAttrs, 0);
          attrBestPredIt = &SampleDomainBuff[8 * numAttrs];
          attrBestPredTransformIt = &transformBuf[8 * numAttrs];
        }
        else if (attr::isIntra(predMode)) {
          attrBestPredIt = attrPredIntra;
          attrBestPredTransformIt = attrPredIntraTransformIt;
        }
        else {
          attrBestPredIt = attrPredInter;
          attrBestPredTransformIt = attrPredInterTransformIt;
        }

        // subtract transformed prediction (should skipping DC, but ok)
        if (!attr::isNull(predMode)) {
          for (int t = 0; t < numAttrs * 8; t++)
            transformBuf[t] -= attrBestPredTransformIt[t];
        }

        // ----------scan blocks -------
        int coefftoCode[8][numAttrs] = {};
        int64_t zeroDist2Total = 0; // Distorsion when coeffs are all set to zero
        int64_t recDist2Total = 0; // Distorsion when coeffs are direct quantization of transformed prediction residuals

        for (int idxB = 0; idxB < 8; idxB++) {
          if ((idxB == 0 || weights[24 + idxB]) // there is always the DC coefficient (empty blocks are not transformed)
            && !(inheritDc && !idxB)) {  // skip the DC coefficient unless at the root of the tree

            int64_t quantizedValues[numAttrs];
            blockIndex[nodelvlSum] = idxB;
            auto quantizers = qpset.quantizers(nodeQp[idxB]);

            if (!haarFlag) {
              for (int k = 0; k < numAttrs; k++) {
                const int k8idx = 8 * k + idxB;
                int64_t ChromaPred = 0;
                if (numAttrs == 3) {
                  if (k == 1 && CCRPFlag) {
                    ChromaPred = fpReduce<kCCRPFiltPrecisionbits>(quantizedValues[0] * ccrpFilter.getYCbFilt());
                    transformBuf[k8idx] -= fpExpand<kFPFracBits>(ChromaPred);
                  }
                  if (k == 2) {
                    if (enableCrCCRP) {
                      ChromaPred = fpReduce<kCCRPFiltPrecisionbits>(quantizedValues[0] * ccrpFilter.getYCrFilt());
                    }
                    else if (enableCCCP) {
                      ChromaPred = fpReduce<kCCCPFiltPrecisionbits>(quantizedValues[1] * CccpCoeff);
                    }
                    transformBuf[k8idx] -= fpExpand<kFPFracBits>(ChromaPred);
                  }
                }

                if (is420 && k > 0)
                  transformBuf[k8idx] = 0;

                int64_t coeff = fpReduce<kFPFracBits>(transformBuf[k8idx]);
                zeroDist2Total += coeff * coeff >> (k ? 1 : 0);

                auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
                int64_t Qcoeff = q.quantize(coeff << kFixedPointAttributeShift, (k > 0) + 2 * attr::isInter(predMode), Wsum);
                auto recCoeff = q.scale(Qcoeff, (k > 0) + 2 * attr::isInter(predMode), Wsum);
                quantizedValues[k] = recCoeff + ChromaPred;
                recDist2Total += (coeff - recCoeff) * (coeff - recCoeff) >> (k ? 1 : 0);
              }
            }

            // The RAHT transform

            for (int k = 0; k < numAttrs; k++) {
              const int k8idx = 8 * k + idxB;

              if (is420 && k > 0)
                transformBuf[k8idx] = 0;
              auto coeff = fpReduce<kFPFracBits>(transformBuf[k8idx]);
              assert(coeff <= INT_MAX && coeff >= INT_MIN);
              auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

              quantizedValues[k] = coeff;
              if (!haarFlag) {
                coeff = q.quantize(coeff << kFixedPointAttributeShift, (k > 0) + 2 * attr::isInter(predMode), Wsum); // does nothing for Haar
                quantizedValues[k] = q.scale(coeff, (k > 0) + 2 * attr::isInter(predMode), Wsum);
              }

              //skipInverseTransform = skipInverseTransform && (quantizedValues[k] == 0);
              if (numAttrs == 3)
                CoeffRecBuf[nodelvlSum * numAttrs + k] = quantizedValues[k];
              RecBuf[k8idx] = fpExpand<kFPFracBits>(quantizedValues[k]);
              coefftoCode[nodelvlSum][k] = coeff;
            } //end loop on attributes

            nodelvlSum++;
          }
        }// end of scan block

        numCoeffNot0 += parentIt->infoParent & 7;
        numCoeffTotal += ((parentIt->infoParent >> 3) & 7) + 1;

        bool eligibleSkipCoeff = nodelvlSum >= 1;

        bool doSkipRDO = false;
        auto skipRDO = encoder.arithmeticEncoder.makeRDO(dLambda,
          encoder.getAttrCtx(), numCoeffTotal, parentIt->CntNonZero, numCoeffNot0,
          recDist2Total);
        if (eligibleSkipCoeff) {
          skipCoeff = true;

          int sac = 0;
          for (int c = 0; c < nodelvlSum; c++)
            for (int k = 0; k < numAttrs; k++) {
              skipCoeff = skipCoeff && (coefftoCode[c][k] == 0);
              sac += abs(coefftoCode[c][k]);
            }

          // do not try RDOQ if coefficients are too big
          doSkipRDO = enableskipCoeffRDOQ && !skipCoeff && inheritDc && sac <= numAttrs * nodelvlSum;
          if (doSkipRDO) {
            skipRDO.start();
            // start testing to force skip (quantize all coeff to zero)
            skipRDO.startAlternative();
          }
          if (skipCoeff || doSkipRDO) {
            encoder.encodeZeroBlock(true, nodelvlSum, numCoeffNot0, numCoeffTotal, (parentIt->infoParent >> 6) & 1, enableAveragePrediction);
          }
          if (doSkipRDO) {
            skipRDO.finishAlternative(zeroDist2Total);
            // start testing to not force quantizing for skip
            skipRDO.startAlternative();
          }
          if (!skipCoeff) {
            encoder.encodeZeroBlock(false, nodelvlSum, numCoeffNot0, numCoeffTotal, (parentIt->infoParent >> 6) & 1, enableAveragePrediction);
          }
        }

        if (!skipCoeff) {
          bool foundNotZero = !eligibleSkipCoeff;
          int existsNoZeroInBlock = 0;
          for (int c = 0; c < nodelvlSum; c++) {
            bool doFlagZeroRDO = enableRDOQ && (inheritDc || c > 0);
            auto flagZeroRDO = encoder.arithmeticEncoder.makeRDO(dLambda,
              encoder.getAttrCtx(), numCoeffTotal, foundNotZero, existsNoZeroInBlock,
                parentIt->CntNonZero, numCoeffNot0);
            // coeff entropy coding
            bool flagZero = numAttrs == 3 ? !coefftoCode[c][0] && !coefftoCode[c][1] && !coefftoCode[c][2] : !coefftoCode[c][0];
            int64_t forceZeroDist = 0;
            if (foundNotZero || c < nodelvlSum - 1) {
              doFlagZeroRDO = doFlagZeroRDO && !flagZero;
              if (doFlagZeroRDO) {
                // Do not allow to force quantize last non-zero coefficient
                if (!foundNotZero) {
                  doFlagZeroRDO = false;
                  for (int c2 = c + 1; !doFlagZeroRDO && c2 < nodelvlSum; c2++)
                    doFlagZeroRDO = numAttrs == 3 ? coefftoCode[c2][0] || coefftoCode[c2][1] || coefftoCode[c2][2] : coefftoCode[c2][0];
                }
                int sac = 0;
                for (int k = 0; k < numAttrs; k++) {
                  skipCoeff = skipCoeff && (coefftoCode[c][k] == 0);
                  sac += abs(coefftoCode[c][k]);
                }
                // do not try RDOQ if coefficients are too big
                doFlagZeroRDO = doFlagZeroRDO && sac <= 3;
              }
              if (doFlagZeroRDO) {
                flagZeroRDO.start();
                // start testing to force flagZero (quantize coeff for all components to zero)
                flagZeroRDO.startAlternative();
              }
              if (flagZero || doFlagZeroRDO) {
                encoder.encodeZeroCoeffs(true, numCoeffNot0, numCoeffTotal, blockIndex[c], existsNoZeroInBlock, enableAveragePrediction);
                numCoeffTotal += 4;
              }
              if (doFlagZeroRDO) {
                for (int k = 0; k < numAttrs; k++) {
                  const int k8idx = 8 * k + blockIndex[c];
                  int64_t errCoeff = fpReduce<kFPFracBits>(transformBuf[k8idx] - RecBuf[k8idx]);
                  int64_t errZero = fpReduce<kFPFracBits>(transformBuf[k8idx]);
                  forceZeroDist += errZero * errZero - errCoeff * errCoeff >> (k ? 1 : 0);
                }
                flagZeroRDO.finishAlternative(forceZeroDist);
                // start testing to not force quantizing for flagZero
                flagZeroRDO.startAlternative();
              }
              if (!flagZero/*|| doFlagZeroRDO*/) {
                encoder.encodeZeroCoeffs(false, numCoeffNot0, numCoeffTotal, blockIndex[c], existsNoZeroInBlock, enableAveragePrediction);
                numCoeffTotal += 4;
              }
            } else {
              // no RDO is made on this coefficient, as there is no choise
              // Coefficients are coded anyway and can't be zero (otherwize it should be skip block)
              flagZero = false;
              doFlagZeroRDO = false;
              numCoeffTotal += 4;
            }

            if (!flagZero/*|| doFlagZeroRDO*/) {
              foundNotZero = true;
              parentIt->CntNonZero++;
              numCoeffNot0 += 4;
              if (numAttrs == 3)
                encoder.encode(coefftoCode[c][0], coefftoCode[c][1], coefftoCode[c][2], enableAveragePrediction, numCoeffNot0, numCoeffTotal, is420);
              else
                encoder.encode(coefftoCode[c][0]);
              existsNoZeroInBlock = std::max(1, existsNoZeroInBlock);
              if (std::abs(coefftoCode[c][0]) > 1 || std::abs(coefftoCode[c][1]) > 1 || std::abs(coefftoCode[c][2]) > 1)
                existsNoZeroInBlock = 2;
            }

            if (doFlagZeroRDO) {
              auto resFlagZeroRDO = flagZeroRDO.finishAlternative();
              flagZeroRDO.finish();
              bool zeroForced = !resFlagZeroRDO.second;
              if (zeroForced) {
                zeroForcedMask |= 1 << c;
                recDist2Total += forceZeroDist;
              }
            }
          }
        }

        if (doSkipRDO) {
          auto resSkipRDO = skipRDO.finishAlternative(recDist2Total);
          skipRDO.finish();
          skipCoeff = skipForced = !resSkipRDO.second;
          if (skipForced) {
            // update rec Distortion to be skipped distorsion
            recDist2Total = zeroDist2Total;
          }
        }

        skipInverseTransform = skipInverseTransform && skipCoeff;

        if (useRDO) {
          auto resRDO = RDO.finishAlternative(recDist2Total);
          if (resRDO.second)
            modeTakenRDO = predMode;
        }
      } // end of RDO

      if (useRDO) {
        RDO.finish();
        predMode = PredMode(modeTakenRDO);
      }

      // Update according to RDO

      if (skipForced) {
        for (int c = 0; c < nodelvlSum; ++c) {
          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + blockIndex[c];
            // force coefficients to be 0
            transformBuf[k8idx] = 0;
            RecBuf[k8idx] = 0;
            if (numAttrs == 3)
              CoeffRecBuf[c * numAttrs + k] = 0;
          } //end loop on attributes
        } // end of scan block
      } else if (zeroForcedMask) {
        for (int c = 0; c < nodelvlSum; c++) {
          if ((zeroForcedMask >> c) & 1)
            for (int k = 0; k < numAttrs; k++) {
              const int k8idx = 8 * k + blockIndex[c];
              // force coefficients to be 0
              transformBuf[k8idx] = 0;
              RecBuf[k8idx] = 0;
              if (numAttrs == 3)
                CoeffRecBuf[c * numAttrs + k] = 0;
            } //end loop on attributes
        }
      }

      // Coding ends here
      if ((enableInterPred || enableACInterPred) && nodeCnt > 1) {
        parentIt->child_mode = predMode;
        uint8_t W = predMode == PredMode::Inter ? 0 : 128;
        parentIt->childsWintra = W;
        // store pred mode in child nodes, to determine best mode at next depth
        for (int n = 0; n < nodeCnt; n++) {
          auto& childNode = traversal.nodeIt[childTable[n]];
          childNode->Wintra = W;
          childNode->mode = predMode;
        }
      }

      uint8_t infoForChild = parentIt->CntNonZero + (nodeCnt - inheritDc << 3) + (skipCoeff << 6);
      for (int t = 0; t < 8; t++) {
        if (traversal.nodeOccupancy >> t & 1)
          traversal.nodeIt[t]->infoParent = infoForChild;
      }

      // Transform Domain Pred for Lossless case
      if (haarFlag)
        for (int t = 0; t < numAttrs * 8; t++)
          RecBuf[t] += attrBestPredTransformIt[t];

      if (numAttrs == 3 && !haarFlag) {
        for (int c = 0; c < nodelvlSum; ++c) {
          if (CCRPFlag)
            CoeffRecBuf[c * numAttrs + 1] += fpReduce<kCCRPFiltPrecisionbits>(CoeffRecBuf[c * numAttrs + 0] * ccrpFilter.getYCbFilt());
          if(enableCrCCRP)
            CoeffRecBuf[c * numAttrs + 2] += fpReduce<kCCRPFiltPrecisionbits>(CoeffRecBuf[c * numAttrs + 0] * ccrpFilter.getYCrFilt());
          else if (enableCCCP)
            CoeffRecBuf[c * numAttrs + 2] += fpReduce<kCCCPFiltPrecisionbits>(CoeffRecBuf[c * numAttrs + 1] * CccpCoeff);

          if (CCRPFlag) {
            CCRPcorr[0] += CoeffRecBuf[c * numAttrs + 0] * CoeffRecBuf[c * numAttrs + 0];
            CCRPcorr[1] += CoeffRecBuf[c * numAttrs + 0] * CoeffRecBuf[c * numAttrs + 1];
            CCRPcorr[2] += CoeffRecBuf[c * numAttrs + 0] * CoeffRecBuf[c * numAttrs + 2];
          }
          if (enableCCCP) {
            CCCPcorr[0] += CoeffRecBuf[c * numAttrs + 1] * CoeffRecBuf[c * numAttrs + 2];
            CCCPcorr[1] += CoeffRecBuf[c * numAttrs + 1] * CoeffRecBuf[c * numAttrs + 1];
          }

          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + blockIndex[c];
            RecBuf[k8idx] = fpExpand<kFPFracBits>(CoeffRecBuf[c * numAttrs + k]);
          } //end loop on attributes
        } // end of scan block

        // ------- update CCRP and CCCP ----------
        if (CCRPFlag && !skipInverseTransform)
          ccrpFilter.update(CCRPcorr);

        if (nodeCnt > 1 && enableCCCP && !skipInverseTransform)
          CccpCoeff = curlevelCccp.computeCrossChromaComponentPredictionCoeff(CCCPcorr);
      }



      // -------  handle DC coefficient   ----------
      int64_t PredDC[3] = { 0,0,0 };
      int64_t rsqrtweightsum;
      if (!haarFlag) {
        rsqrtweightsum = fastIrsqrt(sumweights);
        //compute DC of prediction signals
        for (int n = 0; n < nodeCnt; n++) {
          const int childIdx = childTable[n];

          int64_t normSqrtW = sqrtweightsbuf[childIdx] * rsqrtweightsum >> 40;
          normsqrtsbuf[childIdx] = normSqrtW;
          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + childIdx;
            PredDC[k] += fpReduce<kFISqrtFracBits>(normSqrtW * attrBestPredIt[k8idx]);
          }
        }

        // replace DC coefficient with parent minus pred if inheritable
        if (inheritDc)
          for (int k = 0; k < numAttrs; k++)
            RecBuf[8 * k] = parentIt->reconstructedRahtDC[k] - PredDC[k];
      }

      // for Haar
      if (haarFlag && inheritDc) { // replace DC coefficient with parent if
        for (int k = 0; k < numAttrs; k++)
          RecBuf[8 * k] = parentIt->reconstructedAttr[k];
      }

      // ---------- reconstructed attribute above root node   ---------------
      if (!inheritDc && blockBoundaries && boundariesPosMask & 0x38) {
        int64_t DC[3];
        for (int k = 0; k < numAttrs; k++) {
          DC[k] = RecBuf[8 * k];
          if (!haarFlag) {
            DC[k] += PredDC[k];
            DC[k] *= rsqrtweightsum >> 40 - kFPFracBits;
            DC[k] >>= kFPFracBits;
          }
        }
        for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
          if (boundariesPosMask & 8 << faceIdx) {
            assert(levelD + 1 < blockBoundaries->layerEnd[faceIdx].size());
            int64_t rootPos = nodesParent[0].pos;
            // propagate to any upper layer
            for (int L = blockBoundaries->layerEnd[faceIdx].size() - 1; L >= levelD + 1; --L) {
              if (boundaryCurrParentNodeIdx[faceIdx] >= 0) {
                auto& nodeBound = blockBoundaries->nodes[faceIdx][boundaryCurrParentNodeIdx[faceIdx]];
                nodeBound.occupancy = 1;
                nodeBound.CntNonZero = 1;
                nodeBound.child_mode = predMode;
                nodeBound.childsWintra = predMode == PredMode::Inter ? 0 : 128;
              }

              blockBoundaries->nodes[faceIdx].emplace_back(rootPos & mask_pos_boundary_rso[faceIdx]);
              auto& bNode = blockBoundaries->nodes[faceIdx].back();
              bNode.mode = predMode;
              for (int k = 0; k < numAttrs; k++)
                bNode.reconstructedAttr[k] = DC[k];
              blockBoundaries->layerEnd[faceIdx][L] = blockBoundaries->nodes[faceIdx].size();

              ++boundaryCurrParentNodeIdx[faceIdx];
            }
          }
        }
      }

      // ------------ inverse transform -----------
      if (haarFlag) {
        InvTransformBlock222<haarFlag>
          ::template apply<numAttrs>(&RecBuf[0], weights);
      }
      else {
        if (skipInverseTransform) {
          for (int n = nodeCnt - 1; n >= 0; n--) {
            int childIdx = childTable[n];
            int64_t sq = normsqrtsbuf[childIdx];
            for (int k = 0; k < numAttrs; k++)
              RecBuf[8 * k + childIdx] = fpReduce<kFISqrtFracBits>(sq * RecBuf[8 * k]);
          }
        }
        else {
          InvTransformBlock222<haarFlag>
          ::template apply<numAttrs>(&RecBuf[0], weights);
        }
      }

      // ------------- add DC domain lossy best predictor --------------------
      if (!haarFlag) {
        for (int n = 0; n < nodeCnt; n++) {
          const int childIdx = childTable[n];
          auto& childNode = traversal.nodeIt[childIdx];

          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + childIdx;
            RecBuf[k8idx] += attrBestPredIt[k8idx];
            RecBuf[k8idx] = RecBuf[k8idx] >= 0 ? RecBuf[k8idx] : 0;
            childNode->reconstructedRahtDC[k] = RecBuf[k8idx];
          }
        }
      }

      // ------------- compute  statistics for a posteriror mode  -------------
      int64_t numeratorLMS = 0;
      int64_t denominatorLMS = 0;
      if (enableAveragePrediction) {
        for (int n = 0; n < nodeCnt; n++) {
          const int childIdx = childTable[n];

          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + childIdx;
            int64_t InterIntra = fpReduce<kFPFracBits>(attrPredInter[k8idx] - predIntraCopyBuff[k8idx]);
            int64_t InterDec = fpReduce<kFPFracBits>(attrPredInter[k8idx] - RecBuf[k8idx]);
            numeratorLMS += InterDec * InterIntra;
            denominatorLMS += InterIntra * InterIntra;
          }
        }
      }


      // --------- obtain child mean attribute values ---------------
      for (int n = 0; n < nodeCnt; ++n) {
        const int nodeIdx = childTable[n];
        auto& childNode = traversal.nodeIt[nodeIdx];

        //lossy scale values for next level
        if (!haarFlag) {
          if (weights[nodeIdx] > 1) {
            uint64_t w = weights[nodeIdx];
            int shift = (15 - kFPFracBits + 1 >> 1) + 5 * ((w > 1024) + (w > 1048576));
            int64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            for (int k = 0; k < numAttrs; k++) {
              const int k8idx = 8 * k + nodeIdx;
              RecBuf[k8idx] = fpReduce<kFPFracBits>((RecBuf[k8idx] >> shift) * rsqrtWeight);
            }
          }
        }

        for (int k = 0; k < numAttrs; k++) {
          const int k8idx = 8 * k + nodeIdx;
          if (!haarFlag)
            RecBuf[k8idx] = RecBuf[k8idx] <= clipMaxAttrRec ? RecBuf[k8idx] : clipMaxAttrRec;
          childNode->reconstructedAttr[k] = RecBuf[8 * k + nodeIdx];
        }
      }


      // --------- determine a posteriori best mode based on decoded block ---------
      if (enableAveragePrediction) {
        int mu = 64;
        if (denominatorLMS > 0) {
          numeratorLMS = std::max(int64_t(0), std::min(denominatorLMS, numeratorLMS));
          mu = divApprox(128 * numeratorLMS, denominatorLMS, 0);
        }

        auto aposterioriMode = (mu >= 64)  ? PredMode::Intra : PredMode::Inter;
        for (int n = 0; n < nodeCnt; n++) {
          auto& childNode = traversal.nodeIt[childTable[n]];
          childNode->Wintra = mu;
          childNode->mode = aposterioriMode;
        }
        parentIt->childsWintra = mu;
      }

      // ----------  nodes at boundary ---------------
      if (blockBoundaries && boundariesPosMask & 0x38) {
        // current parent node may have childs added from here
        for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
          if (boundariesPosMask & 8 << faceIdx) {
            auto& nodeBound = blockBoundaries->nodes[faceIdx][boundaryCurrParentNodeIdx[faceIdx]];
            nodeBound.occupancy = occupancy;
            nodeBound.CntNonZero = parentIt->CntNonZero;
            nodeBound.mode = parentIt->mode;
            nodeBound.child_mode = parentIt->child_mode;
            nodeBound.childsWintra = parentIt->childsWintra;

            static constexpr int64_t maskRowIdx[3] =
              { RSO_RAHT::maskX, RSO_RAHT::maskX, RSO_RAHT::maskY };

            auto rowPos = traversal.nodePos & maskRowIdx[faceIdx];
            if (rowPos != boundaryPrevRowIdx[faceIdx]) {
              blockBoundaries->flushNextRow(faceIdx);
              boundaryPrevRowIdx[faceIdx] = rowPos;
            }
          }
        }

        for (int n = 0; n < nodeCnt; ++n) {
          const int nodeIdx = childTable[n];
          auto& childNode = traversal.nodeIt[nodeIdx];

          static constexpr int maskFaces[3][2] = {
            // Z
            { 1 << 1 | 1 << 3 | 1 << 5 | 1 << 7,
              1 << 0 | 1 << 2 | 1 << 4 | 1 << 6 },
            // Y
            { 1 << 2 | 1 << 3 | 1 << 6 | 1 << 7,
              1 << 0 | 1 << 1 | 1 << 4 | 1 << 5 },
            // X
            { 1 << 4 | 1 << 5 | 1 << 6 | 1 << 7,
              1 << 0 | 1 << 1 | 1 << 2 | 1 << 3 }
          };

          int64_t pos = childNode->pos;
          for (int faceIdx = 0; faceIdx < 3; ++faceIdx)
            if (boundariesPosMask & 8 << faceIdx
                && 1 << nodeIdx & maskFaces[faceIdx][
                  // for non powers of 2 block sizes
                  !(bBoundsParentChecker.highBoundChildLevel & (1LL << faceIdx * RSO_RAHT::numBitsPerDim))]) {
              RahtBoundaryNode* pNode;
              auto facePos = pos & mask_pos_boundary_rso[faceIdx];
              if (nodeIdx & 4 >> (faceIdx > 1)) {
                blockBoundaries->nodesNextRow[faceIdx].emplace_back(facePos);
                pNode = &blockBoundaries->nodesNextRow[faceIdx].back();
              } else {
                blockBoundaries->nodes[faceIdx].emplace_back(facePos);
                pNode = &blockBoundaries->nodes[faceIdx].back();
              }
              pNode->mode = childNode->mode;
              for (int k = 0; k < numAttrs; k++)
                pNode->reconstructedAttr[k] = RecBuf[8 * k + nodeIdx];
            }
        }

        for (int faceIdx = 0; faceIdx < 3; ++faceIdx)
          if (boundariesPosMask & (8 << faceIdx)) {
            ++boundaryCurrParentNodeIdx[faceIdx];
          }
      }
      traversal.next();
      if (rahtPredParams.intra_prediction_enabled_flag)
        neighbours.next();
    }

    if (blockBoundaries) {
      for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
        blockBoundaries->flushNextRow(faceIdx);
        blockBoundaries->layerEnd[faceIdx][levelD] = blockBoundaries->nodes[faceIdx].size();
      }
    }

    nodesParent.clear();
  }

  // -------------- process duplicate points at level 0 --------------
  if (flagNoDuplicate) { // write-back reconstructed attributes

    auto attrOut = attributes;
    for (auto& attrNode : nodesLfStack[0]) {
      for (int k = 0; k < numAttrs; k++) {
        auto v = attrNode.reconstructedAttr[k] + kFPOneHalf >> kFPFracBits;
        *attrOut++ = attr_t(PCCClip<int64_t>(v, 0, clipMax));
      }
    }
    return;
  }
  else {
    throw std::runtime_error("Duplicates are not supported currently");
  }
#if 0 // duplictae code is not fucntional anymore
  // case there are duplicates
  std::swap(reconstructedAttr, reconstructedAttrParent);
  auto reconstructedAttrParentIt = reconstructedAttrParent.cbegin();
  auto attrsHfIt = nodesHf.cbegin();

  std::vector<UrahtNodeRSOEncoder>& nodesLf = nodesLfStack[0];
  for (int i = 0, out = 0, iEnd = nodesLf.size(); i < iEnd; i++) {
    int weight = nodesLf[i].weight;

    // unique points have weight = 1
    if (weight == 1) {
      for (int k = 0; k < numAttrs; k++)
        reconstructedAttr[out++] = *reconstructedAttrParentIt++;
      continue;
    }

    // duplicates
    Qps nodeQp = {
      nodesLf[i].qp[0] >> regionQpShift,
      nodesLf[i].qp[1] >> regionQpShift};

    int64_t attrSum[3];
    int64_t attrRecDc[3];
    int64_t sqrtWeight = fastIsqrt(uint64_t(weight));

    int64_t sumCoeff = 0;
    for (int k = 0; k < numAttrs; k++) {
      attrSum[k] = int64_t(nodesLf[i].sumAttr[k]);
      attrRecDc[k] = *reconstructedAttrParentIt++;
      if (!haarFlag) {
        attrRecDc[k] = fpReduce<kFISqrtFracBits>(
          attrRecDc[k] * sqrtWeight);
      }
    }

    int64_t rsqrtWeight;
    for (int w = weight - 1; w > 0; w--) {
      RahtKernel kernel(w, 1);
      HaarKernel haarkernel(w, 1);
      int shift = (15 - kFPFracBits + 1 >> 1) + 5 * ((w > 1024) + (w > 1048576));
      rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;

      auto quantizers = qpset.quantizers(nodeQp);
      for (int k = 0; k < numAttrs; k++) {
        auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

        int64_t transformBuf[2];

        // invert the initial reduction (sum)
        // NB: read from (w-1) since left side came from attrsLf.
        transformBuf[1] = fpExpand<kFPFracBits>(int64_t(attrsHfIt[w - 1].sumAttr[k]));
        if (haarFlag) {
          attrSum[k] -= transformBuf[1] >> 1;
          transformBuf[1] += attrSum[k];
          transformBuf[0] = attrSum[k];
        } else {
          attrSum[k] -= transformBuf[1];
          transformBuf[0] = attrSum[k];

          // NB: weight of transformBuf[1] is by construction 1.
          transformBuf[0] = fpReduce<kFPFracBits>(
            (transformBuf[0] >> shift) * rsqrtWeight);
        }

        if (haarFlag) {
          haarkernel.fwdTransform(transformBuf[0], transformBuf[1]);
        } else {
          kernel.fwdTransform(transformBuf[0], transformBuf[1]);
        }

        auto coeff = fpReduce<kFPFracBits>(transformBuf[1]);
        assert(coeff <= INT_MAX && coeff >= INT_MIN);
        *coeffBufItK[k]++ = coeff = q.quantize(coeff << kFixedPointAttributeShift); // quantize does nothing for Haar
        transformBuf[1] = fpExpand<kFPFracBits>(q.scale(coeff));
        sumCoeff += std::abs(q.quantize(coeff << kFixedPointAttributeShift)); // quantize does nothing for Haar

        // inherit the DC value
        transformBuf[0] = attrRecDc[k];

        if (haarFlag) {
          haarkernel.invTransform(transformBuf[0], transformBuf[1]);
        } else {
          kernel.invTransform(transformBuf[0], transformBuf[1]);
        }

        attrRecDc[k] = transformBuf[0];
        reconstructedAttr[out + w * numAttrs + k] = transformBuf[1];
        if (w == 1)
          reconstructedAttr[out + k] = transformBuf[0];
      }
    }

    attrsHfIt += (weight - 1);
    out += weight * numAttrs;
  }

  // write-back reconstructed attributes
  assert(reconstructedAttr.size() == numAttrs * numPoints);
  auto attrOut = attributes;
  for (auto attr : reconstructedAttr) {
    auto v = attr + kFPOneHalf >> kFPFracBits;
    *attrOut++ = PCCClip(v, 0, std::numeric_limits<attr_t>::max());
  }
#endif
}

//============================================================================
/*
 * RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 *
 * Outputs:
 * weights = list of 'voxelCount' weights associated with each transform coefficient
 * coefficients = quantized transformed attributes array, in column-major order
 * binaryLayer = binary layer where each coefficient was generated
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalTransformRSO(
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
  RAHT::BlockRefBoundaries* blockRefBoundaries)
{
  if (aps.lossless_flag || attribCount != 3)
    abh.is420 = false;

  int clipMax = (1 << desc.bitdepth) - 1 << desc.internalBitdepth - desc.bitdepth;
  switch (attribCount) {
  case 3:
    if (!aps.lossless_flag)
      uraht_rso_process_encoder<false,3>(
        aps, attrEncParams, abh,qpset, pointQpOffsets, voxelCount, positions,
        attributes, attributes_mc, encoder,  desc.internalBitdepth - desc.bitdepth, clipMax, blockStart,
        blockSizeMinusOne, blockBoundaries, blockRefBoundaries);
    else
      uraht_rso_process_encoder<true,3>(
        aps, attrEncParams, abh, qpset, pointQpOffsets, voxelCount, positions,
        attributes, attributes_mc, encoder, desc.internalBitdepth - desc.bitdepth, clipMax, blockStart,
        blockSizeMinusOne, blockBoundaries, blockRefBoundaries);
    break;
  default:
    throw std::runtime_error("attribCount != 3 not tested yet");
  }
}

//============================================================================

}  // namespace pcc
