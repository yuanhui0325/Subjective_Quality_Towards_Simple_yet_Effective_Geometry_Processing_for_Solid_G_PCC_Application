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

using pcc::attr::PredMode;
using namespace pcc::RAHT;

namespace pcc {

#if defined _MSC_VER && defined _DEBUG
// These variables are not used but must be instanciated for MSVC in debug mode
int UrahtNodeDecoderHaar::weight;
std::array<int16_t, 2> UrahtNodeDecoderHaar::qp;
int64_t UrahtNodeDecoderHaar::reconstructedRahtDC[3];
#endif

//============================================================================
// remove any non-unique leaves from a level in the uraht tree

template<int numAttrs, typename Node>
int
reduceUniqueDecoder(
  const bool isInter,
  int numNodes,
  std::vector<Node>* weightsIn,
  std::vector<Node>* weightsOut)
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
    if (isInter)
      for (int k = 0; k < numAttrs; k++)
        (weightsInWrIt - 1)->sumAttrInter[k] += node.sumAttrInter[k];

    weightsOut->push_back(node);
  }

  // number of nodes in next level
  return std::distance(weightsIn->begin(), weightsInWrIt);
}

//============================================================================

template <int numAttrs, typename Node>
struct reduceDepthDecoder_ {
  template <typename It>
  static inline void process(
    const bool isInter, const int level,
    const OneLevelRAHTNodesTraversal<It>& traversal,
    Node& outNode)
  {
    throw std::runtime_error("Not implemented");
  }
};


template <int numAttrs>
struct reduceDepthDecoder_<numAttrs, UrahtNodeDecoder> {
  template <typename It>
  static inline void process(
    const bool isInter, const int level,
    const OneLevelRAHTNodesTraversal<It>& traversal,
    UrahtNodeDecoder& outNode)
  {
    outNode.weight = 0;
    for (int k = 0; k < numAttrs; k++) {
      outNode.sumAttrInter[k] = 0;
    }

    for (int nodeIdx = 0; nodeIdx < 8; ++nodeIdx)
      if (outNode.occupancy & (1 << nodeIdx)) {
        const auto& node = *traversal.nodeIt[nodeIdx];
        outNode.weight += node.weight;

        if (isInter)
          for (int k = 0; k < numAttrs; k++) {
            outNode.sumAttrInter[k] += node.sumAttrInter[k];
      }
    }

    // TODO: local QP
    outNode.qp = traversal.nodeIt[0]->qp;
  }
};

template <int numAttrs>
struct reduceDepthDecoder_<numAttrs, UrahtNodeDecoderHaar> {
  template <typename It>
  static inline void process(
    const bool isInter, const int level,
    const OneLevelRAHTNodesTraversal<It>& traversal,
    UrahtNodeDecoderHaar& outNode)
  {
    if (!isInter)
      return;

    struct HaarNode {
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
          auto temp = nodeB.sumAttrInter[k] - nodeA.sumAttrInter[k];
          hN.attrInter[k] = nodeA.sumAttrInter[k] + (temp >> 1);
        }
      } else if (pairOccupancy) {
        // single node
        const auto& node = *traversal.nodeIt[2 * nodePairIdx + (pairOccupancy >> 1)];
        auto& hN = haarNode[nodePairIdx];

        for (int k = 0; k < numAttrs; k++) {
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
          auto temp = nodeB.attrInter[k] - nodeA.attrInter[k];
          hN.attrInter[k] = nodeA.attrInter[k] + (temp >> 1);
        }
      } else if (pairOccupancy) {
        // single node
        const auto& node = haarNode[2 * nodePairIdx + (pairOccupancy >> 1)];
        auto& hN = haarNode[nodePairIdx];

        for (int k = 0; k < numAttrs; k++) {
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
          auto temp = nodeB.attrInter[k] - nodeA.attrInter[k];
          hN.sumAttrInter[k] = nodeA.attrInter[k] + (temp >> 1);
        }
      } else {
        // single node
        const auto& node = haarNode[2 * nodePairIdx + (pairOccupancy >> 1)];
        auto& hN = outNode;

        for (int k = 0; k < numAttrs; k++) {
          hN.sumAttrInter[k] = node.attrInter[k];
        }
      }
    }
  }
};

template<typename TreeNode, bool haarFlag, int numAttrs>
int
reduceDepthDecoder(
  int level,
  int numNodes,
  std::vector<TreeNode>& weightsIn,
  std::vector<TreeNode>& weightsOut,
  const bool isInter)
{
  // level of child nodes
  const int levelIn = level - 1;

  OneLevelRAHTNodesTraversal<decltype(weightsIn.cbegin())>
    traversal(levelIn, weightsIn.cbegin(), weightsIn.cend());

  while (!traversal.finished()) {
    traversal.determineOccupancyAndPosition();

    typename decltype(weightsOut.begin())::value_type outNode;

    outNode.pos = traversal.nodePos;
    outNode.occupancy = traversal.nodeOccupancy;

    reduceDepthDecoder_<numAttrs, TreeNode>::process(
      isInter, level, traversal, outNode);

    weightsOut.push_back(outNode);

    traversal.next();
  }
  // number of nodes in next level
  return weightsOut.size();
}

//============================================================================
// Core transform process (for decoder)

template<bool haarFlag, int numAttrs, typename Node>
inline void
uraht_process_decoder(
  const AttributeParameterSet& aps,
  AttributeBrickHeader& abh,
  const QpSet& qpset,
  const Qps* pointQpOffsets,
  const int numPoints,
  int64_t* positions,
  attr_t* attributes,
  const attr_t* attributes_mc,
  PCCResidualsDecoder& decoder,
  int clipMax,
  point_t blockStart,
  point_t blockSizeMinusOne,
  BlockBoundaries* blockBoundaries,
  BlockRefBoundaries* blockRefBoundaries)
{
  // coeff entropy coder;
  const RahtPredictionParams& rahtPredParams = aps.rahtPredParams;

  const bool isInter = attributes_mc;

  // --------- ascend tree per depth  -----------------
  // create leaf nodes
  int regionQpShift = 4;
  std::vector<Node> nodesHf;
  std::vector<std::vector<Node>> nodesLfStack;

  nodesLfStack.emplace_back();
  nodesLfStack.back().reserve(numPoints);
  auto nodesLfRef = &nodesLfStack.back();
  auto attrPredictor = attributes_mc;

  for (int i = 0; i < numPoints; i++) {
    Node node;
    node.pos = positions[i];
    if (!haarFlag) {
      node.weight = 1;
      node.qp = {
        int16_t(pointQpOffsets[i][0] << regionQpShift),
        int16_t(pointQpOffsets[i][1] << regionQpShift)};
    }
    if (isInter)
      for (int k = 0; k < numAttrs; k++)
        node.sumAttrInter[k] = *attrPredictor++;
    nodesLfRef->emplace_back(node);
  }

  // -----------  bottom up per depth  --------------------
  int numNodes = nodesLfRef->size();
#if 0 // duplicates won't work anymore with haar Nodes
      // Do we still need that with GeS-TM, if yes,
      // TODO: fix duplicates
  // for duplicates, skipable if it is known there is no duplicate
  numNodes =
    reduceUniqueDecoder<numAttrs>(isInter, numNodes, nodesLfRef, &nodesHf);
#endif

  const bool flagNoDuplicate = nodesHf.size() == 0;
  const bool singleNode = numNodes == 1;
  int numDepth = 0;
  for (int level = 1; numDepth == 0 || numNodes > 1; ++level) {
    // one depth reduction
    nodesLfStack.emplace_back();
    nodesLfStack.back().reserve(numNodes / 3);
    nodesLfRef = &nodesLfStack.back();

    auto& nodesLfRefold = nodesLfStack[nodesLfStack.size() - 2];
    numNodes = reduceDepthDecoder<Node, haarFlag, numAttrs>(
      level, numNodes, nodesLfRefold, *nodesLfRef, isInter);
    numDepth++;
  }

  // --------- initialize stuff ----------------
  // root node
  auto& rootNode = nodesLfStack.back()[0];
  if (!haarFlag)
    assert(rootNode.weight == numPoints);

  const bool enableACInterPred = aps.inter_prediction_enabled_flag && isInter;
  const bool enableAveragePredictionBlock = rahtPredParams.enable_average_prediction && enableACInterPred;
  const bool CCRPEnabled = rahtPredParams.cross_component_residual_prediction_flag;
  const int maxlevelCCRPenabled = rahtPredParams.chroma_pred_mode_layer_threshold - 1;

  // reconstruction buffers
  int64_t clipMaxAttrRec = int64_t(clipMax) << kFPFracBits;

  // Prediction buffers
  int64_t attrPredNull[numAttrs * 8];
  memset(attrPredNull, 0, sizeof(int64_t) * numAttrs * 8);
  int64_t attrPredIntra[numAttrs * 8];
  int64_t predIntraCopyBuff[numAttrs * 8];
  int64_t attrPredInter[numAttrs * 8];
  int64_t* attrBestPred;

  // -------------- descend tree, loop on depth --------------
  int boundaryCurrParentNodeIdx[3] = {-1, -1, -1};
  for (int levelD = numDepth, isAttRootNode = 1; levelD > 0; /*nop*/) {
    // references
    std::vector<Node>& nodesParent = nodesLfStack[levelD];
    std::vector<Node>& nodesLf = nodesLfStack[levelD - 1];
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

    //CCCP parameters
    int CccpCoeff = 0;
    PCCRAHTComputeCCCP curlevelCccp;

    //CCRP parameters
    const bool CCRPFlag = !is420 && !haarFlag && CCRPEnabled && (!enableACInterPred || levelD <= maxlevelCCRPenabled);
    CCRPFilter ccrpFilter;

    int64_t boundaryPrevRowIdx[3] = {0, 0, 0};

    // -------------- loop on nodes of the depth --------------
    OneLevelRAHTNodesTraversal<decltype(nodesLf.begin())>
     traversal(levelD, nodesLf.begin(), nodesLf.end());
    for (auto parentIt = nodesParent.begin();
        parentIt < nodesParent.end();
        parentIt++) {

      traversal.setOccupancyAndPosition(parentIt->occupancy, parentIt->pos);

      typename std::conditional<haarFlag, bool, int64_t>::type weights[8 + 8 + 8 + 8 + 24] = {};
      int64_t interPredictor[8 * 3];
      int childTable[8];

      // generate weights, occupancy mask, and fwd transform
      // for all siblings of the current node.
      Qps nodeQp[8] = {};
      const uint8_t occupancy = parentIt->occupancy;
      const int nodeCnt = popcnt(parentIt->occupancy);

      for (int t = 0, j0 = 0; t < nodeCnt; t++, j0++) {
        while (!(traversal.nodeOccupancy >> j0 & 1)) ++j0;
        int nodeIdx = j0;
        auto& node = *traversal.nodeIt[j0];
        childTable[t] = nodeIdx;
        if (haarFlag) {
          weights[nodeIdx] = true;
        } else {
          weights[nodeIdx] = node.weight;
          nodeQp[nodeIdx][0] = node.qp[0] >> regionQpShift;
          nodeQp[nodeIdx][1] = node.qp[1] >> regionQpShift;
        }

        // inter predictor
        if (isInter) {
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
      }
      parentIt->decoded = true;

      const bool enableIntraPred = rahtPredParams.intra_prediction_enabled_flag && (nodeCnt > 1 || singleNode);
      const bool enableInterPred = enableACInterPred && (nodeCnt > 1 || singleNode);

      int boundariesPosMask = bBoundsParentChecker.computeBoundariesPosMask(parentIt->pos);

      // ------- compute enable intra/iter and neighbor counts ----------
      int voteMode[3] = { 0,0,0 }; // Null, intra, inter
      int voteInterWeight = 1, voteIntraWeight = 1;
      int numCoeffNot0 = 0;
      int numCoeffTotal = 0;

      if (enableIntraPred) {
        neighbours.getNeighborsMode(occupancy, voteInterWeight,
          voteIntraWeight, voteMode, numCoeffNot0, numCoeffTotal, boundariesPosMask);
      }

      // ---------- determine best prediction mode -------------
      PredMode predMode = PredMode::Null;
      for (int n = 0; n < nodeCnt; n++) {
        auto& childNode = traversal.nodeIt[childTable[n]];
        childNode->Wintra = parentIt->Wintra;
        childNode->mode = parentIt->mode;
      }
      parentIt->childsWintra = parentIt->Wintra;

      bool noChooseMode = (haarFlag ? (1 << 2 * levelD) : parentIt->weight) < rahtPredParams.min_weight_for_mode_selection;
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
          int cousinPred = 0;
          int cousinPredW = voteMode[1] + voteMode[2];
          if (cousinPredW)
            cousinPred = (divApprox(6 * voteMode[2], cousinPredW, 5) + 16) >> 5;
          int extreme = (voteMode[1] == 0) + 2 * (voteMode[2] == 0);

          int predCtxMode = attr::getPredCtxMode(nodeCnt, cousinPred, cousinPredW, extreme);
          predMode = decoder.decodePredMode(predCtxMode);
        }

        if (nodeCnt > 1) {
          parentIt->child_mode = predMode;
          uint8_t W = predMode == PredMode::Inter ? 0 : 128;
          for (int n = 0; n < nodeCnt; n++) {
            auto& childNode = traversal.nodeIt[childTable[n]];
            childNode->Wintra = W;
            childNode->mode = predMode;
          }
          parentIt->childsWintra = W;
        }
      }

      const bool enableAveragePrediction = enableAveragePredictionBlock && noChooseMode && enableIntraPred && enableInterPred;

      // ---------- prepare null, intra and inter predictors  -------------
      if (attr::isNull(predMode)) { // null mode
        attrBestPred = attrPredNull;
      }
      if (attr::isIntra(predMode) || enableAveragePrediction) { // intra or average mode
        neighbours.template intraDcPred<haarFlag, numAttrs>(occupancy, nodesParent.begin(), nodesLf.begin(),
          attrPredIntra, rahtPredParams, boundariesPosMask, blockRefBoundaries, inheritDc);
        attrBestPred = attrPredIntra;
      }
      if (attr::isInter(predMode) || enableAveragePrediction) {
        memcpy(attrPredInter, interPredictor, 8 * numAttrs * sizeof(int64_t));
        if (attr::isInter(predMode))
          attrBestPred = attrPredInter;
      }


      // ------- best predictor: sample domain for Harr; denormalize for RAHT    ----------
      int64_t sqrtweightsbuf[8] = {};
      using Kernel = typename std::conditional<haarFlag, HaarKernel, RahtKernel>::type;
      mkWeightTree<haarFlag>::template apply<true, Kernel>(weights);

      if (!attr::isNull(predMode) && haarFlag) { // only if not NULL // HAAR
        if (enableAveragePrediction)
          for (int t = 0; t < 8 * numAttrs; t++)
            predIntraCopyBuff[t] = attrPredIntra[t];

        FwdTransformBlock222<haarFlag>
          ::template apply<numAttrs>(attrBestPred, weights);

        if (enableAveragePrediction)
          FwdTransformBlock222<haarFlag>
            ::template apply<numAttrs>(attrPredInter, weights);
      }

      if (!haarFlag) { // only if RAHT
        // normalise predicted attribute values
        for (int n = 0; n < nodeCnt; n++) {
          const int childIdx = childTable[n];
          int64_t w = weights[childIdx];
          if (w <= 1) {
            sqrtweightsbuf[childIdx] = 1 << kFISqrtFracBits;
            continue;
          }

          int64_t sqrtWeight = fastIsqrt(uint64_t(w));
          sqrtweightsbuf[childIdx] = sqrtWeight;
          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + childIdx;
            attrBestPred[k8idx] = fpReduce<kFISqrtFracBits>(attrBestPred[k8idx] * sqrtWeight);
            if (enableAveragePrediction)
              attrPredInter[k8idx] = fpReduce<kFISqrtFracBits>(attrPredInter[k8idx] * sqrtWeight);
          }
        }

        if (enableAveragePrediction)
          memcpy(predIntraCopyBuff, attrPredIntra, 8 * numAttrs * sizeof(int64_t));
      }

      // -------  compute average predictor    ----------
      if (enableAveragePrediction) {
        int64_t weightIntra = divApprox(voteIntraWeight << kFPFracBits, voteInterWeight + voteIntraWeight, 0);
        int64_t weightInter = (1 << kFPFracBits) - weightIntra;

        for (int t = 0; t < 8 * numAttrs; t++) {
          attrBestPred[t] = fpReduce<kFPFracBits>(attrPredInter[t] * weightInter + attrBestPred[t] * weightIntra);
          if (haarFlag)
            attrBestPred[t] &= kFPIntMask;
        }
      }

      //  ---- decode residual coefficients   ------------
      int coeffCnt = nodeCnt - inheritDc;
      numCoeffNot0 += parentIt->infoParent & 7;
      numCoeffTotal += ((parentIt->infoParent >> 3) & 7) + 1;

      bool eligibleSkipCoeff = coeffCnt >= 1;
      bool skipCoeff = false;
      if (eligibleSkipCoeff)
        skipCoeff = decoder.decodeZeroBlock(coeffCnt, numCoeffNot0, numCoeffTotal, (parentIt->infoParent >> 6) & 1, enableAveragePrediction);

      bool skipInverseTransform = !haarFlag && inheritDc && skipCoeff;

      int coeffPosition[8];
      int idxDecOrder = 0;
      for (int idxB = 0; idxB < 8; idxB++) {
        if (weights[24 + idxB] && (!inheritDc || idxB))
          coeffPosition[idxDecOrder++] = idxB;
      }

      int coefftoCode[8][3] = {0};
      if (!skipCoeff) {
        bool foundNotZero = !eligibleSkipCoeff;
        int existsNoZeroInBlock = 0;
        for (int c = 0; c < coeffCnt; c++) {
          // local coeff entropy decoder
          bool flagZero = false;
          if (foundNotZero || c < coeffCnt - 1)
            flagZero = decoder.decodeZeroCoeffs(numCoeffNot0, numCoeffTotal, coeffPosition[c], existsNoZeroInBlock, enableAveragePrediction);
          numCoeffTotal += 4;

          if (!flagZero)  {
            foundNotZero = true;
            parentIt->CntNonZero++;
            numCoeffNot0 += 4;
            if (numAttrs == 3)
              decoder.decode(coefftoCode[c], enableAveragePrediction, numCoeffNot0, numCoeffTotal, is420);
            else
              coefftoCode[c][0] = decoder.decode();
            existsNoZeroInBlock = std::max(1, existsNoZeroInBlock);
            if (std::abs(coefftoCode[c][0]) > 1 || std::abs(coefftoCode[c][1]) > 1 || std::abs(coefftoCode[c][2]) > 1)
              existsNoZeroInBlock = 2;
          }
        }
      }

      uint8_t infoForChild = parentIt->CntNonZero + (coeffCnt << 3) + (skipCoeff << 6);
      for (int n = 0; n < nodeCnt; n++) {
        auto& childNode = traversal.nodeIt[childTable[n]];
        childNode->infoParent = infoForChild;
      }

      // ----------scan coefficients -------
      int64_t RecBuf[8 * numAttrs] = { 0 };

      if (!skipInverseTransform) {
        int64_t dequantizedCoeff[numAttrs];
        bool enableCrCCRP = CCRPFlag && (levelD <= maxlevelCCRPenabled || CccpCoeff == 0);
        bool enableCCCP = !haarFlag && !is420 && rahtPredParams.cross_chroma_component_prediction_flag && !enableACInterPred;
        Vec3<int64_t> CCRPcorr = 0;
        Vec3<int64_t> CCCPcorr = 0;
        int Wsum = haarFlag ? 0 : parentIt->weight;

        for (int coeffNum = 0; coeffNum < coeffCnt; coeffNum++) {
          int idxB = coeffPosition[coeffNum];
          auto quantizers = qpset.quantizers(nodeQp[idxB]);

          // dequantize coeff
          for (int k = 0; k < numAttrs; k++) {
            auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];
            int64_t coeff = coefftoCode[coeffNum][k];
            if (!haarFlag)
              dequantizedCoeff[k] = q.scale(coeff, (k > 0) + 2 * attr::isInter(predMode), Wsum);
            else
              dequantizedCoeff[k] = coeff;
          }

          // cross-channel component
          if (numAttrs == 3 && !haarFlag) {
            if (CCRPFlag)
              dequantizedCoeff[1] += fpReduce<kCCRPFiltPrecisionbits>(dequantizedCoeff[0] * ccrpFilter.getYCbFilt());
            if (enableCrCCRP)
              dequantizedCoeff[2] += fpReduce<kCCRPFiltPrecisionbits>(dequantizedCoeff[0] * ccrpFilter.getYCrFilt());
            else if (enableCCCP)
              dequantizedCoeff[2] += fpReduce<kCCCPFiltPrecisionbits>(dequantizedCoeff[1] * CccpCoeff);

            if (CCRPFlag) {
              CCRPcorr[0] += dequantizedCoeff[0] * dequantizedCoeff[0];
              CCRPcorr[1] += dequantizedCoeff[0] * dequantizedCoeff[1];
              CCRPcorr[2] += dequantizedCoeff[0] * dequantizedCoeff[2];
            }
            if (enableCCCP) {
              CCCPcorr[0] += dequantizedCoeff[1] * dequantizedCoeff[2];
              CCCPcorr[1] += dequantizedCoeff[1] * dequantizedCoeff[1];
            }
          }

          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + idxB;
            RecBuf[k8idx] = fpExpand<kFPFracBits>(dequantizedCoeff[k]);
            if (haarFlag)  // Transform Domain Pred for Lossless case
              RecBuf[k8idx] += attrBestPred[k8idx];
          }
        } // end scan

        // ------- update CCRP and CCCP ----------
        if (numAttrs == 3 && !haarFlag) {
          if (CCRPFlag)
            ccrpFilter.update(CCRPcorr);

          if (nodeCnt > 1 && enableCCCP)
            CccpCoeff = curlevelCccp.computeCrossChromaComponentPredictionCoeff(CCCPcorr);
        }
      }


      // -------  handle DC coefficient   ----------
      // for RAHT
      int64_t PredDC[3] = { 0, 0, 0 };
      int64_t rsqrtweightsum;
      int64_t normsqrtsbuf[8] = {};
      if (!haarFlag) {
        rsqrtweightsum = fastIrsqrt(parentIt->weight);
        // Compute DC of best pred
        for (int n = 0; n < nodeCnt; n++) {
          int childIdx = childTable[n];
          int64_t normalizedsqrtweight = sqrtweightsbuf[childIdx] * rsqrtweightsum >> 40;
          normsqrtsbuf[childIdx] = normalizedsqrtweight;
          for (int k = 0; k < numAttrs; k++)
            PredDC[k] += fpReduce<kFISqrtFracBits>(normalizedsqrtweight * attrBestPred[8 * k + childIdx]);
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

      // ----------   reconstructed attribute above root node ---------------
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
        int64_t pos = parentIt->pos;
        for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
          if (boundariesPosMask & 8 << faceIdx) {
            assert(levelD + 1 < blockBoundaries->layerEnd[faceIdx].size());
            // propagate to any upper layer
            for (int L = blockBoundaries->layerEnd[faceIdx].size() - 1; L >= levelD + 1; --L) {
              if (boundaryCurrParentNodeIdx[faceIdx] >= 0) {
                auto& nodeBound = blockBoundaries->nodes[faceIdx][boundaryCurrParentNodeIdx[faceIdx]];
                nodeBound.occupancy = 1;
                nodeBound.CntNonZero = 1;
                nodeBound.child_mode = predMode;
                nodeBound.childsWintra = predMode == PredMode::Inter ? 0 : 128;
              }

              blockBoundaries->nodes[faceIdx].emplace_back(pos & mask_pos_boundary_rso[faceIdx]);
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
        ::template apply<numAttrs>(RecBuf, weights);
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
          ::template apply<numAttrs, true>(RecBuf, weights);
        }
      }

      // ------------- add DC domain lossy best predictor --------------------
      if (!haarFlag) {
        for (int n = 0; n < nodeCnt; n++) {
          const int childIdx = childTable[n];
          auto& childNode = traversal.nodeIt[childIdx];

          for (int k = 0; k < numAttrs; k++) {
            const int k8idx = 8 * k + childIdx;
            RecBuf[k8idx] += attrBestPred[k8idx];
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
            int64_t InterIntra = fpReduce<kFPFracBits>((haarFlag ? interPredictor[k8idx] : attrPredInter[k8idx]) - predIntraCopyBuff[k8idx]);
            int64_t InterDec = fpReduce<kFPFracBits>((haarFlag ? interPredictor[k8idx] : attrPredInter[k8idx]) - RecBuf[k8idx]);
            numeratorLMS += InterDec * InterIntra;
            denominatorLMS += InterIntra * InterIntra;
          }
        }
      }

      // --------- obtain child mean attribute values ---------------
      for (int n = 0; n < nodeCnt; n++) {
        const int childIdx = childTable[n];
        auto& childNode = traversal.nodeIt[childIdx];

        //lossy scale values for next level
        if (!haarFlag) {
          uint64_t w = weights[childIdx];
          if (w > 1) {
            int shift = (15 - kFPFracBits + 1 >> 1) + 5 * ((w > 1024) + (w > 1048576));
            int64_t rsqrtWeight = fastIrsqrt(w) >> 40 - shift - kFPFracBits;
            for (int k = 0; k < numAttrs; k++) {
              const int k8idx = 8 * k + childIdx;
              RecBuf[k8idx] = fpReduce<kFPFracBits>((RecBuf[k8idx] >> shift) * rsqrtWeight);
            }
          }
        }

        for (int k = 0; k < numAttrs; k++) {
          const int k8idx = 8 * k + childIdx;
          if (!haarFlag)
            RecBuf[k8idx] = RecBuf[k8idx] <= clipMaxAttrRec ? RecBuf[k8idx] : clipMaxAttrRec;
          childNode->reconstructedAttr[k] = RecBuf[k8idx];
        }
      }

      // --------- determine a posteriori best mode based on decoded block ---------
      if (enableAveragePrediction) {
        int mu = 64;
        if (denominatorLMS > 0) {
          numeratorLMS = std::max(int64_t(0), std::min(denominatorLMS, numeratorLMS));
          mu = divApprox(128 * numeratorLMS, denominatorLMS, 0);
        }

        auto aposterioriMode = (mu >= 64) ? PredMode::Intra : PredMode::Inter;
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
        for (int faceIdx = 0; faceIdx<3; ++faceIdx) {
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
          for (int faceIdx = 0; faceIdx<3; ++faceIdx)
            if (boundariesPosMask & 8 << faceIdx
                && 1 << nodeIdx & maskFaces[faceIdx][
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
        // we've processed all the childs of the current parent node
        for (int faceIdx = 0; faceIdx < 3; ++faceIdx)
          if (boundariesPosMask & 8 << faceIdx)
            ++boundaryCurrParentNodeIdx[faceIdx];
      }

      traversal.next();
      if (rahtPredParams.intra_prediction_enabled_flag)
        neighbours.next();
    } // end loop on nodes of depth

    if (blockBoundaries) {
      for (int faceIdx = 0; faceIdx < 3; ++faceIdx) {
        blockBoundaries->flushNextRow(faceIdx);
        blockBoundaries->layerEnd[faceIdx][levelD] = blockBoundaries->nodes[faceIdx].size();
      }
    }

    nodesParent.clear();
  } // end loop on depth


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

#if 0 // This code won't work with haar Nodes
      // Do we still need that with GeS-TM, if yes,
      // TODO: fix duplicates
      // One way to fix it would be to keep the duplicates count in a sevarate
      // buffer instead of using weight (which is not required for the rest of
      // haar)
  // case there are duplicates
  std::swap(reconstructedAttr, reconstructedAttrParent);
  auto reconstructedAttrParentIt = reconstructedAttrParent.cbegin();

  std::vector<int64_t> attrsHf;
  attrsHf.resize(nodesHf.size()* numAttrs);
  auto attrsHfIt = attrsHf.cbegin();

  std::vector<Node>& nodesLf = nodesLfStack[0];
  for (int i = 0, out = 0, iEnd = nodesLf.size(); i < iEnd; i++) {
    // unique points have weight = 1
    int weight = nodesLf[i].weight;
    if (weight == 1) {
      for (int k = 0; k < numAttrs; k++)
        reconstructedAttr[out++] = *reconstructedAttrParentIt++;
      continue;
    }

    // duplicates
    Qps nodeQp = {
      nodesLf[i].qp[0] >> regionQpShift,
      nodesLf[i].qp[1] >> regionQpShift};

    int64_t attrRecDc[3];
    int64_t sqrtWeight = fastIsqrt(uint64_t(weight));

    for (int k = 0; k < numAttrs; k++) {
      attrRecDc[k] = *reconstructedAttrParentIt++;
      if (!haarFlag) {
        attrRecDc[k] = fpReduce<kFISqrtFracBits>(
          attrRecDc[k] * sqrtWeight);
      }
    }

    if (haarFlag) {
      for (int w = weight - 1; w > 0; w--) {
        HaarKernel haarkernel(w, 1);

        for (int k = 0; k < numAttrs; k++) {

          int64_t transformBuf[2];
          transformBuf[1] = *coeffBufItK[k]++;
          transformBuf[0] = attrRecDc[k]; // inherit the DC value

          haarkernel.invTransform(transformBuf[0], transformBuf[1]);

          attrRecDc[k] = transformBuf[0];
          reconstructedAttr[out + w * numAttrs + k] = transformBuf[1];
          if (w == 1)
            reconstructedAttr[out + k] = transformBuf[0];
        }
      }
    } else {
      for (int w = weight - 1; w > 0; w--) {
        RahtKernel kernel(w, 1);

        auto quantizers = qpset.quantizers(nodeQp);
        for (int k = 0; k < numAttrs; k++) {
          auto& q = quantizers[std::min(k, int(quantizers.size()) - 1)];

          int64_t transformBuf[2];
          int64_t coeff = *coeffBufItK[k]++;
          transformBuf[1] = fpExpand<kFPFracBits>q.scale(coeff));
          transformBuf[0] = attrRecDc[k]; // inherit the DC value

          kernel.invTransform(transformBuf[0], transformBuf[1]);

          attrRecDc[k] = transformBuf[0];
          reconstructedAttr[out + w * numAttrs + k] = transformBuf[1];
          if (w == 1)
            reconstructedAttr[out + k] = transformBuf[0];
        }
      }

    }

    attrsHfIt += (weight - 1) * numAttrs;
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


/*
 * inverse RAHT Fixed Point
 *
 * Inputs:
 * quantStepSizeLuma = Quantization step
 * mortonCode = list of 'voxelCount' Morton codes of voxels, sorted in ascending Morton code order
 * attribCount = number of attributes (e.g., 3 if attributes are red, green, blue)
 * voxelCount = number of voxels
 * coefficients = quantized transformed attributes array, in column-major order
 *
 * Outputs:
 * attributes = 'voxelCount' x 'attribCount' array of attributes, in row-major order
 *
 * Note output weights are typically used only for the purpose of
 * sorting or bucketing for entropy coding.
 */
void
regionAdaptiveHierarchicalInverseTransform(
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
  RAHT::BlockRefBoundaries* blockRefBoundaries)
{
  int clipMax = (1 << desc.bitdepth) - 1 << desc.internalBitdepth - desc.bitdepth;
  switch (attribCount) {
  case 3:
    if (!aps.lossless_flag)
      uraht_process_decoder<false, 3, UrahtNodeDecoder>(
        aps, abh,qpset, pointQpOffsets, voxelCount, positions,
        attributes, attributes_mc, decoder, clipMax,
        blockStart, blockSizeMinusOne, blockBoundaries,
        blockRefBoundaries);
    else
      uraht_process_decoder<true, 3, UrahtNodeDecoderHaar>(
        aps, abh, qpset, pointQpOffsets, voxelCount, positions,
        attributes, attributes_mc, decoder, clipMax,
        blockStart, blockSizeMinusOne, blockBoundaries,
        blockRefBoundaries);
    break;
  default:
    throw std::runtime_error("attribCount != 3 not tested yet");
  }
}

//============================================================================

}  // namespace pcc
