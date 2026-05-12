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

#include "PCCTMC3Encoder.h"

#include <cassert>
#include <limits>
#include <set>
#include <stdexcept>
#include <numeric>

#include "Attribute.h"
#include "AttributeCommon.h"
#include "geometry_params.h"
#include "hls.h"
#include "pointset_processing.h"
#include "geometry.h"
#include "geometry_octree.h"
#include "io_hls.h"
#include "osspecific.h"
#include "partitioning.h"
#include "pcc_chrono.h"
#include "ply.h"
#include "TMC3.h"

namespace pcc {

//============================================================================

PCCPointSet3
getPartition(const PCCPointSet3& src, const std::vector<int32_t>& indexes);

PCCPointSet3 getPartition(
  const PCCPointSet3& src,
  const SrcMappedPointSet& map,
  const std::vector<int32_t>& indexes);

//============================================================================

PCCTMC3Encoder3::PCCTMC3Encoder3() : _frameCounter(-1)
{
  _ctxtMemOctreeGeom.reset(new GeometryOctreeContexts);
}

//----------------------------------------------------------------------------

PCCTMC3Encoder3::~PCCTMC3Encoder3() = default;

//============================================================================

int
PCCTMC3Encoder3::compress(
  const PCCPointSet3& inputPointCloud,
  EncoderParams* params,
  PCCTMC3Encoder3::Callbacks* callback,
  CloudFrame* reconCloud)
{
  // start of frame
  _frameCounter++;

  if (_frameCounter == 0) {
    // Angular predictive geometry coding needs to determine spherical
    // positions.  To avoid quantization of the input disturbing this:
    //  - sequence scaling is replaced by decimation of the input
    //  - any user-specified global scaling is honoured
    _inputDecimationScale = 1.;

    deriveParameterSets(params);
    fixupParameterSets(params);

    _srcToCodingScale = params->codedGeomScale;

    // Determine input bounding box (for SPS metadata) if not manually set
    Box3<int> bbox;
    if (params->autoSeqBbox)
      bbox = inputPointCloud.computeBoundingBox();
    else {
      bbox.min = params->sps.seqBoundingBoxOrigin;
      bbox.max = bbox.min + params->sps.seqBoundingBoxSize - 1;
    }

    if (params->trisoup.alignToNodeGrid && !params->trisoupNodeSizes.empty()) {
      int nodeSizeCommon = lcm_all(
        params->trisoupNodeSizes.begin(),
        params->trisoupNodeSizes.end());
      bbox.min = (bbox.min / nodeSizeCommon) * nodeSizeCommon;
      bbox.max = ((bbox.max + nodeSizeCommon - 1) / nodeSizeCommon) * nodeSizeCommon;
    }

    // Note whether the bounding box size is defined
    // todo(df): set upper limit using level
    bool bboxSizeDefined = params->sps.seqBoundingBoxSize > 0;
    if (!bboxSizeDefined)
      params->sps.seqBoundingBoxSize = (1 << 21) - 1;

    // Then scale the bounding box to match the reconstructed output
    for (int k = 0; k < 3; k++) {
      auto min_k = bbox.min[k];
      auto max_k = bbox.max[k];

      // the sps bounding box is in terms of the conformance scale
      // not the source scale.
      // NB: plus one to convert to range
      min_k = std::round(min_k * params->seqGeomScale);
      max_k = std::round(max_k * params->seqGeomScale);
      params->sps.seqBoundingBoxOrigin[k] = min_k;
      params->sps.seqBoundingBoxSize[k] = max_k - min_k + 1;

      // Compensate the sequence origin such that source point (0,0,0) coded
      // as P_c is reconstructed as (0,0,0):
      //   0 = P_c * globalScale + seqOrigin
      auto gs = Rational(params->sps.globalScale);
      int rem = params->sps.seqBoundingBoxOrigin[k] % gs.numerator;
      rem += rem < 0 ? gs.numerator : 0;
      params->sps.seqBoundingBoxOrigin[k] -= rem;
      params->sps.seqBoundingBoxSize[k] += rem;

      // Convert the origin to coding coordinate system
      _originInCodingCoords[k] = params->sps.seqBoundingBoxOrigin[k];
      _originInCodingCoords[k] /= double(gs);
    }

    // derive local motion parameters
    if (params->gps.inter_prediction_enabled_flag)
      deriveMotionParams(params);

    // Determine the number of bits to signal the bounding box
    params->sps.sps_bounding_box_offset_bits =
      numBits(params->sps.seqBoundingBoxOrigin.abs().max());

    params->sps.sps_bounding_box_size_bits = bboxSizeDefined
      ? numBits(params->sps.seqBoundingBoxSize.abs().max())
      : 0;

    // Allocate storage for attribute contexts and dual motion
    _ctxtMemAttrs.resize(params->sps.attributeSets.size());
    _ctxtMemPredMode.resize(params->sps.attributeSets.size());
    _ctxtMemDualMotion.resize(params->sps.attributeSets.size());
  }

  // placeholder to "activate" the parameter sets
  _sps = &params->sps;
  _gps = &params->gps;
  _aps.clear();
  for (const auto& aps : params->aps) {
    _aps.push_back(&aps);
  }

  // initial geometry IDs
  _tileId = 0;
  _sliceId = 0;
  _sliceOrigin = Vec3<int>{0};
  _firstSliceInFrame = true;

  // Configure output coud
  if (reconCloud) {
    reconCloud->setParametersFrom(*_sps, params->outputFpBits);
    reconCloud->frameNum = _frameCounter;
  }

  // Partition the input point cloud into tiles
  //  - quantize the input point cloud (without duplicate point removal)
  //  - inverse quantize the cloud above to get the initial-sized cloud
  //  - if tile partitioning is enabled,partitioning function produces
  //    vectors tileMaps which map tileIDs to point indexes.
  //    Compute the tile metadata for each partition.
  //  - if not, regard the whole input cloud as a single tile to facilitate
  //    slice partitioning subsequent
  //  todo(df):
  PartitionSet partitions;
  SrcMappedPointSet quantizedInput = quantization(inputPointCloud);

  // write out all parameter sets prior to encoding
  callback->onOutputBuffer(write(*_sps));
  callback->onOutputBuffer(write(*_sps, *_gps));
  for (const auto aps : _aps) {
    callback->onOutputBuffer(write(*_sps, *aps));
  }

  std::vector<std::vector<int32_t>> tileMaps;
  if (params->partition.tileSize) {
    tileMaps = tilePartition(params->partition, quantizedInput.cloud);

    // To tag the slice with the tile id there must be sufficient bits.
    // todo(df): determine sps parameter from the paritioning?
    assert(numBits((int32_t)tileMaps.size() - 1) <= _sps->slice_tag_bits);

    // Default is to use implicit tile ids (ie list index)
    partitions.tileInventory.tile_id_bits = 0;

    // all tileOrigins are relative to the sequence bounding box
    partitions.tileInventory.origin = _sps->seqBoundingBoxOrigin;

    // Get the bounding box of current tile and write it into tileInventory
    partitions.tileInventory.tiles.resize(tileMaps.size());

    // Convert tile bounding boxes to sequence coordinate system.
    // A position in the box must remain in the box after conversion
    // irrispective of how the decoder outputs positions (fractional | integer)
    //   => truncate origin (eg, rounding 12.5 to 13 would not allow all
    //      decoders to find that point).
    //   => use next integer for upper coordinate.
    double gs = Rational(_sps->globalScale);

    for (int t = 0; t < tileMaps.size(); t++) {
      Box3<int32_t> bbox = quantizedInput.cloud.computeBoundingBox(
        tileMaps[t].begin(), tileMaps[t].end());

      auto& tileIvt = partitions.tileInventory.tiles[t];
      tileIvt.tile_id = t;
      for (int k = 0; k < 3; k++) {
        auto origin = std::trunc(bbox.min[k] * gs);
        auto size = std::ceil(bbox.max[k] * gs) - origin + 1;
        tileIvt.tileOrigin[k] = origin;
        tileIvt.tileSize[k] = size;
      }
    }
  } else {
    tileMaps.emplace_back();
    auto& tile = tileMaps.back();
    for (int i = 0; i < quantizedInput.cloud.getPointCount(); i++)
      tile.push_back(i);
  }

  if (partitions.tileInventory.tiles.size() > 1) {
    auto& inventory = partitions.tileInventory;
    assert(inventory.tiles.size() == tileMaps.size());
    std::cout << "Tile number: " << tileMaps.size() << std::endl;
    inventory.ti_seq_parameter_set_id = _sps->sps_seq_parameter_set_id;
    inventory.ti_origin_bits_minus1 =
      numBits(inventory.origin.abs().max()) - 1;

    // The inventory comes into force on the first frame
    inventory.ti_frame_ctr_bits = _sps->frame_ctr_bits;
    inventory.ti_frame_ctr = _frameCounter & ((1 << _sps->frame_ctr_bits) - 1);

    // Determine the number of bits for encoding tile sizes
    int maxValOrigin = 1;
    int maxValSize = 1;
    for (const auto& entry : inventory.tiles) {
      maxValOrigin = std::max(maxValOrigin, entry.tileOrigin.max());
      maxValSize = std::max(maxValSize, entry.tileSize.max() - 1);
    }
    inventory.tile_origin_bits_minus1 = numBits(maxValOrigin) - 1;
    inventory.tile_size_bits_minus1 = numBits(maxValSize) - 1;

    callback->onOutputBuffer(write(*_sps, partitions.tileInventory));
  }

  // Partition the input point cloud
  //  - get the partitial cloud of each tile
  //  - partitioning function produces a list of point indexes, origin and
  //    optional tile metadata for each partition.
  //  - encode any tile metadata
  //  NB: the partitioning method is required to ensure that the output
  //      slices conform to any codec limits.
  //  todo(df): consider requiring partitioning function to sort the input
  //            points and provide ranges rather than a set of indicies.

  // use the largest trisoup node size as a partitioning boundary for
  // consistency between slices with different trisoup node sizes.

  int32_t partitionBoundary = 1;
  if (!params->trisoupNodeSizes.empty())
    partitionBoundary = lcm_all(
      params->trisoupNodeSizes.begin(),
      params->trisoupNodeSizes.end());

  do {
    for (int t = 0; t < tileMaps.size(); t++) {
      const auto& tile = tileMaps[t];
      auto tile_id = partitions.tileInventory.tiles.empty()
        ? 0
        : partitions.tileInventory.tiles[t].tile_id;

      // Get the point cloud of current tile and compute their bounding boxes
      PCCPointSet3 tileCloud = getPartition(quantizedInput.cloud, tile);
      Box3<int32_t> bbox = tileCloud.computeBoundingBox();

      if (params->partition.safeTrisoupPartionning) {
        bbox.min[0] -= (bbox.min[0] % partitionBoundary);
        bbox.min[1] -= (bbox.min[1] % partitionBoundary);
        bbox.min[2] -= (bbox.min[2] % partitionBoundary);
      }
      // Move the tile cloud to coodinate origin
      // for the convenience of slice partitioning
      for (int i = 0; i < tileCloud.getPointCount(); i++)
        tileCloud[i] -= bbox.min;

      // don't partition if partitioning would result in a single slice.
      auto partitionMethod = params->partition.method;
      if (tileCloud.getPointCount() < params->partition.sliceMaxPoints)
        partitionMethod = PartitionMethod::kNone;

      //Slice partition of current tile
      std::vector<Partition> curSlices;
      switch (partitionMethod) {
      case PartitionMethod::kNone:
        curSlices = partitionNone(params->partition, tileCloud, tile_id);
        break;

      case PartitionMethod::kUniformGeom:
        curSlices = partitionByUniformGeom(
          params->partition, tileCloud, tile_id, partitionBoundary);
        break;

      case PartitionMethod::kUniformSquare:
        curSlices = partitionByUniformSquare(
          params->partition, tileCloud, tile_id, partitionBoundary);
        break;

      case PartitionMethod::kOctreeUniform:
        curSlices =
          partitionByOctreeDepth(params->partition, tileCloud, tile_id);
        break;

      case PartitionMethod::kNpoints:
        curSlices = partitionByNpts(params->partition, tileCloud, tile_id);
        break;
      }
      // Map slice indexes to tile indexes(the original indexes)
      for (int i = 0; i < curSlices.size(); i++) {
        for (int p = 0; p < curSlices[i].pointIndexes.size(); p++) {
          curSlices[i].pointIndexes[p] = tile[curSlices[i].pointIndexes[p]];
        }
      }

      partitions.slices.insert(
        partitions.slices.end(), curSlices.begin(), curSlices.end());
    }
    std::cout << "Slice number: " << partitions.slices.size() << std::endl;
  } while (0);

  // Encode each partition:
  //  - create a pointset comprising just the partitioned points
  //  - compress
  for (const auto& partition : partitions.slices) {
    // create partitioned point set
    PCCPointSet3 sliceCloud =
      getPartition(quantizedInput.cloud, partition.pointIndexes);

    PCCPointSet3 sliceSrcCloud =
      getPartition(inputPointCloud, quantizedInput, partition.pointIndexes);

    _sliceId = partition.sliceId;
    _tileId = partition.tileId;
    _sliceOrigin = sliceCloud.computeBoundingBox().min;
    if (!params->partition.fixedSliceOrigin.empty()) {
      int idx = std::min(
        _sliceId, int(params->partition.fixedSliceOrigin.size()) - 1);
      _sliceOrigin[0] = std::min<uint32_t>(
        _sliceOrigin[0], params->partition.fixedSliceOrigin[idx][0]);
      _sliceOrigin[1] = std::min<uint32_t>(
        _sliceOrigin[1], params->partition.fixedSliceOrigin[idx][1]);
      _sliceOrigin[2] = std::min<uint32_t>(
        _sliceOrigin[2], params->partition.fixedSliceOrigin[idx][2]);
    }

    if (params->partition.safeTrisoupPartionning) {
      _sliceOrigin[0] -= (_sliceOrigin[0] % partitionBoundary);
      _sliceOrigin[1] -= (_sliceOrigin[1] % partitionBoundary);
      _sliceOrigin[2] -= (_sliceOrigin[2] % partitionBoundary);
    }

    if (params->trisoup.alignToNodeGrid && !params->trisoupNodeSizes.empty()) {
      int nodeSizeCommon = lcm_all(
        params->trisoupNodeSizes.begin(),
        params->trisoupNodeSizes.end());
      _sliceOrigin = ((_sliceOrigin / nodeSizeCommon) * nodeSizeCommon);
    }

    compressPartition(sliceCloud, sliceSrcCloud, params, callback, reconCloud);
  }

  if (_sps->inter_frame_prediction_enabled_flag) {
    // buffer the current frame for potential use in prediction
    _refFrame = *reconCloud;
  }

  // Apply global scaling to reconstructed point cloud
  if (reconCloud)
    scaleGeometry(
      reconCloud->cloud, _sps->globalScale, reconCloud->outputFpBits);

  return 0;
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::deriveParameterSets(EncoderParams* params)
{
  // fixup extGeomScale in the case that we're coding metres
  if (params->sps.seq_geom_scale_unit_flag == ScaleUnit::kMetre)
    params->extGeomScale = 0;

  if (params->extGeomScale == 0.)
    params->extGeomScale = params->srcUnitLength;

  // Derive the sps scale factor:  The sequence scale is normalised to an
  // external geometry scale of 1.
  //  - Ie, if the user specifies that extGeomScale=2 (1 seq point is equal
  //    to 2 external points), seq_geom_scale is halved.
  //
  //  - In cases where the sequence scale is in metres, the external system
  //    is defined to have a unit length of 1 metre, and srcUnitLength must
  //    be used to define the sequence scale.
  //  - The user may define the relationship to the external coordinate system.
  //
  // NB: seq_geom_scale is the reciprocal of unit length
  params->sps.seqGeomScale = params->seqGeomScale / params->extGeomScale;

  // Global scaling converts from the coded scale to the sequence scale
  // NB: globalScale is constrained, eg 1.1 is not representable
  // todo: consider adjusting seqGeomScale to make a valid globalScale
  // todo: consider adjusting codedGeomScale to make a valid globalScale
  params->sps.globalScale =
    Rational(params->seqGeomScale / params->codedGeomScale);
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::fixupParameterSets(EncoderParams* params)
{
  // fixup parameter set IDs
  params->sps.sps_seq_parameter_set_id = 0;
  params->gps.gps_seq_parameter_set_id = 0;
  params->gps.gps_geom_parameter_set_id = 0;
  for (int i = 0; i < params->aps.size(); i++) {
    params->aps[i].aps_seq_parameter_set_id = 0;
    params->aps[i].aps_attr_parameter_set_id = i;
  }

  // development level / header
  params->sps.profile.main_profile_compatibility_flag = 0;
  params->sps.profile.reserved_profile_compatibility_21bits = 0;
  params->sps.level = 0;

  // constraints
  params->sps.profile.unique_point_positions_constraint_flag = false;
  params->sps.profile.slice_reordering_constraint_flag =
    params->sps.entropy_continuation_enabled_flag;

  // use one bit to indicate frame boundaries
  params->sps.frame_ctr_bits = 1;

  // number of bits for slice tag (tileid) if tiles partitioning enabled
  // NB: the limit of 64 tiles is arbritrary
  params->sps.slice_tag_bits = params->partition.tileSize > 0 ? 6 : 0;

  // slice origin parameters used by this encoder implementation
  params->gps.geom_box_log2_scale_present_flag = true;
  params->gps.gps_geom_box_log2_scale = 0;

  // fixup attribute parameters
  for (auto it : params->attributeIdxMap) {
    //auto& attr_sps = params->sps.attributeSets[it.second];
    auto& attr_aps = params->aps[it.second];
    //auto& attr_enc = params->attr[it.second];

    // this encoder does not (yet) support variable length attributes
    // todo(df): add variable length attribute support
    attr_aps.raw_attr_variable_len_flag = 0;
  }
}


//============================================================================

const AttributeDescription*
findColAttrDesc(const std::vector<AttributeDescription>& attrDescs)
{
  // todo(df): don't assume that there is only one in the sps
  for (const auto& desc : attrDescs) {
    if (desc.attributeLabel == KnownAttributeLabel::kColour)
      return &desc;
  }
  return nullptr;
}

//---------------------------------------------------------------------------
// motion parameter derivation
// Global : Setup the block sizes b
// local : Setup the block sizes based on a predefined mode:
//  1 => Large scale sparse point clouds (eg cat3)
//  2 => Small voxelised point clouds (eg cat2)
//  3  => TriSoup
void
PCCTMC3Encoder3::deriveMotionParams(EncoderParams* params)
{
  auto scaleFactor = params->codedGeomScale;

  // local
  auto& gps = params->gps;
  auto& motionGPS = params->gps.motion;
  auto& motionGEnc = params->motion;
  auto& skipGEnc = params->skip;
  auto presetMode = params->motionPreset;
  int TriSoupSize = 1;
  int TriSoupSizeLog2 = 0;
  double QstepTriSoup = 1.;

  auto desc_col = findColAttrDesc(params->sps.attributeSets);

  // Set Geometry parameters

  // some parameters are already set from command line

  switch (presetMode) {
  case MotionPreset::kNoMotion:
    break;

  case MotionPreset::kMotionPresetForOctree:
    std::cout << "presetMode  for octree dense \n";

    if (motionGPS.motion_block_size < 0)
      motionGPS.motion_block_size = 1 << ilog2(uint32_t(std::max(32, int(std::round(128 * scaleFactor)))));

    if (motionGPS.motion_min_pu_size < 0)
      motionGPS.motion_min_pu_size = motionGPS.motion_block_size >> 1;

    // search parameters

    if (motionGEnc.Amotion0 < 0)
      motionGEnc.Amotion0 = 1;

    if (motionGEnc.lambda < 0)
      motionGEnc.lambda = 2.0 / std::sqrt(scaleFactor);

    if (motionGEnc.d_color_factor < 0)
      motionGEnc.d_color_factor = 1.0 / 8;

    if (motionGEnc.d_geom_factor < 0)
      motionGEnc.d_geom_factor = 1.0;

    if (motionGEnc.decimate < 0)
      motionGEnc.decimate = 7;

    break;

  case MotionPreset::kMotionPresetForTrisoup:
    if (params->gps.trisoup_enabled_flag)
      TriSoupSize = *std::max_element(
        params->trisoupNodeSizes.begin(),
        params->trisoupNodeSizes.end());
    std::cout << "presetMode for Trisoup dense of size " << TriSoupSize << " \n";

    TriSoupSizeLog2 = ilog2(uint32_t(TriSoupSize - 1)) + 1;

    // use fixedPoint approximation of exp2 to avoid lib math implementation
    // dependencies and potential encoding mismatches.
    QstepTriSoup = fpToDouble<32>(fpExp2<32>(fpExpand<32>(params->gbh.trisoup_QP - 12LL) / 6));

    if (motionGPS.motion_block_size < 0)
      motionGPS.motion_block_size = std::min(256, 16 << TriSoupSizeLog2);

    if (motionGPS.motion_min_pu_size < 0)
      motionGPS.motion_min_pu_size = std::max(1 << TriSoupSizeLog2, motionGPS.motion_block_size >> 1);

    // search parameters

    if (motionGEnc.Amotion0 < 0)
      motionGEnc.Amotion0 = 1;

    if (motionGEnc.lambda < 0)
      motionGEnc.lambda = 8 * QstepTriSoup * QstepTriSoup;

    if (motionGEnc.d_color_factor < 0)
      motionGEnc.d_color_factor = 1.0 / 48;

    if (motionGEnc.d_geom_factor < 0)
      motionGEnc.d_geom_factor = 1.0;

    if (motionGEnc.decimate < 0)
      motionGEnc.decimate = 9;

    if (gps.trisoup_early_skip_coding_mode_enabled_flag) {
      gps.trisoup_early_skip_coding_mode_node_size = 4 << TriSoupSizeLog2;
      skipGEnc.subsampleStep = 8;
      skipGEnc.QstepTriSoup = QstepTriSoup;
      skipGEnc.occRate = 2.; // 8. * 0.25; // 8bits intra and inter pred
      skipGEnc.triRate = 3.5; // (3 * 0.25 + 3 * 0.5) + (1 * 0.25 + 2 * 0.5);// 3 edges with one occupied, 1 centro  intra and inter pred
      skipGEnc.lambda = 256. * QstepTriSoup * QstepTriSoup;
      skipGEnc.strength = params->trisoupEarlySkipStrength;
    }

    break;

  default:
    throw std::runtime_error("unsupported motion preset");
  }

  if (desc_col)
    motionGEnc.d_color_factor /= 1 << desc_col->internalBitdepth - desc_col->bitdepth;

  std::cout << "params->sps.geomPreScale : " << scaleFactor << "\n";
  std::cout << "motion_block_size(geom) : " << motionGPS.motion_block_size << "\n";
  std::cout << "motion_min_pu_size(geom) : " << motionGPS.motion_min_pu_size << "\n";

  // set attributes parameters
  // for each attribute

  const int log2slabThickness = ilog2(uint32_t(std::max(
      params->sps.localized_attributes_slab_thickness_minus1,
      params->sps.localized_attributes_slab_block_size_minus1))) + 1;

  for (const auto& it : params->attributeIdxMap) {
    int attrIdx = it.second;
    const auto& attr_sps = params->sps.attributeSets[attrIdx];
    auto& attr_aps = params->aps[attrIdx];
    auto& attr_enc = params->attr[attrIdx];
    const auto& label = attr_sps.attributeLabel;

    if (!attr_aps.dual_motion_field_flag)
      continue;

    auto& motionAPS = attr_aps.motion;
    auto& motionAEnc = attr_enc.motion;

    int QPluma = params->aps[0].init_qp_minus4 + 4;
    double Qstep = fpToDouble<32>(fpExp2<32>(fpExpand<32>(int64_t(QPluma - 12)) / 6));

    switch (presetMode) {
    case MotionPreset::kNoMotion:
      break;

    case MotionPreset::kMotionPresetForOctree:
      if (motionAPS.motion_block_size < 0)
        motionAPS.motion_block_size = motionGPS.motion_block_size;

      if (motionAPS.motion_min_pu_size < 0)
        motionAPS.motion_min_pu_size =
          std::min(motionAPS.motion_block_size, motionGPS.motion_min_pu_size);

      // search parameters

      if (motionAEnc.Amotion0 < 0)
        motionAEnc.Amotion0 = 4;

      if (motionAEnc.lambda < 0)
        motionAEnc.lambda = 0.25 * Qstep * Qstep;

      if (motionAEnc.d_color_factor < 0)
        motionAEnc.d_color_factor = 2.0;

      if (motionAEnc.d_geom_factor < 0)
        motionAEnc.d_geom_factor = 1.0;

      if (motionAEnc.decimate < 0)
        motionAEnc.decimate = 8;

      break;

    case MotionPreset::kMotionPresetForTrisoup:
      if (motionAPS.motion_block_size < 0)
        motionAPS.motion_block_size =
          std::max(1 << TriSoupSizeLog2, motionGPS.motion_block_size / 2);

      if (motionAPS.motion_min_pu_size < 0)
        motionAPS.motion_min_pu_size =
          std::max(1 << TriSoupSizeLog2, std::min(
            motionAPS.motion_block_size, motionGPS.motion_min_pu_size / 2));

      // search parameters

      if (motionAEnc.Amotion0 < 0)
        motionAEnc.Amotion0 = 4;

      if (motionAEnc.lambda < 0)
        motionAEnc.lambda = 0.25 * Qstep * Qstep;

      if (motionAEnc.d_color_factor < 0)
        motionAEnc.d_color_factor = 2.0;

      if (motionAEnc.d_geom_factor < 0)
        motionAEnc.d_geom_factor = 1.0;

      if (motionAEnc.decimate < 0)
        motionAEnc.decimate = 8;

      break;

    default:
      throw std::runtime_error("unsupported motion preset");
    }

    motionAEnc.d_color_factor /= 1 << attr_sps.internalBitdepth - attr_sps.bitdepth;

    if (attr_aps.dual_motion_field_flag
        && ilog2(uint32_t(motionAPS.motion_block_size - 1)) + 1 > log2slabThickness) {
      motionAPS.motion_block_size = std::max(
        motionAPS.motion_min_pu_size, 1 << log2slabThickness);
    }

    std::cout << "motion_block_size(" << label << ") : " << motionAPS.motion_block_size << "\n";
    std::cout << "motion_min_pu_size(" << label << ") : " << motionAPS.motion_min_pu_size << "\n";
  }
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::compressPartition(
  const PCCPointSet3& inputPointCloud,
  const PCCPointSet3& originPartCloud,
  EncoderParams* params,
  PCCTMC3Encoder3::Callbacks* callback,
  CloudFrame* reconCloud)
{
  // geometry compression consists of the following stages:
  //  - prefilter/quantize geometry (non-normative)
  //  - encode geometry (single slice, id = 0)
  //  - recolour

  pointCloud.clear();
  pointCloud = inputPointCloud;

  // apply a custom trisoup node size
  params->gbh.trisoup_node_size = 1;
  if (_gps->trisoup_enabled_flag) {
    int idx = std::min(_sliceId, int(params->trisoupNodeSizes.size()) - 1);
    params->gbh.trisoup_node_size =
      params->trisoupNodeSizes[idx];
  }

  // Offset the point cloud to account for (preset) _sliceOrigin.
  // The new maximum bounds of the offset cloud
  const size_t pointCount = pointCloud.getPointCount();
  for (size_t i = 0; i < pointCount; ++i) {
    pointCloud[i] -= _sliceOrigin;
  }
  // the encoded coordinate system is non-negative.
  Box3<int32_t> bbox = pointCloud.computeBoundingBox();

  // todo(df): don't update maxBound if something is forcing the value?
  // NB: size is max - min + 1
  _sliceBoxWhd = bbox.max + 1;

  // geometry encoding
  attrInterPredParams.compensatedPointCloud.clear();
  attrInterPredParams.compensatedPointCloud.addRemoveAttributes(false, false);
  attrInterPredParams.attributes_mc.clear();
  attrInterPredParams.mvField.clear();

  payload_geom = PayloadBuffer(PayloadType::kGeometryBrick);

  payload_attr = std::vector<PayloadBuffer>(
      params->attributeIdxMap.size(), PayloadType::kAttributeBrick);

  clock_user_geom = pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock>();

  clock_user_attr =
    std::vector<pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock>>(
      params->attributeIdxMap.size());

  _abh = std::vector<AttributeBrickHeader>(params->attributeIdxMap.size());

  attrEncoder =
    std::vector<decltype(makeAttributeEncoder())>(params->attributeIdxMap.size());

  type = params->localized_attributes_encoding
      ? Type::kLocalEncoder
      : Type::kGlobalEncoder;
  currSlabIdx = -1;
  slabBlockSize = {
    _sps->localized_attributes_slab_thickness_minus1 + 1,
    _sps->localized_attributes_slab_block_size_minus1 + 1,
    _sps->localized_attributes_slab_block_size_minus1 + 1
  };

  startAndNumPointsPerSlabBlock.clear();
  origin = _originInCodingCoords + _sliceOrigin;
  targetToSourceScaleFactor = 1.0 / _srcToCodingScale;
  bBoxOrigin = originPartCloud.computeBoundingBox();
  this->params = params;
  this->originPartCloud = &originPartCloud;

  // init

  // some parameters used by attributes may be derived...
  startEncodeGeometryBrick(params);

  for (int i = 0; i < params->attributeIdxMap.size(); ++i) {
    attrEncoder[i] = makeAttributeEncoder();
  }

  int deltaQPhierarchical = 0;
  if (params->hierarchicalRAHTQP.size() >= 2) {
    int frameInGop = _frameCounter % params->randomAccessPeriod;
    bool isIntra = frameInGop == 0;
    int frameInGopStructure = (frameInGop - 1) % (params->hierarchicalRAHTQP.size() - 1);
    deltaQPhierarchical = params->hierarchicalRAHTQP[isIntra ? 0 : frameInGopStructure + 1];
  }

  for (const auto& attributeIdxMap : params->attributeIdxMap) {
    int attrIdx = attributeIdxMap.second;
    const auto& attr_sps = _sps->attributeSets[attrIdx];
    const auto& attr_aps = *_aps[attrIdx];
    const auto& attr_enc = params->attr[attrIdx];

    // todo(df): move elsewhere?
    AttributeBrickHeader& abh = _abh[attrIdx];
    abh.attr_attr_parameter_set_id = attr_aps.aps_attr_parameter_set_id;
    abh.attr_sps_attr_idx = attrIdx;
    abh.attr_geom_slice_id = _sliceId;
    abh.attr_qp_delta_luma = attr_aps.inter_prediction_enabled_flag ? deltaQPhierarchical : 0;
    abh.attr_qp_delta_chroma = attr_aps.inter_prediction_enabled_flag ? deltaQPhierarchical : 0;
    abh.is420 = attr_enc.abh.is420;

    // NB: regionQpOrigin/regionQpSize use the STV axes, not XYZ.
    if (false) {
      abh.qpRegions.emplace_back();
      auto& region = abh.qpRegions.back();
      region.regionOrigin = 0;
      region.regionSize = 0;
      region.attr_region_qp_offset = {0, 0};
      abh.attr_region_bits_minus1 = -1
        + numBits(
          std::max(region.regionOrigin.max(), region.regionSize.max()));
    }
    // Number of regions is constrained to at most 1.
    assert(abh.qpRegions.size() <= 1);

    // TODO: this is not ok if more than one attribute, fix it:
    attrInterPredParams.enableAttrInterPred =
      attr_aps.inter_prediction_enabled_flag && _gbh.slice_inter_prediction_flag;

    auto& ctxtMemAttr = _ctxtMemAttrs.at(abh.attr_sps_attr_idx);
    auto& ctxtMemPredMode = _ctxtMemPredMode.at(abh.attr_sps_attr_idx);
    auto& ctxtMemDualMotion = _ctxtMemDualMotion.at(abh.attr_sps_attr_idx);
    attrEncoder[attrIdx]->startEncode(
      *_sps,*_gps, attr_sps, attr_aps, abh, ctxtMemAttr, ctxtMemPredMode, ctxtMemDualMotion,
      pointCloud.getPointCount());
  }

  // geometry encoding
  if (1) {
    PayloadBuffer& payload = payload_geom;

    auto& clock_user = clock_user_geom;
    clock_user.start();

    encodeGeometryBrick(params, &payload, attrInterPredParams);

    clock_user.stop();

    auto motionBits = getMotionBits();
    std::cout << "positions motion field size (approximated): "
              << motionBits / 65536.0 << " bits" << std::endl;

    double bpp = double(8 * payload.size()) / inputPointCloud.getPointCount();
    std::cout << "positions bitstream size " << payload.size() << " B (" << bpp
              << " bpp)\n";

    auto total_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << "positions processing time (user): "
              << total_user.count() / 1000.0 << " s" << std::endl;

    callback->onOutputBuffer(payload);
  }

  // verify that the per-level slice constraint has been met
  // todo(df): avoid hard coded value here (should be level dependent)
  if (params->enforceLevelLimits)
    if (pointCloud.getPointCount() > 5000000)
      throw std::runtime_error(
        std::string("level slice point count limit (5000000) exceeded: ")
        + std::to_string(pointCloud.getPointCount()));

  // non local attributes encoding
  if (!params->localized_attributes_encoding) {
    // global recolouring
    // NB: recolouring is required if points are added / removed
    if (!_gps->geom_duplicated_points_flag || _gps->trisoup_enabled_flag) {
      for (const auto& attr_sps : _sps->attributeSets) {
        recolour(
          attr_sps, params->recolour, originPartCloud, _srcToCodingScale,
          _originInCodingCoords + _sliceOrigin, &pointCloud);
      }
    }

    // dump recoloured point cloud
    // todo(df): this needs to work with partitioned clouds
    callback->onPostRecolour(pointCloud);

    // attributeCoding

    // for localized attributes
    uint32_t startIdx = 0;
    PCCPointSet3 slabBlockPointCloud;
    point_t slabBlockSize = {
      _sps->localized_attributes_slab_thickness_minus1 + 1,
      _sps->localized_attributes_slab_block_size_minus1 + 1,
      _sps->localized_attributes_slab_block_size_minus1 + 1
    };
    attrInterPredParams.setSlabBlockSize(slabBlockSize);

    for (const auto& startAndNumPts : startAndNumPointsPerSlabBlock) {
      point_t slabBlockStart = startAndNumPts.first;
      int numPts = startAndNumPts.second;
      slabBlockPointCloud.clear();
      slabBlockPointCloud.appendPartition(pointCloud, startIdx, startIdx+numPts);
      attrInterPredParams.setSlabBlockStart(slabBlockStart);
      // for each attribute
      for (const auto& it : params->attributeIdxMap) {
        int attrIdx = it.second;
        const auto& attr_sps = _sps->attributeSets[attrIdx];
        const auto& attr_aps = *_aps[attrIdx];
        const auto& attr_enc = params->attr[attrIdx];
        const auto motionEncoder = attrEncoder[attrIdx]->getMotionEncoder();

        auto& clock_user = clock_user_attr[attrIdx];
        clock_user.start();

        // local motion search performed by slab
        if (attrInterPredParams.enableAttrInterPred) {
          attrInterPredParams.prepareEncodeMotion(
            attr_aps.dual_motion_field_flag ? attr_aps.motion : _gps->motion,
            *_gps, _gbh, slabBlockPointCloud, slabBlockStart);
        }
        if (attrInterPredParams.enableAttrInterPred && attr_aps.dual_motion_field_flag) {
          attrInterPredParams.findMotion(
            params, attr_enc.motion, attr_aps.motion, *_gps, _gbh, *motionEncoder, slabBlockPointCloud);
        } else if (attrInterPredParams.enableAttrInterPred && _gbh.slice_inter_prediction_flag) {
          if(startAndNumPointsPerSlabBlock.size() > 1)
            attrInterPredParams.extractMotionForSlabBlock(slabBlockStart, slabBlockSize);
          else
            attrInterPredParams.copyMotion();
        }
        // local attributes coding
        attrEncoder[attrIdx]->encodeSlabBlock(
          *_sps, *_gps, attr_sps, attr_aps, _gbh, slabBlockPointCloud, nullptr,
          attr_enc, attrInterPredParams);

        clock_user.stop();
      }
      // set reconstructed colors values
      pointCloud.setFromPartition(slabBlockPointCloud, 0, numPts, startIdx);

      startIdx += numPts;
    }
  }
  for (const auto& it : params->attributeIdxMap) {
    int attrIdx = it.second;
    const auto& attr_sps = _sps->attributeSets[attrIdx];
    const auto& attr_aps = *_aps[attrIdx];
    const auto& label = attr_sps.attributeLabel;

    PayloadBuffer& payload = payload_attr[attrIdx];

    auto& clock_user = clock_user_attr[attrIdx];

    auto& ctxtMemAttr = _ctxtMemAttrs.at(_abh[attrIdx].attr_sps_attr_idx);
    auto& ctxtMemPredMode = _ctxtMemPredMode.at(_abh[attrIdx].attr_sps_attr_idx);
    auto& ctxtMemDualMotion = _ctxtMemDualMotion.at(_abh[attrIdx].attr_sps_attr_idx);
    attrEncoder[attrIdx]->finishEncode(
      *_sps, *_gps, attr_sps, attr_aps, ctxtMemAttr, ctxtMemPredMode, ctxtMemDualMotion, &payload);

    auto motionBits = attrEncoder[attrIdx]->getMotionBits();
    std::cout << label << "s motion field size (approximated): "
              << motionBits / 65536.0 << " bits" << std::endl;

    int coded_size = int(payload.size());
    double bpp = double(8 * coded_size) / inputPointCloud.getPointCount();
    std::cout << label << "s bitstream size " << coded_size << " B (" << bpp
              << " bpp)\n";

    auto time_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << label
              << "s processing time (user): " << time_user.count() / 1000.0
              << " s" << std::endl;

    callback->onOutputBuffer(payload);
  }

  // Note the current slice id for loss detection with entropy continuation
  _prevSliceId = _sliceId;

  // prevent re-use of this sliceId:  the next slice (geometry + attributes)
  // should be distinguishable from the current slice.
  _sliceId++;
  _firstSliceInFrame = false;

  if (reconCloud)
    appendSlice(reconCloud->cloud);
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::processNextSlabAttributes(
  PCCPointSet3& slabPointCloud,
  const point_t& slabStart,
  std::vector<int32_t>& numPointsSlabBlocks,
  bool isLast
)
{
  if (!slabPointCloud.getPointCount())
    return;

  clock_user_geom.stop();
  ++currSlabIdx;
  //const point_t slabBlockSize = _sps->localized_attributes_slab_thickness_minus1 + 1;
  if (!currSlabIdx)
    attrInterPredParams.setSlabBlockSize(slabBlockSize);

  const point_t numSlabBlocksPerDim = {
    (_gbh.maxRootNodeSize + slabBlockSize[0] - 1) / slabBlockSize[0],
    (_gbh.maxRootNodeSize + slabBlockSize[1] - 1) / slabBlockSize[1],
    (_gbh.maxRootNodeSize + slabBlockSize[2] - 1) / slabBlockSize[2]
  };

  if (type == Type::kLocalEncoder) { // encoder
    // TODO: origin cloud could be ordered in lexicographic order so that
    // extraction for a slab becomes quite fast.
    PCCPointSet3 localOriginPartCloud;
    std::vector<int> index_local_part;
    if ((!_gps->geom_duplicated_points_flag || _gps->trisoup_enabled_flag) && originPartCloud->size()) {
      localOriginPartCloud.reserve(originPartCloud->size());
      index_local_part.reserve(originPartCloud->size());
      // extract slab from original point cloud for recoloring
      const int _margin = 0;
      const int xStart = slabStart[0];
      const int xEnd = xStart + slabBlockSize[0];
      int xStartOrigin = currSlabIdx == 0 ? bBoxOrigin.min[0] :
        std::round((xStart + origin[0] - _margin) * targetToSourceScaleFactor);
      int xEndOrigin = isLast ? bBoxOrigin.max[0] + 1 :
        std::round((xEnd + origin[0] + _margin) * targetToSourceScaleFactor);
      do {
        localOriginPartCloud.clear();
        index_local_part.clear();

        for (int i = 0; i < originPartCloud->size(); ++i) {
          const auto& pt = (*originPartCloud)[i];
          if(pt[0] >= xStartOrigin && pt[0] < xEndOrigin) {
            index_local_part.push_back(i);
          }
        }
        --xStartOrigin;
        ++xEndOrigin;
      }
      while(index_local_part.size() < 1);

      localOriginPartCloud.appendPartition(*originPartCloud, index_local_part, true);
    }

    for (const auto& attributeIdxMap : params->attributeIdxMap) {
      int attrIdx = attributeIdxMap.second;
      const auto& attr_sps = _sps->attributeSets[attrIdx];
      const auto& attr_aps = *_aps[attrIdx];
      const auto& attr_enc = params->attr[attrIdx];
      const auto& label = attr_sps.attributeLabel;
      const auto motionEncoder = attrEncoder[attrIdx]->getMotionEncoder();


      if (! (label == KnownAttributeLabel::kColour
          || label == KnownAttributeLabel::kReflectance))
        continue;

      // local recolor
      // NB: recolouring is required if points are added / removed
      if (!_gps->geom_duplicated_points_flag || _gps->trisoup_enabled_flag) {
        recolour(
          attr_sps, params->recolour, localOriginPartCloud, _srcToCodingScale,
          origin, &slabPointCloud);
      }

      // TODO: store local recolored point cloud for later
      //    dump of recoloured point cloud
      //    using: callback->onPostRecolour(recoloredPointCloud);

      // recoloring was not part of the color processing time but motion search
      clock_user_attr[attrIdx].start();

      PCCPointSet3 slabBlockPointCloud;
      slabBlockPointCloud.addRemoveAttributes(slabPointCloud);
      slabBlockPointCloud.reserve(slabPointCloud.size());

      int slabBlockIdx = 0;
      uint32_t pointStartIndex = 0;
      for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
        for (int z = 0; z < numSlabBlocksPerDim[2]; ++z, ++slabBlockIdx)
          if (numPointsSlabBlocks[slabBlockIdx]) {
            slabBlockPointCloud.clear();
            slabBlockPointCloud.appendPartition(
              slabPointCloud, pointStartIndex, pointStartIndex + numPointsSlabBlocks[slabBlockIdx]);

            point_t slabBlockStart = {
              slabStart[0], y * slabBlockSize[1], z * slabBlockSize[2] };
            attrInterPredParams.setSlabBlockStart(slabBlockStart);
            // local motion search performed by slab
            if (attrInterPredParams.enableAttrInterPred) {
              attrInterPredParams.prepareEncodeMotion(
                attr_aps.dual_motion_field_flag ? attr_aps.motion : _gps->motion,
                *_gps, _gbh, slabBlockPointCloud, slabBlockStart);
            }
            if (attrInterPredParams.enableAttrInterPred && attr_aps.dual_motion_field_flag) {
              attrInterPredParams.findMotion(
                params, attr_enc.motion, attr_aps.motion, *_gps, _gbh, *motionEncoder, slabBlockPointCloud);
            } else if (attrInterPredParams.enableAttrInterPred && _gbh.slice_inter_prediction_flag) {
              if (!_gbh.uniqueSlabBlock)
                attrInterPredParams.extractMotionForSlabBlock(slabBlockStart, slabBlockSize);
              else
                attrInterPredParams.copyMotion();
            }
            // local attributes coding
            attrEncoder[attrIdx]->encodeSlabBlock(
              *_sps, *_gps, attr_sps, attr_aps, _gbh, slabBlockPointCloud, nullptr,
              attr_enc, attrInterPredParams);
            // retrieve coded result (assuming points order have changed)
            slabPointCloud.setFromPartition(
              slabBlockPointCloud, 0, numPointsSlabBlocks[slabBlockIdx], pointStartIndex);
            pointStartIndex += numPointsSlabBlocks[slabBlockIdx];
          }

      clock_user_attr[attrIdx].stop();
    }
  } else if (type == Type::kGlobalEncoder) { // encoder
    // assume point order will not change before attributes processing
    // register slab for size for later processing
    int slabBlockIdx = 0;
    for (int y = 0; y < numSlabBlocksPerDim[1]; ++y)
      for (int z = 0; z < numSlabBlocksPerDim[2]; ++z, ++slabBlockIdx)
        if (numPointsSlabBlocks[slabBlockIdx]) {
          point_t slabBlockStart = {
            slabStart[0], y * slabBlockSize[1], z * slabBlockSize[2] };

          startAndNumPointsPerSlabBlock.push_back(std::make_pair(slabBlockStart, numPointsSlabBlocks[slabBlockIdx]));
        }
  }
  clock_user_geom.start();
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::startEncodeGeometryBrick(
  const EncoderParams* params)
{
  GeometryBrickHeader& gbh = _gbh; /* put directly into _gbh which is also
    used by local attributes */
  gbh = GeometryBrickHeader(); // reset content

  gbh.geom_geom_parameter_set_id = _gps->gps_geom_parameter_set_id;
  gbh.geom_slice_id = _sliceId;
  gbh.prev_slice_id = _prevSliceId;
  // NB: slice_tag could be set to some other (external) meaningful value
  gbh.slice_tag = std::max(0, _tileId);
  gbh.frame_ctr_lsb = _frameCounter & ((1 << _sps->frame_ctr_bits) - 1);
  gbh.geomBoxOrigin = _sliceOrigin;
  gbh.geom_box_origin_bits_minus1 = numBits(gbh.geomBoxOrigin.max()) - 1;
  gbh.geom_box_log2_scale = 0;
  gbh.geom_stream_cnt_minus1 = params->gbh.geom_stream_cnt_minus1;
  gbh.trisoup_node_size = params->gbh.trisoup_node_size;
  gbh.trisoup_QP = params->gbh.trisoup_QP;
  gbh.qu_size_log2 = params->gbh.qu_size_log2;
  gbh.trisoup_centroid_vertex_residual_flag =
    params->gbh.trisoup_centroid_vertex_residual_flag;
  gbh.trisoup_face_vertex_flag =
    params->gbh.trisoup_face_vertex_flag;
  gbh.trisoup_halo_flag =
    params->gbh.trisoup_halo_flag;
  gbh.trisoup_vertex_merge_flag = params->gbh.trisoup_vertex_merge_flag;
  gbh.trisoup_vertex_consistency_flag = params->gbh.trisoup_vertex_consistency_flag;
  gbh.trisoup_thickness = params->gbh.trisoup_thickness;
  gbh.slice_inter_prediction_flag = _codeCurrFrameAsInter;
  gbh.uniqueSlabBlock = params->singleSlabOptimizedImplementationEnabled;

  // Entropy continuation is not permitted in the first slice of a frame
  gbh.entropy_continuation_flag = false;
  if (_sps->entropy_continuation_enabled_flag)
    gbh.entropy_continuation_flag = !_firstSliceInFrame;

  int nodeSizeCommon, nodeSizeLog2Common;
  if (params->trisoup.alignToNodeGrid && !params->trisoupNodeSizes.empty()) {
    // TODO: do we really want/need to align upper boundaries with trisoup grid ?
    nodeSizeCommon = lcm_all(
      params->trisoupNodeSizes.begin(),
      params->trisoupNodeSizes.end());

    nodeSizeLog2Common = ilog2(uint32_t(nodeSizeCommon - 1)) + 1;
  } else {
    nodeSizeCommon = gbh.trisoupNodeSize(*_gps);
    nodeSizeLog2Common = gbh.trisoupNodeSizeLog2(*_gps);
  }

  // inform the geometry coder what the root node size is
  // NB: A minimum whd of 2 means there is always at least 1 tree level
  gbh.rootNodeSizeLog2 = ceillog2(std::max(2,
    ((_sliceBoxWhd.max() + nodeSizeCommon - 1) / nodeSizeCommon)))
    + nodeSizeLog2Common; // dicing if any

  // The root node size cannot be smaller than the trisoup node size
  // since this is how the root node size is defined at the decoder.
  gbh.rootNodeSizeLog2 = std::max(gbh.trisoupNodeSizeLog2(*_gps), gbh.rootNodeSizeLog2);

  // forget (reset) all saved context state at boundary
  if (!gbh.entropy_continuation_flag) {
    if (
      !_gps->gof_geom_entropy_continuation_enabled_flag
      || !gbh.slice_inter_prediction_flag) {
      _ctxtMemOctreeGeom->reset();
      _ctxtMemMotion.reset();
    }
    if (!gbh.slice_inter_prediction_flag) {
      for (auto& ctxtMem : _ctxtMemAttrs)
        ctxtMem.reset();
      for (auto& ctxtMem : _ctxtMemPredMode)
        ctxtMem.reset();
      for (auto& ctxtMem : _ctxtMemDualMotion)
        ctxtMem.reset();
    }
  }
}

void
PCCTMC3Encoder3::encodeGeometryBrick(
  const EncoderParams* params,
  PayloadBuffer* buf,
  AttributeInterPredParams& attrInterPredParams)
{
  GeometryBrickHeader& gbh = _gbh;

  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = int(pointCloud.getPointCount()) * 3 * 4 + 1024;

  // allocate entropy streams
  std::vector<std::unique_ptr<EntropyEncoder>> arithmeticEncoders;
  for (int i = 0; i < 1 + gbh.geom_stream_cnt_minus1; i++) {
    arithmeticEncoders.emplace_back(new EntropyEncoder(maxAcBufLen, nullptr));
    auto& aec = arithmeticEncoders.back();
    aec->setBypassBinCodingWithoutProbUpdate(_sps->bypass_bin_coding_without_prob_update);
    aec->start();
  }

  if (!_gps->trisoup_enabled_flag) {
    encodeGeometryOctree(
      *params, *_gps, gbh, pointCloud, *_ctxtMemOctreeGeom, _ctxtMemMotion,
      arithmeticEncoders, _refFrame, *_sps, attrInterPredParams, *this);
  }
  else
  {
    // limit the number of points to the slice limit
    // todo(df): this should be derived from the level
    gbh.footer.geom_num_points_minus1 = params->partition.sliceMaxPointsTrisoup - 1;
    encodeGeometryTrisoup(
      *params, *_gps, gbh, pointCloud,
      *_ctxtMemOctreeGeom, _ctxtMemMotion, arithmeticEncoders, _refFrame,
      *_sps, attrInterPredParams, *this);
  }

  // signal the actual number of points coded
  gbh.footer.geom_num_points_minus1 = pointCloud.getPointCount() - 1;

  // assemble data unit
  //  - record the position of each aec buffer for chunk concatenation
  std::vector<std::pair<size_t, size_t>> aecStreams;
  write(*_sps, *_gps, gbh, buf);
  for (auto& arithmeticEncoder : arithmeticEncoders) {
    auto aecLen = arithmeticEncoder->stop();
    auto aecBuf = arithmeticEncoder->buffer();
    aecStreams.emplace_back(buf->size(), aecLen);
    buf->insert(buf->end(), aecBuf, aecBuf + aecLen);
  }

  // append the footer
  write(*_gps, gbh, gbh.footer, buf);
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::appendSlice(PCCPointSet3& accumCloud)
{
  // offset current point cloud to be in coding coordinate system
  size_t numPoints = pointCloud.getPointCount();
  for (size_t i = 0; i < numPoints; i++)
    for (int k = 0; k < 3; k++)
      pointCloud[i][k] += _sliceOrigin[k];

  accumCloud.append(pointCloud);
}

//----------------------------------------------------------------------------
// translates and scales inputPointCloud, storing the result in
// this->pointCloud for use by the encoding process.

SrcMappedPointSet
PCCTMC3Encoder3::quantization(const PCCPointSet3& src)
{
  // Currently the sequence bounding box size must be set
  assert(_sps->seqBoundingBoxSize != Vec3<int>{0});

  // Clamp all points to [clampBox.min, clampBox.max] after translation
  // and quantisation.
  Box3<int32_t> clampBox(0, std::numeric_limits<int32_t>::max());

  // When using predictive geometry, sub-sample the point cloud and let
  // the predictive geometry coder quantise internally.
  if (_inputDecimationScale != 1.)
    return samplePositionsUniq(
      _inputDecimationScale, _srcToCodingScale, _originInCodingCoords, src);

  if (!_gps->geom_duplicated_points_flag)
    return quantizePositionsUniq(
      _srcToCodingScale, _originInCodingCoords, clampBox, src);

  SrcMappedPointSet dst;
  quantizePositions(
    _srcToCodingScale, _originInCodingCoords, clampBox, src, &dst.cloud);
  return dst;
}

//----------------------------------------------------------------------------
// get the partial point cloud according to required point indexes

PCCPointSet3
getPartition(const PCCPointSet3& src, const std::vector<int32_t>& indexes)
{
  PCCPointSet3 dst;
  dst.addRemoveAttributes(src);

  int partitionSize = indexes.size();
  dst.resize(partitionSize);

  for (int i = 0; i < partitionSize; i++) {
    int inputIdx = indexes[i];
    dst[i] = src[inputIdx];

    if (src.hasColors())
      dst.setColor(i, src.getColor(inputIdx));

    if (src.hasReflectances())
      dst.setReflectance(i, src.getReflectance(inputIdx));

    if (src.hasLaserAngles())
      dst.setLaserAngle(i, src.getLaserAngle(inputIdx));
  }

  return dst;
}

//----------------------------------------------------------------------------
// get the partial point cloud according to required point indexes

PCCPointSet3
getPartition(
  const PCCPointSet3& src,
  const SrcMappedPointSet& map,
  const std::vector<int32_t>& indexes)
{
  // Without the list, do nothing
  if (map.idxToSrcIdx.empty())
    return {};

  // work out the destination size.
  // loop over each linked list until an element points to itself
  int size = 0;
  for (int idx : indexes) {
    int prevIdx, srcIdx = map.idxToSrcIdx[idx];
    do {
      size++;
      prevIdx = srcIdx;
      srcIdx = map.srcIdxDupList[srcIdx];
    } while (srcIdx != prevIdx);
  }

  PCCPointSet3 dst;
  dst.addRemoveAttributes(src);
  dst.resize(size);

  int dstIdx = 0;
  for (int idx : indexes) {
    int prevIdx, srcIdx = map.idxToSrcIdx[idx];
    do {
      dst[dstIdx] = src[srcIdx];

      if (src.hasColors())
        dst.setColor(dstIdx, src.getColor(srcIdx));

      if (src.hasReflectances())
        dst.setReflectance(dstIdx, src.getReflectance(srcIdx));

      if (src.hasLaserAngles())
        dst.setLaserAngle(dstIdx, src.getLaserAngle(srcIdx));

      dstIdx++;
      prevIdx = srcIdx;
      srcIdx = map.srcIdxDupList[srcIdx];
    } while (srcIdx != prevIdx);
  }

  return dst;
}

//============================================================================

}  // namespace pcc
