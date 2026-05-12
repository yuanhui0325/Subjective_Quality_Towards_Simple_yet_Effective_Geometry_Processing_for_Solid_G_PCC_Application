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

#include "TMC3.h"

#include <memory>

#include "PCCTMC3Encoder.h"
#include "PCCTMC3Decoder.h"
#include "constants.h"
#include "ply.h"
#include "pointset_processing.h"
#include "program_options_lite.h"
#include "io_tlv.h"
#include "version.h"

using namespace std;
using namespace pcc;

//============================================================================

enum class OutputSystem
{
  // Output after global scaling, don't convert to external system
  kConformance = 0,

  // Scale output to external coordinate system
  kExternal = 1,
};

//----------------------------------------------------------------------------

struct Parameters {
  bool isDecoder;

  // Scale factor to apply when loading the ply before integer conversion.
  // Eg, If source point positions are in fractional metres converting to
  // millimetres will allow some fidelity to be preserved.
  double inputScale;

  // Length of the output point clouds unit vectors.
  double outputUnitLength;

  // output mode for ply writing (binary or ascii)
  bool outputBinaryPly;

  // Fractional fixed-point bits retained in conformance output
  int outputFpBits;

  // Output coordinate system to use
  OutputSystem outputSystem;

  // when true, configure the encoder as if no attributes are specified
  bool disableAttributeCoding;

  // Frame number of first file in input sequence.
  int firstFrameNum;

  // Number of frames to process.
  int frameCount;

  // use optimized implementation when full brick is encompassed by a single
  // slab block
  bool singleSlabOptimizedImplementationEnabled;

  std::string uncompressedDataPath;
  std::string compressedStreamPath;
  std::string reconstructedDataPath;

  // Filename for saving recoloured point cloud (encoder).
  std::string postRecolorPath;

  // Filename for saving pre inverse scaled point cloud (decoder).
  std::string preInvScalePath;

  pcc::EncoderParams encoder;
  pcc::DecoderParams decoder;

  // perform attribute colourspace conversion on ply input/output.
  bool convertColourspace;
};

//----------------------------------------------------------------------------

class SequenceCodec {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceCodec(Parameters* params) : params(params) {}

  // Perform conversions and write output point cloud
  //  \params cloud  a mutable copy of reconFrame.cloud
  void writeOutputFrame(
    const std::string& postInvScalePath,
    const std::string& preInvScalePath,
    const CloudFrame& reconFrame,
    PCCPointSet3& cloud);

  // determine the output ply scale factor
  double outputScale(const CloudFrame& cloud) const;

  // the output ply origin, scaled according to output coordinate system
  Vec3<double> outputOrigin(const CloudFrame& cloud) const;

  void scaleAttributesForInput(
    const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

  void scaleAttributesForOutput(
    const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

protected:
  Parameters* params;
};

//----------------------------------------------------------------------------

class SequenceEncoder
  : public SequenceCodec
  , PCCTMC3Encoder3::Callbacks {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceEncoder(Parameters* params);

  int compress(Stopwatch* clock);

protected:
  int compressOneFrame(Stopwatch* clock);

  void onOutputBuffer(const PayloadBuffer& buf) override;
  void onPostRecolour(const PCCPointSet3& cloud) override;

private:
  ply::PropertyNameMap _plyAttrNames;

  PCCTMC3Encoder3 encoder;

  std::ofstream bytestreamFile;

  int frameNum;
};

//----------------------------------------------------------------------------

class SequenceDecoder
  : public SequenceCodec
  , PCCTMC3Decoder3::Callbacks {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceDecoder(Parameters* params);

  int decompress(Stopwatch* clock);

protected:
  void onOutputCloud(const CloudFrame& cloud) override;

  int getOuputFrameNumber(const CloudFrame& frame) override;

private:
  PCCTMC3Decoder3 decoder;

  std::ofstream bytestreamFile;

  Stopwatch* clock;
};

//============================================================================

void convertToGbr(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

void convertFromGbr(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud);

//============================================================================

int
main(int argc, char* argv[])
{
  cout << "MPEG PCC tmc3 version " << ::pcc::version << endl;

  Parameters params;

  try {
    if (!ParseParameters(argc, argv, params))
      return 1;
  }
  catch (df::program_options_lite::ParseFailure& e) {
    std::cerr << "Error parsing option \"" << e.arg << "\" with argument \""
              << e.val << "\"." << std::endl;
    return 1;
  }

  // Timers to count elapsed wall/user time
  pcc::chrono::Stopwatch<std::chrono::steady_clock> clock_wall;
  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;

  clock_wall.start();

  int ret = 0;
  if (params.isDecoder) {
    ret = SequenceDecoder(&params).decompress(&clock_user);
  } else {
    ret = SequenceEncoder(&params).compress(&clock_user);
  }

  clock_wall.stop();

  using namespace std::chrono;
  auto total_wall = duration_cast<milliseconds>(clock_wall.count()).count();
  auto total_user = duration_cast<milliseconds>(clock_user.count()).count();
  std::cout << "Processing time (wall): " << total_wall / 1000.0 << " s\n";
  std::cout << "Processing time (user): " << total_user / 1000.0 << " s\n";

  return ret;
}

//---------------------------------------------------------------------------

std::array<const char*, 3>
axisOrderToPropertyNames(AxisOrder order)
{
  static const std::array<const char*, 3> kAxisOrderToPropertyNames[] = {
    {"z", "y", "x"}, {"x", "y", "z"}, {"x", "z", "y"}, {"y", "z", "x"},
    {"z", "y", "x"}, {"z", "x", "y"}, {"y", "x", "z"}, {"x", "y", "z"},
  };

  return kAxisOrderToPropertyNames[int(order)];
}

//---------------------------------------------------------------------------
// :: Command line / config parsing helpers

template<typename T>
static std::istream&
readUInt(std::istream& in, T& val)
{
  unsigned int tmp;
  in >> tmp;
  val = T(tmp);
  return in;
}

namespace pcc {
static std::istream&
operator>>(std::istream& in, ScaleUnit& val)
{
  try {
    readUInt(in, val);
  }
  catch (...) {
    in.clear();
    std::string str;
    in >> str;

    val = ScaleUnit::kDimensionless;
    if (str == "metre")
      val = ScaleUnit::kMetre;
    else if (!str.empty())
      throw std::runtime_error("Cannot parse unit");
  }
  return in;
}
}  // namespace pcc

static std::istream&
operator>>(std::istream& in, OutputSystem& val)
{
  return readUInt(in, val);
}

namespace pcc {
static std::istream&
operator>>(std::istream& in, ColourMatrix& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, AxisOrder& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, PartitionMethod& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, MotionPreset& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, TriStateEncParam& val)
{
  return readUInt(in, val);
}
}  // namespace pcc


static std::ostream&
operator<<(std::ostream& out, const OutputSystem& val)
{
  switch (val) {
  case OutputSystem::kConformance: out << "0 (Conformance)"; break;
  case OutputSystem::kExternal: out << "1 (External)"; break;
  }
  return out;
}

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const ScaleUnit& val)
{
  switch (val) {
  case ScaleUnit::kDimensionless: out << "0 (Dimensionless)"; break;
  case ScaleUnit::kMetre: out << "1 (Metre)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const ColourMatrix& val)
{
  switch (val) {
  case ColourMatrix::kIdentity: out << "0 (Identity)"; break;
  case ColourMatrix::kBt709: out << "1 (Bt709)"; break;
  case ColourMatrix::kUnspecified: out << "2 (Unspecified)"; break;
  case ColourMatrix::kReserved_3: out << "3 (Reserved)"; break;
  case ColourMatrix::kUsa47Cfr73dot682a20:
    out << "4 (Usa47Cfr73dot682a20)";
    break;
  case ColourMatrix::kBt601: out << "5 (Bt601)"; break;
  case ColourMatrix::kSmpte170M: out << "6 (Smpte170M)"; break;
  case ColourMatrix::kSmpte240M: out << "7 (Smpte240M)"; break;
  case ColourMatrix::kYCgCo: out << "8 (kYCgCo)"; break;
  case ColourMatrix::kBt2020Ncl: out << "9 (Bt2020Ncl)"; break;
  case ColourMatrix::kBt2020Cl: out << "10 (Bt2020Cl)"; break;
  case ColourMatrix::kSmpte2085: out << "11 (Smpte2085)"; break;
  default: out << "Unknown"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const AxisOrder& val)
{
  switch (val) {
  case AxisOrder::kZYX: out << "0 (zyx)"; break;
  case AxisOrder::kXYZ: out << "1 (xyz)"; break;
  case AxisOrder::kXZY: out << "2 (xzy)"; break;
  case AxisOrder::kYZX: out << "3 (yzx)"; break;
  case AxisOrder::kZYX_4: out << "4 (zyx)"; break;
  case AxisOrder::kZXY: out << "5 (zxy)"; break;
  case AxisOrder::kYXZ: out << "6 (yxz)"; break;
  case AxisOrder::kXYZ_7: out << "7 (xyz)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const PartitionMethod& val)
{
  switch (val) {
  case PartitionMethod::kNone: out << "0 (None)"; break;
  case PartitionMethod::kUniformGeom: out << "2 (UniformGeom)"; break;
  case PartitionMethod::kOctreeUniform: out << "3 (UniformOctree)"; break;
  case PartitionMethod::kUniformSquare: out << "4 (UniformSquare)"; break;
  case PartitionMethod::kNpoints: out << "5 (NPointSpans)"; break;
  default: out << int(val) << " (Unknown)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const MotionPreset& val)
{
  switch (val) {
  case MotionPreset::kNoMotion: out << "0 (NoMotion)"; break;
  case MotionPreset::kMotionPresetForOctree: out << "2 (MotionPresetForOctree)"; break;
  case MotionPreset::kMotionPresetForTrisoup: out << "3 (MotionPresetForTrisoup)"; break;
  default: out << int(val) << " (Unknown)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const TriStateEncParam& val)
{
  switch (val) {
  case TriStateEncParam::kNotSet: out << "-1 (NotSet)"; break;
  case TriStateEncParam::kFalse: out << "0"; break;
  case TriStateEncParam::kTrue: out << "1"; break;
  default: out << int(val) << " (Unknown)"; break;
  }
  return out;
}
}  // namespace pcc


namespace df {
namespace program_options_lite {
  template<typename T>
  struct option_detail<pcc::Vec3<T>> {
    static constexpr bool is_container = true;
    static constexpr bool is_fixed_size = true;
    typedef T* output_iterator;

    static void clear(pcc::Vec3<T>& container){};
    static output_iterator make_output_iterator(pcc::Vec3<T>& container)
    {
      return &container[0];
    }
  };
}  // namespace program_options_lite
}  // namespace df

//---------------------------------------------------------------------------
// :: Command line / config parsing

void sanitizeEncoderOpts(
  Parameters& params, df::program_options_lite::ErrorReporter& err);

void sanitizeDecoderOpts(
  Parameters& params, df::program_options_lite::ErrorReporter& err);

//---------------------------------------------------------------------------

bool
ParseParameters(int argc, char* argv[], Parameters& params)
{
  namespace po = df::program_options_lite;

  struct {
    AttributeDescription desc;
    AttributeParameterSet aps;
    EncoderAttributeParams encoder;
  } params_attr;

  bool print_help = false;

  // a helper to set the attribute
  std::function<po::OptionFunc::Func> attribute_setter =
    [&](po::Options&, const std::string& name, po::ErrorReporter) {
      // copy the current state of parsed attribute parameters
      //
      // NB: this does not cause the default values of attr to be restored
      // for the next attribute block.  A side-effect of this is that the
      // following is allowed leading to attribute foo having both X=1 and
      // Y=2:
      //   "--attr.X=1 --attribute foo --attr.Y=2 --attribute foo"
      //

      // NB: insert returns any existing element
      const auto& it = params.encoder.attributeIdxMap.insert(
        {name, int(params.encoder.attributeIdxMap.size())});

      if (it.second) {
        params.encoder.sps.attributeSets.push_back(params_attr.desc);
        params.encoder.aps.push_back(params_attr.aps);
        params.encoder.attr.push_back(params_attr.encoder);
        return;
      }

      // update existing entry
      params.encoder.sps.attributeSets[it.first->second] = params_attr.desc;
      params.encoder.aps[it.first->second] = params_attr.aps;
      params.encoder.attr[it.first->second] = params_attr.encoder;
    };

  /* clang-format off */
  // The definition of the program/config options, along with default values.
  //
  // NB: when updating the following tables:
  //      (a) please keep to 80-columns for easier reading at a glance,
  //      (b) do not vertically align values -- it breaks quickly
  //
  po::Options opts;
  opts.addOptions()
  ("help", print_help, false, "this help text")
  ("config,c", po::parseConfigFile, "configuration file name")

  (po::Section("General"))

  ("mode", params.isDecoder, false,
    "The encoding/decoding mode:\n"
    "  0: encode\n"
    "  1: decode")

  // i/o parameters
  ("firstFrameNum",
     params.firstFrameNum, 0,
     "Frame number for use with interpolating %d format specifiers "
     "in input/output filenames")

  ("frameCount",
     params.frameCount, 1,
     "Number of frames to encode")

  ("reconstructedDataPath",
    params.reconstructedDataPath, {},
    "The ouput reconstructed pointcloud file path (decoder only)")

  ("uncompressedDataPath",
    params.uncompressedDataPath, {},
    "The input pointcloud file path")

  ("compressedStreamPath",
    params.compressedStreamPath, {},
    "The compressed bitstream path (encoder=output, decoder=input)")

  ("postRecolorPath",
    params.postRecolorPath, {},
    "Recolored pointcloud file path (encoder only)")

  ("preInvScalePath",
    params.preInvScalePath, {},
    "Pre inverse scaled pointcloud file path (decoder only)")

  ("convertPlyColourspace",
    params.convertColourspace, true,
    "Convert ply colourspace according to attribute colourMatrix")

  ("outputBinaryPly",
    params.outputBinaryPly, true,
    "Output ply files using binary (or ascii) format")

  ("outputUnitLength",
    params.outputUnitLength, 0.,
    "Length of reconstructed point cloud x,y,z unit vectors\n"
    " 0: use srcUnitLength")

  ("outputScaling",
    params.outputSystem, OutputSystem::kExternal,
    "Output coordnate system scaling\n"
    " 0: Conformance\n"
    " 1: External")

  ("outputPrecisionBits",
    params.outputFpBits, -1,
    "Fractional bits in conformance output (prior to external scaling)\n"
    " 0: integer,  -1: automatic (full)")

  ("singleSlabOptimizedImplementationEnabled",
    params.singleSlabOptimizedImplementationEnabled, true,
    "Use optimized implementation for slab block encompassing "
    "the full frame (encoder/decoder, non normative)."
  )

  // This section controls all general geometry scaling parameters
  (po::Section("Coordinate system scaling"))

  ("srcUnitLength",
    params.encoder.srcUnitLength, 1.,
    "Length of source point cloud x,y,z unit vectors in srcUnits")

  ("srcUnit",
    params.encoder.sps.seq_geom_scale_unit_flag, ScaleUnit::kDimensionless,
    " 0: dimensionless\n 1: metres")

  ("inputScale",
    params.inputScale, 1.,
    "Scale input while reading src ply. "
    "Eg, 1000 converts metres to integer millimetres")

  ("codingScale",
    params.encoder.codedGeomScale, 1.,
    "Scale used to represent coded geometry. Relative to inputScale")

  ("sequenceScale",
    params.encoder.seqGeomScale, 1.,
    "Scale used to obtain sequence coordinate system. "
    "Relative to inputScale")

  ("externalScale",
    params.encoder.extGeomScale, 1.,
    "Scale used to define external coordinate system.\n"
    "Meaningless when srcUnit = metres\n"
    "  0: Use srcUnitLength\n"
    " >0: Relative to inputScale")

  (po::Section("Decoder"))

  ("skipOctreeLayers",
    params.decoder.minGeomNodeSizeLog2, 0,
    "Partial decoding of octree and attributes\n"
    " 0   : Full decode\n"
    " N>0 : Skip the bottom N layers in decoding process")

  ("decodeMaxPoints",
    params.decoder.decodeMaxPoints, 0,
    "Partially decode up to N points")

  (po::Section("Encoder"))

  ("geometry_axis_order",
    params.encoder.sps.geometry_axis_order, AxisOrder::kXYZ,
    "Sets the geometry axis coding order:\n"
    "  0: (zyx)\n  1: (xyz)\n  2: (xzy)\n"
    "  3: (yzx)\n  4: (zyx)\n  5: (zxy)\n"
    "  6: (yxz)\n  7: (xyz)")

  ("autoSeqBbox",
    params.encoder.autoSeqBbox, true,
    "Calculate seqOrigin and seqSizeWhd automatically.")

  // NB: the underlying variable is in STV order.
  //     Conversion happens during argument sanitization.
  ("seqOrigin",
    params.encoder.sps.seqBoundingBoxOrigin, {0},
    "Origin (x,y,z) of the sequence bounding box "
    "(in input coordinate system). "
    "Requires autoSeqBbox=0")

  // NB: the underlying variable is in STV order.
  //     Conversion happens during argument sanitization.
  ("seqSizeWhd",
    params.encoder.sps.seqBoundingBoxSize, {0},
    "Size of the sequence bounding box "
    "(in input coordinate system). "
    "Requires autoSeqBbox=0")

  ("keepDuplicatedPoints",
    params.encoder.gps.geom_duplicated_points_flag, false,
    "Enables removal of duplicated points")

  ("partitionMethod",
    params.encoder.partition.method, PartitionMethod::kUniformSquare,
    "Method used to partition input point cloud into slices/tiles:\n"
    "  0: none\n"
    "  2: n Uniform-geometry partition bins along the longest edge\n"
    "  3: Uniform geometry partition at n octree depth\n"
    "  4: Uniform square partition\n"
    "  5: n-point spans of input")

  ("safeTrisoupPartionning",
    params.encoder.partition.safeTrisoupPartionning, true,
    "Use safer partitioning to not break Trisoup surfaces\n"
    "  This is compatible with partitionMethod 2 and 4, but sliceMaxPoints\n"
    "  may be exceeded.")

  ("partitionOctreeDepth",
    params.encoder.partition.octreeDepth, 1,
    "Depth of octree partition for partitionMethod=4")

  ("sliceMaxPointsTrisoup",
    params.encoder.partition.sliceMaxPointsTrisoup, 5000000,
    "Maximum number of points per slice")

  ("sliceMaxPoints",
    params.encoder.partition.sliceMaxPoints, 5000000,
    "Maximum number of points per slice")

  ("sliceMinPoints",
    params.encoder.partition.sliceMinPoints, 2500000,
    "Minimum number of points per slice (soft limit)")

  ("tileSize",
    params.encoder.partition.tileSize, 0,
    "Partition input into cubic tiles of given size")

  ("fixedSliceOrigin",
    params.encoder.partition.fixedSliceOrigin, {},
    "Explicitely provided slice origin")

  ("entropyContinuationEnabled",
    params.encoder.sps.entropy_continuation_enabled_flag, false,
    "Propagate context state between slices")

  ("bypassBinCodingWithoutProbUpdate",
    params.encoder.sps.bypass_bin_coding_without_prob_update, true,
    "Codes the bypass bins without using probability update")

  ("GoFGeometryEntropyContinuationEnabled",
    params.encoder.gps.gof_geom_entropy_continuation_enabled_flag, false,
    "Propagate context state between P frames in GoF")

  ("disableAttributeCoding",
    params.disableAttributeCoding, false,
    "Ignore attribute coding configuration")

  ("enforceLevelLimits",
    params.encoder.enforceLevelLimits, true,
    "Abort if level limits exceeded")

  ("fasterMotionSearch",
    params.encoder.motion.approximate_nn, TriStateEncParam::kTrue,
    "Enable approximate nearest neighbor for faster motion search")

   ("motionParamPreset",
     params.encoder.motionPreset, MotionPreset::kNoMotion,
    "Derivation of motion compensation parameters:"
    "  0: All Intra\n"
    "  2: defaults for point clouds compressed with octree\n"
    "  3: defaults for point clouds compressed with trisoup")

  ("motionSearchDistanceStart",
    params.encoder.motion.Amotion0, -1,
    "Starting distance for motion iterative search")

  ("motionSearchLambda",
    params.encoder.motion.lambda, -1.,
    "Lambda used for Lagrangian optimization of the motion field encoding")

  ("motionSearchDColorFactor",
    params.encoder.motion.d_color_factor, -1.,
    "Distortion factor for color in Lagrangian optimization of the motion"
    " field encoding")

  ("motionSearchDGeomFactor",
    params.encoder.motion.d_geom_factor, -1.,
    "Distortion factor for geometry in Lagrangian optimization of the motion"
    " field encoding")

  ("motionSearchDecimate",
    params.encoder.motion.decimate, -1,
    "Decimation factor for accelerating the search for motion")

  ("localizedAttributesEncoding",
    params.encoder.localized_attributes_encoding, false,
    "Enforce the possibility to use localized attributes decoding")

  ("localizedAttributesSlabThickness",
    params.encoder.localized_attributes_slab_thickness, 64,
    "Thickness of localized attributes' slabs")

  ("localizedAttributesSlabBlockSize",
    params.encoder.localized_attributes_slab_block_size, -1,
    "Size (width and height) of localized attributes' slab blocks:\n"
    " <= 0: same a slab thickness,\n"
    " > 0: requested size")

  ("hierarchicalRAHTQP",
    params.encoder.hierarchicalRAHTQP, { -4, 1, 0 },
    "Hierarchical QP for {I, [P, P, ...]} frames\n"
    " QP(s) for P frames are cyclically repeated until next I.")

  ("trisoupEarlySkipStrength",
    params.encoder.trisoupEarlySkipStrength, 10,
    "Higher values favor more skip. (0= no skip,..., 100 = strong skip)")

  (po::Section("Geometry"))

  ("numOctreeEntropyStreams",
    // NB: this is adjusted by minus 1 after the arguments are parsed
    params.encoder.gbh.geom_stream_cnt_minus1, 1,
    "Number of entropy streams for octree coding")

  ("trisoupNodeSize",
    params.encoder.trisoupNodeSizes, { 0 },
    "Node size for surface triangulation\n"
    " <4: disabled")

  ("trisoupQuantizationQP",
    params.encoder.gbh.trisoup_QP, 0,
    "Trisoup QP for quantization of position of vertices along edges\n"
    " 0: quarter voxel\n"
    " 12: one voxel\n"
    " 18: two voxels\n"
    " max is QP=54")

  ("quNodeSizeLog2",
    params.encoder.gbh.qu_size_log2, 0,
    "Size of local Quality Units for geometry:\n"
    " 0: disabled\n"
    " >0: size = trisoupNodeSize * (1 << quNodeSizeLog2)")

  ("trisoupCentroidResidualEnabled",
    params.encoder.gbh.trisoup_centroid_vertex_residual_flag, true,
    "Trisoup activate residual position value for the centroid vertex")

  ("trisoupFaceVertexEnabled",
    params.encoder.gbh.trisoup_face_vertex_flag, true,
    "Trisoup activate the face vertex")

  ("trisoupHaloEnabled",
    params.encoder.gbh.trisoup_halo_flag, true,
    "Trisoup activate halo around triangles for ray tracing")

  ("trisoupVertexThreshold",
    params.encoder.trisoup.thVertexDetermination, 1,
    "minimum number of points near an edge to determine the presence of\n"
    "a trisoup vertex (encoder)\n")

  ("trisoupVertexThresholdForSkip",
    params.encoder.gbh.trisoup_early_skip_vertex_determination_threshold, 1,
    "minimum number of points near an edge to determine the presence of\n"
    "a trisoup vertex (for early skip)\n")

  ("trisoupVertexMergeEnabled",
    params.encoder.gbh.trisoup_vertex_merge_flag, true,
    "Trisoup activates vertex merge during vertex determination.")

  ("trisoupVertexFix",
    params.encoder.gbh.trisoup_vertex_consistency_flag, true,
    "Trisoup activates vertex consistency realtive to the orignial poitn cloud")

  ("trisoupImprovedEncoderEnabled",
    params.encoder.trisoup.improvedVertexDetermination, true,
    "Trisoup activate improved determination of vertex position (encoder only)")

  ("trisoupAlignToNodeGrid",
    params.encoder.trisoup.alignToNodeGrid, true,
    "Align slices to a grid of trisoup nodes (encoder only)")

  ("trisoupSkipModeEnabled",
    params.encoder.gps.trisoup_skip_mode_enabled_flag, true,
    "Enables skip mode for trisoup")

  ("trisoupEarlySkipCodingModeEnabled",
    params.encoder.gps.trisoup_early_skip_coding_mode_enabled_flag, true,
    "Enable early skip coding mode for trisoup")

  ("trisoupThickness",
    params.encoder.gbh.trisoup_thickness, 36,
    "Thickness of Trisoup triangles")

  ("randomAccessPeriod",
    params.encoder.randomAccessPeriod, 1,
    "Distance (in pictures) between random access points when "
    "encoding a sequence")

  ("interPredictionEnabled",
    params.encoder.gps.inter_prediction_enabled_flag, false,
    "Enable inter prediciton")

  ("motionBlockSize",
    params.encoder.gps.motion.motion_block_size, -1,
    "Largest block size for motion")

  ("motionMinPUSize",
    params.encoder.gps.motion.motion_min_pu_size, -1,
    "Minimum block size for motion")

  ("pointCountMetadata",
    params.encoder.gps.octree_point_count_list_present_flag, false,
    "Add octree layer point count metadata")

  (po::Section("Attributes"))

  // attribute processing
  //   NB: Attribute options are special in the way they are applied (see above)
  ("attribute",
    attribute_setter,
    "Encode the given attribute (NB, must appear after the"
    "following attribute parameters)")

  // NB: the cli option sets +1, the minus1 will be applied later
  ("attrScale",
    params_attr.desc.params.attr_scale_minus1, 1,
    "Scale factor used to interpret coded attribute values")

  ("attrOffset",
    params_attr.desc.params.attr_offset, 0,
    "Offset used to interpret coded attribute values")

  ("bitdepth",
    params_attr.desc.bitdepth, 8,
    "Attribute bitdepth")

  ("internalBitdepth",
    params_attr.desc.internalBitdepth, 16,
    "Attribute internal bitdepth")

  ("defaultValue",
    params_attr.desc.params.attr_default_value, {},
    "Default attribute component value(s) in case of data omission")

  // todo(df): this should be per-attribute
  ("colourMatrix",
    params_attr.desc.params.cicp_matrix_coefficients_idx, ColourMatrix::kBt709,
    "Matrix used in colourspace conversion\n"
    "  0: none (identity)\n"
    "  1: ITU-T BT.709\n"
    "  8: YCgCo")

  ("dualMotionEnabled",
    params_attr.aps.dual_motion_field_flag, false,
    "Controls the use of a dual motion field for the attributes")

  ("motionSearchDistanceStartForDual",
    params_attr.encoder.motion.Amotion0, -1,
    "Starting distance for dual motion iterative search (encoder)")

  ("motionSearchLambdaForDual",
    params_attr.encoder.motion.lambda, -1.,
    "Lambda used for Lagrangian optimization of the dual motion field "
    "encoding (encoder)")

  ("motionSearchDColorFactorForDual",
    params_attr.encoder.motion.d_color_factor, -1.,
    "Distortion factor for color in Lagrangian optimization of the dual motion"
    " field encoding (encoder)")

  ("motionSearchDGeomFactorForDual",
    params_attr.encoder.motion.d_geom_factor, -1.,
    "Distortion factor for geometry in Lagrangian optimization of the dual "
    "motion field encoding (encoder)")

  ("motionSearchDecimateForDual",
    params_attr.encoder.motion.decimate, -1,
    "Decimation factor for accelerating the search for dual motion "
    "(encoder)")

  ("motionBlockSizeForDual",
    params_attr.aps.motion.motion_block_size, -1,
    "Largest block size for dual motion")

  ("motionMinPUSizeForDual",
    params_attr.aps.motion.motion_min_pu_size, -1,
    "Minimum block size for dual motion")

  ("slabBlockSkipEnabled",
    params_attr.aps.slab_block_skip_enabled_flag, true,
    "Controls the use of skip for slab block")

  ("slabBlockSkipLambdaFactor",
    params_attr.encoder.slabBlockSkipRDO_lambdaFactor, 12800.,
    "Scaling factor for lambda computation with accurate RDO "
    "for slab block skip")

  ("slabBlockSkipFastDecision",
    params_attr.encoder.slabBlockSkipFastRDO, false,
    "Controls the use of an accelerated decision for slab block skip "
    "to be used instead of accurate RDO.")

  ("slabBlockSkipFastDecisionMSEFactor",
    params_attr.encoder.slabBlockSkipFastRDO_mseFactor, 256.,
    "Scaling factor for the squared quantization step to "
    "compare with the MSE of the slab block, for fast decision "
    "on slab block skip.")

  ("integerHaar",
    params_attr.aps.lossless_flag, false,
    "Controls if lossless transform method is used:\n"
    " 0: no (RAHT)\n"
    " 1: yes (Integer Haar Transform)")

  ("rahtPredictionEnabled",
    params_attr.aps.rahtPredParams.intra_prediction_enabled_flag, true,
    "Controls the use of transform-domain prediction")

  ("rahtSubnodePredictionEnabled",
    params_attr.aps.rahtPredParams.subnode_prediction_enabled_flag, true,
    "Controls the use of transform-domain subnode prediction")

  ("rahtPredictionWeights",
    params_attr.aps.rahtPredParams.prediction_weights, {9,3,1,5,2},
    "Prediction weights for neighbours")

  ("rahtInterPredictionEnabled",
    params_attr.aps.inter_prediction_enabled_flag, false,
    "Controls the use of transform-domain prediction")

   ("rahtAveragePredictionEnabled",
    params_attr.aps.rahtPredParams.enable_average_prediction, true,
    "Controls the use of transform-domain average prediction")

  ("rahtMinWeightForModeSelection",
    params_attr.aps.rahtPredParams.min_weight_for_mode_selection, 64,
    "Maximum RAHT weight to activate average prediction instead of mode selection")

  ("rahtCrossChromaComponentPrediction",
    params_attr.aps.rahtPredParams.cross_chroma_component_prediction_flag, true,
    "Controls the use of CCCP (first chroma component to predict the second)")

  ("rahtCrossComponentResidualPrediction",
    params_attr.aps.rahtPredParams.cross_component_residual_prediction_flag, true,
    "Controls the of CCRP (luma residual to predict chroma residual)")

  ("rahtChromaPredictionModeLayerThreshold",
    params_attr.aps.rahtPredParams.chroma_pred_mode_layer_threshold, 3,
    "Number of bottom layers:\n"
    " - CCRP is active on (P frames),\n"
    " - before switching CCRP and CCCP priority if also enabled (I frames)")

  ("rahtRDOLambdaFactor",
    params_attr.encoder.rahtRDO_lambdaFactor, 7200.,
    "Scaling factor for lambda computation with accurate RDO for RAHT")

  ("rahtRDOQ",
    params_attr.encoder.useRahtRDOQ, true,
    "Controls the use of RDOQ with RAHT")

  ("rahtskipCoeffRDOQ",
    params_attr.encoder.useRahtskipCoeffRDOQ, true,
    "Controls the use of RDOQ for skip block decision with RAHT")

  ("qp",
    // NB: this is adjusted with minus 4 after the arguments are parsed
    params_attr.aps.init_qp_minus4, 4,
    "Attribute's luma quantisation parameter")

  ("qpChromaOffset",
    params_attr.aps.aps_chroma_qp_offset, 0,
    "Attribute's chroma quantisation parameter offset (relative to luma)")

  ("rahtQPperSlice",
    params_attr.aps.aps_slice_qp_deltas_present_flag, false,
    "Enable signalling of per-slice QP values")

  ("chromaFormat420",
     params_attr.encoder.abh.is420, false,
     "Chroma downsampling (only for RAHT, no Haar)")

  // This section is just dedicated to attribute recolouring (encoder only).
  // parameters are common to all attributes.
  (po::Section("Recolouring"))

  ("recolourSearchRange",
    params.encoder.recolour.searchRange, 1,
    "")

  ("recolourNumNeighboursFwd",
    params.encoder.recolour.numNeighboursFwd, 8,
    "")

  ("recolourNumNeighboursBwd",
    params.encoder.recolour.numNeighboursBwd, 1,
    "")

  ("recolourUseDistWeightedAvgFwd",
    params.encoder.recolour.useDistWeightedAvgFwd, true,
    "")

  ("recolourUseDistWeightedAvgBwd",
    params.encoder.recolour.useDistWeightedAvgBwd, true,
    "")

  ("recolourSkipAvgIfIdenticalSourcePointPresentFwd",
    params.encoder.recolour.skipAvgIfIdenticalSourcePointPresentFwd, true,
    "")

  ("recolourSkipAvgIfIdenticalSourcePointPresentBwd",
    params.encoder.recolour.skipAvgIfIdenticalSourcePointPresentBwd, false,
    "")

  ("recolourDistOffsetFwd",
    params.encoder.recolour.distOffsetFwd, 4.,
    "")

  ("recolourDistOffsetBwd",
    params.encoder.recolour.distOffsetBwd, 4.,
    "")

  ("recolourMaxGeometryDist2Fwd",
    params.encoder.recolour.maxGeometryDist2Fwd, 1000.,
    "")

  ("recolourMaxGeometryDist2Bwd",
    params.encoder.recolour.maxGeometryDist2Bwd, 1000.,
    "")

  ("recolourMaxAttributeDist2Fwd",
    params.encoder.recolour.maxAttributeDist2Fwd, 1000.,
    "")

  ("recolourMaxAttributeDist2Bwd",
    params.encoder.recolour.maxAttributeDist2Bwd, 1000.,
    "")

  ;
  /* clang-format on */

  po::setDefaults(opts);
  po::ErrorReporter err;
  const list<const char*>& argv_unhandled =
    po::scanArgv(opts, argc, (const char**)argv, err);

  for (const auto arg : argv_unhandled) {
    err.warn() << "Unhandled argument ignored: " << arg << "\n";
  }

  if (argc == 1 || print_help) {
    po::doHelp(std::cout, opts, 78);
    return false;
  }

  // set default output units (this works for the decoder too)
  if (params.outputUnitLength <= 0.)
    params.outputUnitLength = params.encoder.srcUnitLength;
  params.encoder.outputFpBits = params.outputFpBits;
  params.decoder.outputFpBits = params.outputFpBits;

  if (!params.isDecoder)
    sanitizeEncoderOpts(params, err);
  else
    sanitizeDecoderOpts(params, err);

  // check required arguments are specified
  if (!params.isDecoder && params.uncompressedDataPath.empty())
    err.error() << "uncompressedDataPath not set\n";

  if (params.isDecoder && params.reconstructedDataPath.empty())
    err.error() << "reconstructedDataPath not set\n";

  if (params.compressedStreamPath.empty())
    err.error() << "compressedStreamPath not set\n";

  // report the current configuration (only in the absence of errors so
  // that errors/warnings are more obvious and in the same place).
  if (err.is_errored)
    return false;

  // Dump the complete derived configuration
  cout << "+ Effective configuration parameters\n";

  po::dumpCfg(cout, opts, "General", 4);
  if (params.isDecoder) {
    po::dumpCfg(cout, opts, "Decoder", 4);
  } else {
    po::dumpCfg(cout, opts, "Coordinate system scaling", 4);
    po::dumpCfg(cout, opts, "Encoder", 4);
    po::dumpCfg(cout, opts, "Geometry", 4);
    po::dumpCfg(cout, opts, "Recolouring", 4);

    for (const auto& it : params.encoder.attributeIdxMap) {
      // NB: when dumping the config, opts references params_attr
      params_attr.desc = params.encoder.sps.attributeSets[it.second];
      params_attr.aps = params.encoder.aps[it.second];
      params_attr.encoder = params.encoder.attr[it.second];
      cout << "    " << it.first << "\n";
      po::dumpCfg(cout, opts, "Attributes", 8);
    }
  }

  cout << endl;

  return true;
}

//----------------------------------------------------------------------------

void
sanitizeEncoderOpts(
  Parameters& params, df::program_options_lite::ErrorReporter& err)
{
  // Input scaling affects the definition of the source unit length.
  // eg, if the unit length of the source is 1m, scaling by 1000 generates
  // a cloud with unit length 1mm.
  params.encoder.srcUnitLength /= params.inputScale;

  // global scale factor must be positive
  if (params.encoder.codedGeomScale > params.encoder.seqGeomScale) {
    err.warn() << "codingScale must be <= sequenceScale, adjusting\n";
    params.encoder.codedGeomScale = params.encoder.seqGeomScale;
  }

  // fix the representation of various options
  params.encoder.gbh.geom_stream_cnt_minus1--;
  for (auto& attr_sps : params.encoder.sps.attributeSets) {
    attr_sps.params.attr_scale_minus1--;
  }
  for (auto& attr_aps : params.encoder.aps) {
    attr_aps.init_qp_minus4 -= 4;
  }

  if (params.encoder.randomAccessPeriod < 1)
    err.error() << "randomAccessPeriod must be at least 1\n";

  // convert coordinate systems if the coding order is different from xyz
  convertXyzToStv(&params.encoder.sps);
  convertXyzToStv(params.encoder.sps, &params.encoder.gps);
  for (auto& aps : params.encoder.aps)
    convertXyzToStv(params.encoder.sps, &aps);

  // Trisoup is enabled when a node size is specified
  // sanity: don't enable if only node size is 0.
  // todo(df): this needs to take into account slices where it is disabled
  if (params.encoder.trisoupNodeSizes.size() == 1)
    if (params.encoder.trisoupNodeSizes[0] < 4)
      params.encoder.trisoupNodeSizes.clear();

  for (auto trisoupNodeSize : params.encoder.trisoupNodeSizes) {
    if (trisoupNodeSize < 4)
      err.error() << "Trisoup node size must be at least 4\n";
    if (trisoupNodeSize > 32)
      err.error() << "Trisoup node size must be at most 32\n";
  }

  params.encoder.gps.trisoup_enabled_flag =
    !params.encoder.trisoupNodeSizes.empty();

  // Certain coding modes are not available when trisoup is enabled.
  // Disable them, and warn if set (they may be set as defaults).
  if (params.encoder.gps.trisoup_enabled_flag) {
    if (params.encoder.gps.geom_duplicated_points_flag)
      err.warn() << "TriSoup geometry does not preserve duplicated points\n";

    params.encoder.gps.geom_duplicated_points_flag = false;
  }

  if (params.encoder.gps.trisoup_enabled_flag
      && (params.encoder.gbh.trisoup_QP < 0
        || params.encoder.gbh.trisoup_QP > 54))
    err.error() << "trisoupQuantizationQP must belong to the interval [0-54]\n";

  // Disable partitionning changes for Trisoup if Trisoup is not used
  if (!params.encoder.gps.trisoup_enabled_flag) {
    params.encoder.partition.safeTrisoupPartionning = false;
  }

  params.encoder.sps.localized_attributes_slab_thickness_minus1 = -1;
  params.encoder.sps.localized_attributes_slab_block_size_minus1 = -1;
  if (params.encoder.localized_attributes_slab_thickness <= 0) {
    err.warn() << "Localized attributes slab thickness must be greater than 0\n";
    params.encoder.localized_attributes_slab_thickness = 1;
  }
  params.encoder.sps.localized_attributes_slab_thickness_minus1
    = params.encoder.localized_attributes_slab_thickness - 1;

  if (params.encoder.localized_attributes_slab_block_size <= 0) {
    params.encoder.localized_attributes_slab_block_size
      = params.encoder.localized_attributes_slab_thickness;
  }
  params.encoder.sps.localized_attributes_slab_block_size_minus1
    = params.encoder.localized_attributes_slab_block_size - 1;

  if (!params.encoder.gps.inter_prediction_enabled_flag) {
    params.encoder.gps.gof_geom_entropy_continuation_enabled_flag = false;
  }

  params.encoder.sps.inter_frame_prediction_enabled_flag
   = params.encoder.gps.inter_prediction_enabled_flag;

  // ensure inter-frame trisoup parameters are correct
  params.encoder.sps.inter_frame_trisoup_enabled_flag =
    params.encoder.gps.inter_prediction_enabled_flag
      && params.encoder.gps.trisoup_enabled_flag;

  params.encoder.gps.trisoup_skip_mode_enabled_flag =
    params.encoder.gps.trisoup_skip_mode_enabled_flag
    && params.encoder.sps.inter_frame_trisoup_enabled_flag;

  params.encoder.gps.trisoup_early_skip_coding_mode_enabled_flag =
    params.encoder.gps.trisoup_early_skip_coding_mode_enabled_flag
    && params.encoder.sps.inter_frame_trisoup_enabled_flag;

  params.encoder.trisoup.alignToNodeGrid =
    params.encoder.trisoup.alignToNodeGrid
      || params.encoder.gps.trisoup_skip_mode_enabled_flag;

  params.encoder.sps.inter_frame_trisoup_align_slices_flag =
    params.encoder.sps.inter_frame_trisoup_enabled_flag
      && params.encoder.trisoup.alignToNodeGrid;

  params.encoder.sps.inter_frame_trisoup_align_slices_step =
    params.encoder.sps.inter_frame_trisoup_align_slices_flag ?
      lcm_all(
        params.encoder.trisoupNodeSizes.begin(),
        params.encoder.trisoupNodeSizes.end())
      : 1;

  // support disabling attribute coding (simplifies configuration)
  if (params.disableAttributeCoding) {
    params.encoder.attributeIdxMap.clear();
    params.encoder.sps.attributeSets.clear();
    params.encoder.aps.clear();
  }

  // fixup any per-attribute settings
  for (const auto& it : params.encoder.attributeIdxMap) {
    auto& attr_sps = params.encoder.sps.attributeSets[it.second];
    auto& attr_aps = params.encoder.aps[it.second];
    auto& attr_params = params.encoder.attr[it.second];

    // default values for attribute
    attr_sps.attr_instance_id = 0;
    auto& attrMeta = attr_sps.params;
    attrMeta.cicp_colour_primaries_idx = 2;
    attrMeta.cicp_transfer_characteristics_idx = 2;
    attrMeta.cicp_video_full_range_flag = true;
    attrMeta.cicpParametersPresent = false;
    attrMeta.attr_frac_bits = attr_sps.internalBitdepth - attr_sps.bitdepth;
    attrMeta.scalingParametersPresent = false;

    // Enable scaling if a paramter has been set
    //  - pre/post scaling is only currently supported for reflectance
    attrMeta.scalingParametersPresent = attrMeta.attr_offset
      || attrMeta.attr_scale_minus1 || attrMeta.attr_frac_bits;

    //// todo(df): remove this hack when scaling is generalised
    //if (it.first != "reflectance" && attrMeta.scalingParametersPresent) {
    //  err.warn() << it.first << ": scaling not supported, disabling\n";
    //  attrMeta.scalingParametersPresent = 0;
    //}

    if (attr_aps.lossless_flag) {
      if (attr_aps.init_qp_minus4)
        err.warn() << "QP set to 4 for lossless\n";

      if (attr_aps.aps_chroma_qp_offset)
        err.warn() << "QP chroma offset set to 0 for lossless\n";

      if (attr_aps.aps_slice_qp_deltas_present_flag)
        err.warn() << "no slice QP deltas for lossless\n";

      attr_aps.init_qp_minus4 = 0; // QP=4 for lossless
      attr_aps.aps_chroma_qp_offset = 0;
      attr_aps.aps_slice_qp_deltas_present_flag = false;
    }

    if (it.first == "reflectance") {
      // Avoid wasting bits signalling chroma quant step size for reflectance
      attr_aps.aps_chroma_qp_offset = 0;

      // There is no matrix for reflectace
      attrMeta.cicp_matrix_coefficients_idx = ColourMatrix::kUnspecified;
      attr_sps.attr_num_dimensions_minus1 = 0;
      attr_sps.attributeLabel = KnownAttributeLabel::kReflectance;
    }

    if (it.first == "color") {
      attr_sps.attr_num_dimensions_minus1 = 2;
      attr_sps.attributeLabel = KnownAttributeLabel::kColour;
      attrMeta.cicpParametersPresent = true;
    }

    // Assume that YCgCo is actually YCgCoR for now
    // This requires an extra bit to represent chroma (luma will have a
    // reduced range)
    if (attrMeta.cicp_matrix_coefficients_idx == ColourMatrix::kYCgCo) {
      attr_sps.bitdepth++;
      // TODO: do we need to increase internalBitDepth ? Probably yes,
      // or handle case where it becomes smaller than bitdepth.
      attr_sps.internalBitdepth++;
    }

    // Extend the default attribute value to the correct width if present
    if (!attrMeta.attr_default_value.empty())
      attrMeta.attr_default_value.resize(
        attr_sps.attr_num_dimensions_minus1 + 1,
        attrMeta.attr_default_value.back());

    auto& predParams = attr_aps.rahtPredParams;
    if (!predParams.intra_prediction_enabled_flag) {
      predParams.subnode_prediction_enabled_flag = false;
    } else {
      if (predParams.subnode_prediction_enabled_flag) {
        auto& weights = predParams.prediction_weights;
        if (weights.size() < 5) {
          err.warn() << "Five raht prediciton weights to be specified, "
                      << "appending with zeros\n";
          weights.resize(5);
        } else if (weights.size() > 5) {
          err.warn() << "Only five raht prediciton weights to be specified, "
                      << "ignoring others.\n";
          weights.erase(weights.begin() + 5, weights.end());
        }
        predParams.setPredictionWeights();
      }
    }

    if (attr_aps.dual_motion_field_flag) {
      // always true for dual motion: much faster and better coding perfs
      attr_params.motion.approximate_nn = TriStateEncParam::kTrue;
    }

    if (params.encoder.gps.inter_prediction_enabled_flag) {
      // add a qp offset equal to 0 for I and P if not set
      while (params.encoder.hierarchicalRAHTQP.size() < 2)
        params.encoder.hierarchicalRAHTQP.push_back(0);
    }

    if (!params.encoder.gps.inter_prediction_enabled_flag)
      attr_aps.inter_prediction_enabled_flag = false;

    if(!attr_aps.inter_prediction_enabled_flag) {
      attr_aps.dual_motion_field_flag = false;
      attr_aps.slab_block_skip_enabled_flag = false;
    }

    if (attr_aps.lossless_flag) {
      err.warn() << "no slab block skip with lossless\n";
      attr_aps.slab_block_skip_enabled_flag = false;
    }
  }

  // sanity checks
  if (
    params.encoder.partition.sliceMaxPoints
    < params.encoder.partition.sliceMinPoints)
    err.error()
      << "sliceMaxPoints must be greater than or equal to sliceMinPoints\n";

  for (const auto& it : params.encoder.attributeIdxMap) {
    const auto& attr_sps = params.encoder.sps.attributeSets[it.second];
    const auto& attr_aps = params.encoder.aps[it.second];

    if (attr_sps.bitdepth > 16)
      err.error() << it.first << ".bitdepth must be less than 17\n";

    if (attr_sps.bitdepth > 15 && attr_sps.params.cicp_matrix_coefficients_idx == ColourMatrix::kYCgCo)
      err.error() << it.first << ".bitdepth must be less than 16 with YCgCo\n";

    if (attr_sps.internalBitdepth > 16)
      err.error() << it.first << ".internalBitdepth must be less than 17\n";

    if (attr_sps.internalBitdepth > 15 && attr_sps.params.cicp_matrix_coefficients_idx == ColourMatrix::kYCgCo)
      err.error() << it.first << ".internalBitdepth must be less than 16 with YCgCo\n";

    if (attr_sps.bitdepth > attr_sps.internalBitdepth)
      err.error() << it.first << ".bitdepth must be less than, or equal to, internalBitdepth\n";

    if (attr_sps.internalBitdepth - attr_sps.bitdepth > 0
        && attr_sps.internalBitdepth - attr_sps.bitdepth < 2)
       err.warn() << "(0 < " << it.first << ".internalBitdepth - "
        << it.first << ".bitdepth < 2) is not efficient\n";

    if (attr_aps.init_qp_minus4 < 0 || attr_aps.init_qp_minus4 + 4 > 99)
      err.error() << it.first << ".qp must be in the range [4,99]\n";

    if (std::abs(attr_aps.aps_chroma_qp_offset) > 99 - 4) {
      err.error() << it.first
                  << ".qpChromaOffset must be in the range [-95,95]\n";
    }
  }

  params.encoder.singleSlabOptimizedImplementationEnabled =
    params.singleSlabOptimizedImplementationEnabled;
}

void
sanitizeDecoderOpts(
  Parameters& params, df::program_options_lite::ErrorReporter& err)
{
  params.decoder.singleSlabOptimizedImplementationEnabled =
    params.singleSlabOptimizedImplementationEnabled;
}

//============================================================================

SequenceEncoder::SequenceEncoder(Parameters* params) : SequenceCodec(params)
{
  // determine the naming (ordering) of ply properties
  _plyAttrNames.position =
    axisOrderToPropertyNames(params->encoder.sps.geometry_axis_order);
}

//----------------------------------------------------------------------------

int
SequenceEncoder::compress(Stopwatch* clock)
{
  bytestreamFile.open(params->compressedStreamPath, ios::binary);
  if (!bytestreamFile.is_open()) {
    std::cout << "Error: can't open output file!" << std::endl;
    return -1;
  }
  const int lastFrameNum = params->firstFrameNum + params->frameCount;
  for (frameNum = params->firstFrameNum; frameNum < lastFrameNum; frameNum++) {
    this->encoder.setInterForCurrPic(
      params->encoder.gps.inter_prediction_enabled_flag
      && ((frameNum - params->firstFrameNum) % params->encoder.randomAccessPeriod));
    if (compressOneFrame(clock))
      return -1;
  }

  std::cout << "Total bitstream size " << bytestreamFile.tellp() << " B\n";
  bytestreamFile.close();

  return 0;
}

//----------------------------------------------------------------------------

int
SequenceEncoder::compressOneFrame(Stopwatch* clock)
{
  std::cout << "Input frame number: " << frameNum << std::endl;

  std::string srcName{expandNum(params->uncompressedDataPath, frameNum)};
  PCCPointSet3 pointCloud;
  if (
    !ply::read(srcName, _plyAttrNames, params->inputScale, pointCloud)
    || pointCloud.getPointCount() == 0) {
    cout << "Error: can't open input file!" << endl;
    return -1;
  }

  // Sanitise the input point cloud
  // todo(df): remove the following with generic handling of properties
  bool codeColour = params->encoder.attributeIdxMap.count("color");
  if (!codeColour)
    pointCloud.removeColors();
  assert(codeColour == pointCloud.hasColors());

  bool codeReflectance = params->encoder.attributeIdxMap.count("reflectance");
  if (!codeReflectance)
    pointCloud.removeReflectances();
  assert(codeReflectance == pointCloud.hasReflectances());

  clock->start();

  if (params->convertColourspace)
    convertFromGbr(params->encoder.sps.attributeSets, pointCloud);

  scaleAttributesForInput(params->encoder.sps.attributeSets, pointCloud);

  // The reconstructed point cloud
  CloudFrame recon;
  auto* reconPtr =
    params->reconstructedDataPath.empty()
    && !params->encoder.sps.inter_frame_prediction_enabled_flag
    ? nullptr : &recon;

  auto bytestreamLenFrameStart = bytestreamFile.tellp();

  int ret = encoder.compress(pointCloud, &params->encoder, this, reconPtr);
  if (ret) {
    cout << "Error: can't compress point cloud!" << endl;
    return -1;
  }

  auto bytestreamLenFrameEnd = bytestreamFile.tellp();
  int frameLen = bytestreamLenFrameEnd - bytestreamLenFrameStart;
  std::cout << "Total frame size " << frameLen << " B" << std::endl;

  clock->stop();

  if (!params->reconstructedDataPath.empty())
    writeOutputFrame(params->reconstructedDataPath, {}, recon, recon.cloud);

  return 0;
}

//----------------------------------------------------------------------------

void
SequenceEncoder::onOutputBuffer(const PayloadBuffer& buf)
{
  writeTlv(buf, bytestreamFile);
}

//----------------------------------------------------------------------------

void
SequenceEncoder::onPostRecolour(const PCCPointSet3& cloud)
{
  if (params->postRecolorPath.empty()) {
    return;
  }

  // todo(df): don't allocate if conversion is not required
  PCCPointSet3 tmpCloud(cloud);
  CloudFrame frame;
  frame.setParametersFrom(params->encoder.sps, params->encoder.outputFpBits);
  frame.cloud = cloud;
  frame.frameNum = frameNum - params->firstFrameNum;

  writeOutputFrame(params->postRecolorPath, {}, frame, tmpCloud);
}

//============================================================================

SequenceDecoder::SequenceDecoder(Parameters* params)
  : SequenceCodec(params), decoder(params->decoder)
{}

//----------------------------------------------------------------------------

int
SequenceDecoder::decompress(Stopwatch* clock)
{
  ifstream fin(params->compressedStreamPath, ios::binary);
  if (!fin.is_open()) {
    return -1;
  }
  this->clock = clock;
  clock->start();

  PayloadBuffer buf;
  while (true) {
    PayloadBuffer* buf_ptr = &buf;
    readTlv(fin, &buf);

    // at end of file (or other error), flush decoder
    if (!fin)
      buf_ptr = nullptr;

    if (decoder.decompress(buf_ptr, this)) {
      cout << "Error: can't decompress point cloud!" << endl;
      return -1;
    }

    if (!buf_ptr)
      break;
  }

  fin.clear();
  fin.seekg(0, ios_base::end);
  std::cout << "Total bitstream size " << fin.tellg() << " B" << std::endl;

  clock->stop();

  return 0;
}

//----------------------------------------------------------------------------

void
SequenceDecoder::onOutputCloud(const CloudFrame& frame)
{
  clock->stop();

  // copy the point cloud in order to modify it according to the output options
  PCCPointSet3 pointCloud(frame.cloud);
  writeOutputFrame(
    params->reconstructedDataPath, params->preInvScalePath, frame, pointCloud);

  clock->start();
}

//----------------------------------------------------------------------------

int
SequenceDecoder::getOuputFrameNumber(const CloudFrame& frame)
{
  return frame.frameNum + params->firstFrameNum;
}


//============================================================================

double
SequenceCodec::outputScale(const CloudFrame& frame) const
{
  switch (params->outputSystem) {
  case OutputSystem::kConformance: return 1.;

  case OutputSystem::kExternal:
    // The scaling converts from the frame's unit length to configured output.
    // In terms of specification this is the external coordinate system.
    return frame.outputUnitLength / params->outputUnitLength;

  default:
    throw std::runtime_error("Unexpected value for OutputSystem");
  }
}

//----------------------------------------------------------------------------

Vec3<double>
SequenceCodec::outputOrigin(const CloudFrame& frame) const
{
  switch (params->outputSystem) {
  case OutputSystem::kConformance: return 0.;

  case OutputSystem::kExternal: return frame.outputOrigin * outputScale(frame);

  default:
    throw std::runtime_error("Unexpected value for OutputSystem");
  }
}

//----------------------------------------------------------------------------

void
SequenceCodec::writeOutputFrame(
  const std::string& postInvScalePath,
  const std::string& preInvScalePath,
  const CloudFrame& frame,
  PCCPointSet3& cloud)
{
  if (postInvScalePath.empty() && preInvScalePath.empty())
    return;

  scaleAttributesForOutput(frame.attrDesc, cloud);

  if (params->convertColourspace)
    convertToGbr(frame.attrDesc, cloud);

  // the order of the property names must be determined from the sps
  ply::PropertyNameMap attrNames;
  attrNames.position = axisOrderToPropertyNames(frame.geometry_axis_order);

  // offset frame number
  int frameNum = frame.frameNum + params->firstFrameNum;

  // Dump the decoded colour using the pre inverse scaled geometry
  if (!preInvScalePath.empty()) {
    std::string filename{expandNum(preInvScalePath, frameNum)};
    ply::write(cloud, attrNames, 1.0, 0.0, filename, !params->outputBinaryPly);
  }

  auto plyScale = outputScale(frame) / (1 << frame.outputFpBits);
  auto plyOrigin = outputOrigin(frame);
  std::string decName{expandNum(postInvScalePath, frameNum)};
  if (!ply::write(
        cloud, attrNames, plyScale, plyOrigin, decName,
        !params->outputBinaryPly)) {
    cout << "Error: can't open output file!" << endl;
  }
}

//============================================================================

const AttributeDescription*
findColourAttrDesc(const std::vector<AttributeDescription>& attrDescs)
{
  // todo(df): don't assume that there is only one colour attribute in the sps
  for (const auto& desc : attrDescs) {
    if (desc.attributeLabel == KnownAttributeLabel::kColour)
      return &desc;
  }
  return nullptr;
}

//----------------------------------------------------------------------------

void
convertToGbr(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud)
{
  const AttributeDescription* attrDesc = findColourAttrDesc(attrDescs);
  if (!attrDesc)
    return;

  switch (attrDesc->params.cicp_matrix_coefficients_idx) {
  case ColourMatrix::kBt709: convertYCbCrBt709ToGbr(cloud); break;

  case ColourMatrix::kYCgCo:
    // todo(df): select YCgCoR vs YCgCo
    // NB: bitdepth is the transformed bitdepth, not the source
    convertYCgCoRToGbr(attrDesc->bitdepth - 1, cloud);
    break;

  default: break;
  }
}

//----------------------------------------------------------------------------

void
convertFromGbr(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud)
{
  const AttributeDescription* attrDesc = findColourAttrDesc(attrDescs);
  if (!attrDesc)
    return;

  switch (attrDesc->params.cicp_matrix_coefficients_idx) {
  case ColourMatrix::kBt709: convertGbrToYCbCrBt709(cloud); break;

  case ColourMatrix::kYCgCo:
    // todo(df): select YCgCoR vs YCgCo
    // NB: bitdepth is the transformed bitdepth, not the source
    convertGbrToYCgCoR(attrDesc->bitdepth - 1, cloud);
    break;

  default: break;
  }
}

//============================================================================

const AttributeDescription*
findReflAttrDesc(const std::vector<AttributeDescription>& attrDescs)
{
  // todo(df): don't assume that there is only one in the sps
  for (const auto& desc : attrDescs) {
    if (desc.attributeLabel == KnownAttributeLabel::kReflectance)
      return &desc;
  }
  return nullptr;
}

//----------------------------------------------------------------------------

struct AttrFwdScaler {
  template<typename T>
  T operator()(const AttributeParameters& params, T val) const
  {
    int scale = params.attr_scale_minus1 + 1;
    return ((val - params.attr_offset) << params.attr_frac_bits) / scale;
  }
};

//----------------------------------------------------------------------------

struct AttrInvScaler {
  template<typename T>
  T operator()(const AttributeParameters& params, T val) const
  {
    int scale = params.attr_scale_minus1 + 1;
    // N.B. some rounding already done to internal bitdepth precision at end of
    // RAHT. If rounding is done here for 1 frac_bit it will introduce distorsion
    // and losses (RAHT values in the range [X + 0.25; X + 1.25) would be rounded
    // to X+1 because of successive roundings
    int roundOffset = params.attr_frac_bits > 1 ? (1 << params.attr_frac_bits - 1) : 0;
    return ((val * scale + roundOffset) >> params.attr_frac_bits) + params.attr_offset;
  }
};

//----------------------------------------------------------------------------

template<typename Op>
void
scaleAttributes(
  const std::vector<AttributeDescription>& attrDescs,
  PCCPointSet3& cloud,
  Op scaler)
{
  for (const auto& attrDesc : attrDescs) {
    if (!attrDesc.params.scalingParametersPresent)
      continue;

    auto& params = attrDesc.params;

    // Parameters present, but nothing to do
    bool unityScale = !params.attr_scale_minus1 && !params.attr_frac_bits;
    if (unityScale && !params.attr_offset)
      return;

    const auto pointCount = cloud.getPointCount();
    if (attrDesc.attributeLabel == KnownAttributeLabel::kReflectance)
      for (size_t i = 0; i < pointCount; ++i) {
        auto& val = cloud.getReflectance(i);
        val = scaler(params, val);
      }
    else if (attrDesc.attributeLabel == KnownAttributeLabel::kColour)
      for (size_t i = 0; i < pointCount; ++i) {
        auto& val = cloud.getColor(i);
        val = scaler(params, val);
      }
  }
}

//----------------------------------------------------------------------------

void
SequenceCodec::scaleAttributesForInput(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud)
{
  scaleAttributes(attrDescs, cloud, AttrFwdScaler());
}

//----------------------------------------------------------------------------

void
SequenceCodec::scaleAttributesForOutput(
  const std::vector<AttributeDescription>& attrDescs, PCCPointSet3& cloud)
{
  scaleAttributes(attrDescs, cloud, AttrInvScaler());
}

//============================================================================
