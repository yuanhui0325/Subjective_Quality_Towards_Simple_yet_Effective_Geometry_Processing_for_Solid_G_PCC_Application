General Information
===================

Reference software is being made available to provide a reference
implementation of a profile for geometry-based coding of solid dynamic
content for the G-PCC standard being developed by MPEG-3DGC (ISO/IEC SC29 WG7).

One of the main goals of the reference software is to provide a
basis upon which to conduct experiments in order to determine which coding
tools provide desired coding performance. It is not meant to be a
particularly efficient implementation of anything, and one may notice its
apparent unsuitability for a particular use. It should not be construed to
be a reflection of how complex a production-quality implementation of a
profile for geometry-based solid and dynamic content coding with a future
G-PCC standard would be.

This document aims to provide guidance on the usage of the reference
software. It is widely suspected to be incomplete and suggestions for
improvements are welcome. Such suggestions and general inquiries may be
sent to the general MPEG 3DGC email reflector at
<mpeg-3dgc@gti.ssr.upm.es> (registration required).

Releases notes
--------------

The GeS-TM (for **Ge**ometry-based **S**olid content **T**est **M**odel),
is being based on a subset of the G-PCC adoptions for coding dynamic solid
contents. The original implementation started from the experimental model
being studied in EE13.60, and is based on the `release-v20.0-rc1` of G-PCC
test model TMC13.

### `ges-tm-v12.0`

Tag `ges-tm-v12.0` will be released after reported issues are resolved,
if any. It may also include some code cleanups.


### `ges-tm-v12.0-rc1`

Tag `ges-tm-v12.0-rc1` includes adoptions made during the meeting #21 of the
WG7 held in Geneva (October 2025), as well as some code refactoring and cleanup.

- `geom/enc: m73984-7 - avoid redundant sort in octree`

    This commit removes redundant sort operation in encoder, in octree
    portion. It also includes minor code changes in decoder to explicitly
    disable the code that is specific to trisoup in the template code
    generation for octree (this could be already interpreted like this by
    the compiler optimizer but is more explicit now, since template
    parameter is directly used in the conditional blocks)

- `geom/motion: m73984-6 - fix motion for point cloud smaller than minPUSize`

    This commit allows signaling of one motion vector even if the point cloud
    is smaller than min PU size.

- `geom/enc: m73984-5 - fix msOctreeCurr and rootNodeSize misalignment`

    This commit fixes misalignment between `msOctreeCurr` and `rootNodeSize`
    octree decomposition, when `rootNodeSize` is artificially grown for
    aligning with trisoup grid.

    This has no impact on CTCs.

- `cli: (en/dis)able fast implementation for full frame slab block`

    This commit adds an option to allow enabling or disabling the faster
    implementation for full frame slab block of attributes. This allows
    comparing both implementations to ensure they are working the same.

    This option may be used in both encoder and decoder.

- `attr/ctc: m73984-4 - fix full-frame conditions`

    This commits sets `--localizedAttributesSlabThickness=4096` with trisoup
    for full-frame CTCs.
    This ensures the implementation optimized for full frame will be used
    at decoder.

- `attr/enc: m73984-4 - fix rootNodeSize and activation of full-frame implementation`

    This commit makes somme cleanups and fixes on use of the root node size:

    - `gbh.maxRootNodeSize` is derived from the octree depth and trisoup node
      size:
        ```
         gbh.maxRootNodeSize =
           trisoupNodeSise << gbh.tree_depth_minus1 + 1;
        ```
      and replaces the use of `(1 << maxRootNodeDimLog2)` as a finer
      approximation of the size of the point cloud for allocating buffers
      of slab blocks for localized attributes.

    - `gbh.rootNodeSizeLog2` is derived from the diced point cloud:
      - in decoder:
        ```
         gbh.rootNodeSizeLog2 =
           trisoupNodeSiseLog2 + gbh.tree_depth_minus1 + 1
        ```
      - in encoder, rather than later adding one, computation is fixed by
        using:
        ```
         gbh.rootNodeSizeLog2 =
           ceillog2(std::max(2,
             (_sliceBoxWhd.max() + nodeSizeCommon - 1) / nodeSizeCommon))
           + nodeSizeLog2Common; // dicing if any
        ```

    with these root node sizes, the decision for using full frame
    implementation is now derived properly and similarly in encoder and in
    decoder. In trisoup it is based on `maxRootNodeSize`, while in octree is
    based on `rootNodeSizeLog2`.

- `attr/enc: m73984-2 - start dual motion search from geometry motion`

    This commits makes the motions search for the largest PUs of dual motion
    to start from the motion vector found for the geometry.

- `attr/enc: m73984-2 - start motion search from parent motion`

    This commit makes motion search for a sub-PU to start from the best
    motion vector found for "parent" PU (i.e. the PU being split).

- `attr/enc/ctc: m73984-1 - fasterMotionSearch* parameters`

    This commit removes `--fasterMotionSearchForDual` option. It is always
    set to true. When set to false and using dual motion, coding is less
    good and much slower.

    This commit also changes the CTCs to always use `--fasterMotionSearch=1`
    for faster geometry encoding.

- `attr: m73919 - use up to 4 NN for MCAP`

    This commit uses up to 4 approximate nearest neighbors to improve
    the prediction quality of motion compensated attributes projection.
    A single neighbor is kept for lossless for complexity reasons.

- `raht/inter: m73913 - use default context initialization`

    This commit uses default value for initializing the entropy coding
    contexts for RAHT mode, avoiding a particular case with no coding gain.

- `raht/rdo: m73912 - use accurate RDO also with lossless`

    This commit allows using new accurate RDO method also with lossless
    where it was too complex. It thus also removes the old method that was
    still used for lossless and that was redundant. The additional
    complexity is removed by only tracking the modified context for mode
    selection, and adding some early decision based on the coefficients
    amplitude.

    The code is thus simplified, removing several things becoming unused.

- `enc/rdo: m73912 - optimize initialization in _RDO objects`

    In a `_RDO` object, created by a call to `makeRDO()`, when constructing
    the variables' stores `dataStart`, the first time `_RDO::start()` is
    called, and `dataEnd`, the first time `_RDO::startAlternative()` is
    called, the constructor of each variable was first called, when
    creating the tuple with a call to `make_tuple()`, then the tuple values
    where copied from the variables' pointers.

    In this commit, the move constructor is used to construct the variables'
    stores from a tuple directly constructed from the values pointed by the
    variables' pointers.

    It avoids a first initialization of all the variables from their default
    constructor, followed by a copy. It should help reducing the number of
    memory operations for big arrays, for instance.

- `raht/tidy: remove dead code`

    This commit removes unused variables, function parameters, and structure
    definition.

- `raht/inter: m67542/m74709 - contexts reduction for mode coding`

    This commit integrates the reduction on the number of contexts used for
    RAHT mode coding. These simplifications were followed as technology
    under consideration (TUC) until stabilization of the design for RAHT.

    These changes are integrated from proposals 2 and 3 reported in the TUC
    output document. Proposal 1 was already adopted.

- `trisoup/entropy: m74207 - reduce OBUF memory footprint`

    This commit reduces the memory footprint of OBUF for TriSoup.

    1) The number of OBUF instances is reduced to one per syntax element
       (presence, pos1, pos2, pos3 and pos4) by removing `mapIdx`.
       The length of secondary information `ctx2` for these elements is
       also reduced.

    2) In addition, the length of the secondary information for octree
       syntax elements `ctx2` is reduced when trisoup is active.

    Thus, the memory footprint of OBUF for both octree and trisoup, when
    trisoup is used, is roughly the same as the memory footprint of OBUF
    when only octree is used.

- `trisoup/enc: m74017 - fix dual centroid`

    This commit improves/fixes encoder derivation for dual centroid.

- `trisoup: m73979 - major cleaning/refactoring`

    This commit provides many simplifications, and cleaning of the trisoup
    code (see contribution for details). It also includes some variable
    renaming.

    All this work was done to simplify and align with the specification
    text while it has been written.

- `trisoup: m73979 - remove quality information`

    This commit removes computation and usage of quality information for
    reference frame and motion compensated frame. It avoids some complexity
    not bringing much and simplifies the spec.

### `ges-tm-v11.0`

Tag `ges-tm-v11.0` includes a bug fix for the handling of some software
parameters.

- `cli/fix: motion parameters handling`

    This commit fixes some input parameters that were wrongly named,
    or not properly set.

    Note: there is no impact on CTCs results.

### `ges-tm-v11.0-rc2`

Tag `ges-tm-v11.0-rc2` includes bug fixes, code cleanups, removal of unused
HLS, as well as renaming, refactoring and optimizations for trisoup made during
the preparation of the specification text.

- `trisoup/tidy: simplify second pass for centroids`

    Knowing that there are only two vertices in the second pass on centroid
    vertices, this commit simplifies and cleans it.

    This commit also avoids unnecessary copy and intermediate object
    in first pass and second pass.

- `trisoup/tidy: compact code in generateTrianglesInNodeRasterScan`

    This commit compacts some code in `generateTrianglesInNodeRasterScan` by
    making use of const reference.

- `trisoup: optimize std::map usage for triNodeEnable`

    This commit optimizes the current usage of a `std::map` for
    `triNodeEnable`. The search tree is used to count a number of vertex
    near of an edge extremity. An the information is later used to merge
    vertices. The tree was perpetually growing when advancing on trisoup
    nodes. However, when changing of trisoup node slice, the information
    for preceding preceding slice of nodes is not needed anymore.
    This commit thus erases from the tree the too hold information, thus
    leading to acceleration of the searches in the tree.

    The key type and value type are also optimized to consume less memory
    and thus reduce memory operations.

- `trisoup/tidy: use a struct in nonClosedNode`

    `nonClosedNode` was using an array of 6 integer to store some variables.
    Using an array for this is confusing and may lead to errors. In this
    commit, the array is replaced by a structure to explicitly name the
    variables and provide them a more suited type.

- `trisoup/tidy: do not use Vertex where not needed`

    This commit replaces several uses of `Vertex` by simple 3D coordinates
    stored in `Vec3<int32_t>`.

    `Vertex` is now only used within `findDominantAxis` function. The
    structure is moved from `geometry_trisoup.h` to `geometry_trisoup.cpp`,
    and it is put close to `findDominantAxis` function. Unnecessary operators
    are also removed from the structure.

- `trisoup/tidy: clean TrisoupCentroidVertex/cVerts`

    This commit performs some cleanups related to structure
    `TrisoupCentroidVertex` and its instances in `cVerts`.

    - `gravityCenter` is moved in the same structure instead of being stored
      separately but processed together with `cVerts`;
    - `drift` is removed, as it is not used;
    - no longer necessary uses of `std::move` are removed;
    - creation of some `cVerts` elements is cleaned.

- `trisoup: avoid a sort operation with non-closed surface`

    This commit avoids a call to `std::sort()` in the non-closed surface
    determination.

    Note: this change instroduces a minor drift in performances due to
    different rounding.

- `trisoup/tidy: cleaning/renaming/refactoring`

    This commit performs some code cleanups. It removes some dead code and
    redundant functions, a make the code more compact.

    The commit also includes some variables renaming to better align the
    code with the specification text.

- `trisoup: fix bits >= 5 of vertex codeword`

    This commit fixes the coding of the 5th and above bits of the vertex
    codeword for TriSoup.
    Context index could go out of the table. Coding 5 bits is for really
    particular configuration far from CTCs, and having really big context
    table is for this particular configurations is not especially useful.
    The context table is simply removed to fix the issue.

- `trisoup: fix init size of OBUF instances`

    This commit fixes the initialization size used for OBUF instances for
    TriSoup. More tables than needed were initialized.

- `trisoup/tidy: cleaning/renaming/refactoring`

    This commit performs some code cleanups and refactoring to get the code
    more compact and avoid repetitive computation of the same things.

    The commit also includes some variables renaming to better align the
    code with the specification text.

- `octree/tidy: remove unused contexts`

    This commit removes context unused since the integration of commit
    - `octree: m72862 - remove geometry scaling`

- `hls/attr/tidy: syntax grouping`

    It is cleaner to group together syntax elements relative to a same
    condition.

    This commit moves `aps.rahtPredParams.min_weight_for_mode_selection`
    in previous `if (aps.inter_prediction_enabled_flag) { ... }` block.

- `hls/geom/tidy: rename inter prediction flags`

    GPS and slice header had both a flag with same name:
    `interPredictionEnabledFlag`.

    In this commit the flag is renamed:
    - `ìnter_prediction_enabled_flag` for the GPS, and
    - `slice_inter_prediction_flag` for the GBH;
    thus avoiding any confusion and using hls naming conventions.

- `hls/geom/tidy: remove obsolete syntax`

    Occupancy atlas is not needed in GeS, and it has been removed a long
    time ago.

    This commit removes corresponding syntax that is remaining in the GPS
    and that is not used anymore.

- `raht/tidy: cleanups on pred mode`

    In this commit some unnecessary code and tests are removed.

    Also, the size of the context table for the coding of the prediction
    mode is properly derived, avoiding some waste of memory (and
    reducing memory traffic in RDO).

- `raht: fix issue #16 on context derivation for pred mode`

    `extreme` may take value in [0..3] interval, and `(cousinPredW > 5)` is
    boolean.

    In this commit, the context derivation is corrected accordingly.

- `inter: fix issue #15 on MSOctree depth`

    With small point clouds, MSOctree depth was set to 1 when it should be
    set to 0. It was introducing a non necessary particular case, leading
    to unexpected split of the root node, with the effect of increasing the
    runtimes.

    This commit fixes the issue.

    Fixing this issue reduced decoding runtime on attributes by up to 5% on
    some CTC configurations.

### `ges-tm-v11.0-rc1`

Tag `ges-tm-v11.0-rc1` includes adoptions made during the meeting #20 of the
WG7 held in Daejeon (June 2025), as well as some code refactoring and cleanup.

- `ctc/tidy: update cfg to not use deprecated parameter`  

    The `positionQuantizationScale` encoding parameter has become deprecated
    and `sequenceScale` should be used instead.

    In this commit, the cfg files are updated accordingly, and the
    deprecated parameter is removed.

- `ctc: m73761 - add octree raht inter local cfgs`  

    This commit adds configurations C1 and C2 for octree raht inter with
    local attributes.

- `enc: m72853/m72891 - add RDO tools and accurate bitrate estimation`  

    This commit add some tools to help in unifying and simplifying RDO in
    GeS-TM, with accurate bitrate estimation.

    One of the main principle is to effectively encode, in order to estimate
    bitrate from the true number of output bits in the bitstream, as well as
    from the decimal part of "unfinished" bits, and thus get accurate
    estimation, without having to rewrite rate estimation functions.

    There are also helpers to save any kind of variables (including contexts)
    and restore them if needed. And encoding alternatives may successively
    be tested, the best bitstream and entropy coder state being maintained.

    This method is useful for development, by not having to ensure and
    maintain RDO to be properly aligned with the effective coding.
    It will work well when not being to intensively used, or when the
    variable states are not to big in memory (saving them an restoring them
    may become costly otherwise).

    **Notes:**

    This commit includes original tools proposed in contribution m72853 and
    the improvements and cleanups from contribution m72891.

    The Improvements include better memory management by
    - avoiding some intermediate memory copies
    - deferring the construction of stored variables to when they are
      really needed.

- `cli: report input/output frame number`  

    This commits provides reporting of input (at encoder) and output
    (at decoder) frame number. This may help for debugging purposes and/or
    to easy statistics extraction.

- `cli/enc: reporting of motion field size`  

    In this commit, the reporting of motion field size for attributes
    (dual motion field) is cleaned to not use global variable, and to
    be reported for each attribute.
    In addition motion field size is now reported in bits, instead of
    bytes, and with accurate decimal estimation.

    This commit also adds reporting of motion field size for geometry.

- `motion/enc: m72853 - apply accurate RDO to motion search`  

    This commit applies accurate RDO tools within motion search, for
    accurate bitrate estimation of coded motion vector.

- `motion/tidy: m73761 - update motion search formulas`  

    In this commit, a color distortion factor is added and the RDO
    expressions are simplified. Presets are updated, in accordance with the
    simplifications, to keep the exact same performances.

- `motion/fix: m73761 - lambda formula for dual motion presets`  

    In the motion search presets for dual-motion with trisoup, the expression
    for computing the lambda value used in the Lagrangian RD optimization did
    not look correct: the expression for lambda was based on a quantization
    step computed from the opposite of the QP. Hence, a factor of the
    inverse of the effective squared quantization step was used, and thus,
    the rate was further penalized when quality was increasing.

    This commit fixes the lambda value by deriving a more classic
    expression, based on a factor of the squared quantization step.
    The expression is cleaned up and the lambda factor is slightly tuned.

- `motion/tidy: m73761 - remove unused parameters`  

    This commit removes obsolete parameters from `EncodeMotionSearchParams`
    that were not used anymore: `window_size`, `max_prefix_bits` and
    `max_suffix_bits`; as well as their derivation methods.

    The `K` parameter is also removed: the convergence does not rely on a
    maximum number of steps anymore.

- `motion/enc: m73761 - options and presets for motion search`  

    In this commit, the parameters in `EncodeMotionSearchParams` that were
    only set thanks to the choice of a preset derivation method are now also
    exposed as encoding parameters to let tests with motion search be more
    simple for the user.
    Encoder options are added for both geometry and attributes (with
    same name and "ForDual" suffix).

    In this commit, the number of presets is also reduced. The only kept
    presets are:
    - 0 (no preset),
    - 2 (presets derived for dense octree), and
    - 3 (presets deriver for dense trisoup).

    The presets' expressions are also simplified, implying a minor drift in
    the results for CTCs.
    The presets' for dual motion with octree (not used in CTCs because
    sub-optimal) are also now derived with same formulas as for trisoup.

- `attr/enc/tidy: remove dead code`  

    This commit deletes unused, outdated and not needed code for
    coefficients entropy estimation in RAHT. It removes the `struct
    PCCRAHTACCoefficientEntropyEstimate`. This structure contains tables
    which are also present in PCCResidualsEncoder. They are also removed.

- `attr/tidy: remove unused variables`  

- `attr: m72853 - add slab block skip for attributes`  

    This commit adds slab block skip for attributes. RDO for selecting the
    mode may be performed using accurate RDO (by default), or taking a fast
    encoding decision.

- `raht: m72879 - using lexicographic (aka raster scan) order`  

    This commit changes RAHT tree reduction and node processing order to
    become lexicographic order (also known as raster-scan order) instead of
    Morton code order.

    With this order, coding performances are slightly decreased, but the
    neighborhood of the nodes is much more simpler to get, avoiding several
    costly O(log2 N) searches for each node.

    By using simple iterators on the neighbors, some functions, like DC
    prediction, are also simplified, and thus accelerated.

    **Note:**

    It may be expected some additional runtime reduction by avoiding the
    naive call to `std::sort()` to order the points in lexicographic order:
    points should already be partially ordered from the geometry decoding
    order (octree and trisoup cases to be considered separately).

- `raht/enc: m72891 - use accurate RDO in RAHT`  

    This commit uses accurate RDO tools to simplify the readability and
    improve the RDOs in RAHT (mode selection and RDOQs).

    **Notes:**

    Chroma RDOQ is removed with this commit. It has not been ported
    has better gains were already obained without it.

    Mode selection in lossless is still using old method at is it more
    intensively used in lossless and using accurate RDO was increasing too
    much the encoding complexity.

- `raht/tidy: m72885 - cleanup intra/inter conditions`  

    This commit simplifies some conditions related to intra and inter
    prediction in RAHT. This has no impact on the performances (no normative
    change only simplifications based on actual values). Following the
    simplifications, the variables `enableIntraPred` and `enableIntedPred` can
    also be put as const.

    This commit thus makes the code more readable and more simple.

- `raht/hls/tidy: m72885 - remove unused rahtPredictionSkip1/prediction_skip1_flag`  

    In RAHT, intra prediction and inter prediction already became
    enabled provided there is more than one coefficient.
    The `prediction_skip1_flag` no longer has any impact.

    This commit removes `prediction_skip1_flag` from APS and associated
    `rahtPredictionSkip1` encoding parameter.

- `raht/tidy: m72885 - cleaner single point case handling`  

    Single point particular case was handled in RAHT by a specific syntax
    that was looking particularly ugly in the spec.

    This commit handles this particular case within the common RAHT
    framework by updating a few conditions.

- `geom/tidy: remove unused variables`  

- `octree: m72862 - remove QTBT`  

    This commit removes non maintained and inactive code for QTBT in octree,
    thus simplifying and accelerating the codec.

- `octree: m72862 - remove geometry scaling`  

    This commit removes the so called 'geometry scaling' in GPCC, which
    corresponds to adaptive quantization scale in the octree with a single
    point being rendered in the quantized node (not by a cuboid of the size of
    the node). This functionality was introduced for sparse content, and more
    specificaly for LiDAR. Its removal simplifies and accelerate the codec,
    and this feature has no meaning for dense content.

- `octree: m72862 - remove isPred`  

    This commit removes `ìsPred` variable from octree coding with OBUF.
    It had been introduced for sparse content and does not bring any value
    in Solid GPCC.

- `octree: m72863-16 - fix contexts for OBUF`  

    Item #16 of m72863, on fixing octree contexts for OBUF, has been agreed.
    This commit applies the bug fix.

- `octree/trisoup: m72883 - fix inter prediction and remove pred*UnComp*`  

    This commit fixes inter prediction in octree with trisoup diced point
    cloud, where prediction was wrong. Since the total removal of the
    prediction based on `pred*UnComp*` does not affect much the coding
    performances in octree only conditions, `pred*UnComp*` is totally removed.

- `octree/trisoup: m72903 - remove division in dicing with motion compensation`  

    This commit removes the division involved in the dicing of the motion
    compensated point cloud, for its use in the prediction of the occupancy
    in the occupancy tree part, when trisoup is used.

    Multipliers and fixed point precisions are carefully selected to get the
    exact same results as if the division was used.

- `trisoup: m72863-4 - remove non-cubic nodes`  

    Item #4 of m72863, on removing the unmaintained non-cubic nodes from
    GeS, has been agreed. This commits removes the associated HLS, and
    unmaintained code. It also cleans and simplifies the code by considering
    a single node width in each dimension.

- `trisoup: refactoring, cleaning and optimizing`  

    This commits refactors some TriSoup code and makes some cleanups.

    Cleaning is made for faceVertex, multi-node vertices, centroid, ...
    and some renaming is performed.

    `CentroidsResidualEncoder`/`CentroidsResidualDecoder` is introduced and
    the functions `GenerateCentroidsInNodeRasterScanEncoder`/\
    `GenerateCentroidsInNodeRasterScanDecoder` are split in two successively
    applied functions (2 passes) instead of calling two times the same
    function but performing two different things each time.

    This commit also fixes `cVerts` and `gravityCenter` under non-CTCs.
    New values are not `push_back` at the end of the lists, but directly
    set in existing list items, as for code tested under CTCs.

    This commit finally performs additional cleanups and optimizations of
    the code for TriSoup. It also removes some integer divisions:
    - for gCenter,
    - for QPedge, and
    - for vertex position quantization (non also used in decoder for early
      skip).

- `trisoup: m72863-11 - simplified condition for double centroid`  

    Item #11 of m72863, on simplifying the condition for activating the
    double centroid, has been agreed.
    This commit imposes the condition `trisoup_centroid_vertex_residual_flag`
    being true to allow double centroid.

- `trisoup: m72863-18 - fix merge vertex`  

    Item #18 of m72863, on fixing a bug in merge TriSoup vertex, has been
    agreed. This commit fixes the bug, with the effect of getting significant
    coding gains on TriSoup geometry.

- `hls/trisoup: m72863-10 - move skip HLS to header`  

    For item #10 of m72863, on moving to HLS the early skip threshold for
    determining vertices for a skipped node, it has been agreed to move the
    threshold syntax into the header. This commit add a separate encoder
    option for the early skipped nodes threshold and the encoder threshold
    used for non skipped nodes, and move the signaling of the early skipped
    one into the gbh.

- `hls/geom: m72863-3 - duplicated points`  

    Item #3 of m72863 suggesting to signal the duplicated points information
    relative to trisoup being not activated has been agreed. This commit
    applies this change.

- `hls/geom/ctc: m73009 - fix duplicated points`  

    Issue reported in m73009 also applies to GeS. In this commit, duplicated
    points flag is then used instead of unique point flag to align the spec
    and the software.

- `hls/attr: m72863-5 - hls disabled in lossless`  

    Item #5 of m72863, on moving HLS that is meaningless for lossless so
    that it is conditional to not being lossless, has been agreed.

    This commit thus changes the APS HLS to be coded relative to not using
    lossless.

    The code is also cleaned by moving and renaming the parameter
    `aps.rahtPredParams.integer_haar_enable_flag` into `aps.lossless_flag`.

- `hls/attr: m72863-7 - use single inter flag`  

    Item #7 of m72863, on removing redundancy of inter activation, has been
    agreed. This commit updates the APS to keep a single flags for
    activating inter-frame prediction.

- `hls/attr: m72863-15 - remove QPs per layer`  

    Item #15 of m72863, on removing QPs per layer, has been agreed.
    This commit removes the corresponding code for per layer QP derivation,
    the encoder parameters, and the HLS.

- `hls/attr: m72863-13 - remove flag deactivating inter`  

    Item #13 of m72863, on removing the flag for deactivating inter
    prediction, at the attributes slice header level, has been agreed.
    This commit removes the corresponding HLS.

- `hls/attr: m72863-14 - chroma HLS only if dim greater than 1`  

    Item #14 of m72863, on signaling chroma HLS in attribute data unit
    header (attributes brick header) only if attribute dimension is greater
    than 1, has been agreed. This commits adds the conditions to the
    corresponding HLS in abh.

### `ges-tm-v10.0`

Tag `ges-tm-v10.0` only updated revision number after the meeting #20 of the
WG7, as no additional issue had been reported or corrected.

### `ges-tm-v10.0-rc3`

Tag `ges-tm-v10.0-rc3` fixes the test conditions for lossless geometry coding
with local attributes.

- `ctc/local-attr: fix parameter for lossless`  

    This commit fixes the test conditions for lossless geometry coding with
    local attributes (octree-rath-local) where rahtMinWeightForModeSelection
    shall be equal to 1.

    This commit also removes from the yaml configuration files the setting
    of rahtMinWeightForModeSelection, where it is set equal to 64, as it is
    the default value.

### `ges-tm-v10.0-rc2`

Tag `ges-tm-v10.0-rc2` provides some cleaning of the geometry,
including refactoring of trisoup to separate the encoder from the decoder.

- `trisoup: split code into encoder and decoder`  

    The structures and functions for both encoding and decoding, that were
    mostly defined in `geometry_trisoup.h`, are split into common, encoder
    specific and decoder specific structures and functions.
    Things in common are put in `geometry_trisoup.(h/cpp)`, encoding stuff
    is put in `geometry_trisoup_encoder.(h/cpp)` and decoding stuff in
    `geometry_trisoup_decoder.(h/cpp)`.

    Several functions that were defined in `geometry_trisoup_decoder.cpp`
    have been been moved at the suited place.

    **Note:** there are still some encoding/decoding function that could benefit
    of some refactoring, in order to put things in common in
    `geometry_trisoup.(h/cpp)`, as they have more 70% of the code in common.

- `hls: SPS/GPS/APS - move syntax before extension flags`  

    This commit moves the syntax elements that were following the extension
    flags in SPS, GPS and APS. These elements are moved before the extension
    flags, and the flag is set equal to 0.

- `hls/attr: remove RAW attribute and transform type`  

    This commit removes the RAW attribute coding type and associated
    encoding/decoding.
    It also removes the signaling of the transform type as only RAHT
    transform is now supported.

- `hls/octree: remove IDCM`  

    IDCM is not supported from the early developments in GeS-TM. Some code
    and HLS was kept (mostly commented out) just in case it could be useful.
    As a cleaning of the solution while writing the specification of GeS,
    this commit removes all the syntax and code that was remaining for IDCM.

- `cli/doc: fix description for trisoupVertexThreshold`  

    Since the adoption of m71267, trisoupVertexThreshold is not anymore an
    encoder only parameter. This commit updates the documentation and the
    help of the command line interface to reflect this change.

- `enc: remove cmath pow(2, x) calls`  

    This commit removes the platform implementation dependent calls to
    std::pow(2, x), used in the encoder for RDO decisions.
    The fixed point approximation function, fpExp2(), is used instead.

- `tidy: remove unused files and fix files access mode`  

    This commit removes some files that were remaining from TMC13, and fixes
    the access mode to source files (removing execution flag)

### `ges-tm-v10.0-rc1`

Tag `ges-tm-v10.0-rc1` includes adoptions made during the meeting #19 of the
WG7 held online (March 2025).

- `entropy: m72077 - remove bypass stream`  
    This commit removes bypass stream, coming from GPCC-v1, as there is no
    real benefit of keeping it in ges.

- `raht: m71792 - improved quantization`  
    This commit adds the improved quantization for raht studied in EE13.63.

- `raht: m71793 - immproved average prediction`  
    This commit provides localization improvements for the average
    prediction in raht, as studied in EE13.70. Activation is not
    determined based on the depth anymore, but based on local information.

- `raht: m71953 - allow inter regardless of intra`  
    This commit allows inter prediction in RAHT regardless of the enabling of
    intra prediction.

    **Notes:**
    - Minor modifications were required to be compatible with the
    changes introduced by the localization of the average prediction.
   - Also to get similar results as in the contribution when comparing
    all intra with inter with no intra, the parameter
    `rahtMinWeightForModeSelection` must be set to significantly smaller
    value, for instance `8`, in order to enable inter on smaller trisoup nodes
    and use inter mode instead of null mode.

- `local-attr: m71852 - remove dedicated full frame processing`  
    This commit removes redundant and dedicated full frame attributes
    processing. Instead slab-block encompassing the point cloud can be used.

    CTCs are modified to use a slab-block size of 2048 for full frame
    attributes coding evaluation.
    Normative HLS is removed and localized attributes are always used.

    Non-normative optimization is added when detecting full frame processing
    is performed to avoid intermediate buffer management (which is costly
    for big slab block sizes).

- `motion/trisoup: m71794 - improve early skip efficiency`  
    This commit improves the performances of early skip for trisoup, by
    using better mv search, better RDO decision, and better derivation of
    trisoup parameters for neighboring nodes.

- `raht/ctc: m71795 - adjust QP for geom/attr rebalancing`  
    This commit updates the QP values for RAHT to balance geometry versus
    attributes bitrate by optimizing attributes versus total bitrate.

- `trisoup: m72162 - centroid generation improvement`  
    This commit fixes visual artefacts introduced when a centroid vertex is
    generated and makes a junction between the edges coming from two
    different surfaces that should not be linked together.

    **Notes:** Some minor changes were made by software coordinator
    - minor cleanups
    - encoder decision based on a distance computed using a square root
      has been replaced by a decision on the square of the distance, as
      suggested during the review of the contribution.

### `ges-tm-v9.0`
Tag `ges-tm-v9.0` only updated revision number after the meeting #19 of the
WG7, as no additional issue had been reported or corrected.

### `ges-tm-v9.0-rc2`

Tag `ges-tm-v9.0-rc2` includes some optimizations and cleanings of the code,
some bug-fixes and additional (non-normative) encoding parameters/options.

- `doc: fix documentation`  
    This commit removes an unintended copy/paste of an old text within the
    documentation.

- `enc/fix: remove dependencies to exp() and double computations`  
    This commit adds a fixed point approximation of the exp() function,
    and replaces some computations made in double by fixed point operations
    within RDO for early trisoup skip mode selection.

    Note: fpExp2() function is also added to perform fixed point
    approximation of exp2() function. It is used by fpExp() function to
    approximate exp() function.
    This functions are not intended to be used in normative parts of the
    codec and should only be used for encoder decisions.

    Note: with template parameter 32, for 32 bits fixed point precision in
    both input and output of the functions, the approximated exp() and
    exp2() should have a precision around the 4 or 5 most significant digits
    in decimal representation (~14-16 most significant bits accuracy)

- `raht/tidy: add encoding parameters for RDOQs`  
    This commit adds 3 new encoding parameters to enable/disable RDOQs in
    RAHT:
     - `--rahtRDOQ` is for historic RDOQ zeroing coefficients to improve
        coding efficiency, and
     - `--rahtSkipBlockRDOQ` is for RDOQ decision to skip coding of the
       coefficients of a block/node.

    This options may be usefull for particular investigations, to avoid some
    side effects and/or to avoid some traversal oder dependencies on the
    output quality (i.e adaptive schemes).

    An additional option is also added to control per level chromaRDOQ
    scaling factors: `--rahtChromaRDOQlevelDistScale`.

    Note: to remove traversal oder dependencies, for some particular testings,
    one could use for instance:
    `--rahtRDOQ=0 --rahtChromaRDOQ=0 --rahtSkipBlockRDOQ=0
     --rahtCrossChromaComponentPrediction=0
     --rahtCrossComponentResidualPrediction=0
     --rahtSubnodePredictionEnabled=0 --attribute=color`,
    and depending on the needs it may also be useful to disable intra
    prediction and/or average prediction.

- `motion/tidy: set pos0 in dual motion field`  
    Positions `pos0` in nodes of the motion field representation for dual
    motion was not set because it is not used/needed with current
    implementation.

    To avoid any bugs or mistakes in the future, e.g. if one wants to play
    with the motion fields, this commit restores the setting of the position
    of the nodes in the dual motion field's tree.

- `dec/fix: uninitialized hls value`  
    This commit fixes the use of uninitialized hls value under some
    configurations. This should have no impact on the decoding but avoids
    some memory analysis tools to complain.

### `ges-tm-v9.0-rc1`

Tag `ges-tm-v9.0-rc1` includes adoptions made during the meeting #18 of the
WG7 held in Geneva (January 2025), and some additional fixes and cleanups.
Also, following adoption of the localized attributes processing (processing of
the attributes within slab blocks), some dead code - kept an updated in case
the TUCs would have been discarded - has been removed.

- `local-attr/raht: m70746 - flexible slab(-block) size`  
    This commit introduces more versatility in the parameters defining slab-
    block size. It is now parameterized by two values: the slab thickness
    applying to X axis and the slab-block size, applying to both Y and Z.
    Thus smaller or wider blocks could be defined for a fixed thickness.

    This commit also removes the constraints on the necessity to align the
    slab thickness (and the slab-block size) onto the TriSoup node size;
    thus providing more flexibility for the choice of the slab-block size.

    Thus, powers of two slab-block sizes can now be used with any TriSoup
    node size, hence letting RAHT working on more equilibrated decomposition
    tree

    **Note:** adopted test conditions are different than test settings used in
      the EE13.70 and they are set in the next commit.

- `local-attr/raht: m70746 - new CTCs`  
    This commit adds new CTCs enabling local attributes with octree raht
    inter conditions for lossless geometry and lossless attributes coding,
    and with trisoup raht inter conditions (see resolutions for m70746).

    Two new experiment folders are generated by then `gen-cfg.sh` script
    when `--all` option is used:
    - octree-raht-inter-local
    - trisoup-raht-inter-local

    These folders can also be generated individually by adding `--local`
    option when calling `gen-cfg.sh` for a particular test condition.

- `attr: m70744 - increased internal bitdepth precision`  
    This commit allows increasing the internal bitdepth precision for
    attributes, including the reference frame inter-frame coding.

    The internal bitdepth is set by default to 16 bits for CTCs to ensure
    compatibility with up to 16 bits content.

    **Note:** internal bitdepth must be set equal to the input bitdepth for
    lossless coding.

    **Note:** it is not recommended to use and internal bitdepth increasing by
      only one or two bits as compared to the input bitdepth. Because of
      successive rounding operations it would introduce some noise and would
      decrease the coding performances.

    Original commit: c3f768315fc82aca06ee3bf54f06a568670d4285

- `motion: m70745 - fix contexts management`  
    This commit puts the same contexts management methods for the motion
    fields coding as for the geometry and attributes coding.

    The contexts can be reused between slab blocks (in case of dual-motion)
    and between frames, as already done for geometry and attributes, rather
    than using locally defined contexts.

    **Note:** It is expected that this better context management will allow to
      design more advanced contexts for motion field coding and will permit
      introducing better RDO methods.


    **Note:** as compared to original commit, the initialization of the contexts
      has been changed by using same initialization method as for attributes
      and geometry contexts. Thus if contexts are added in the future, they
      will be automatically reinitialized without having to modify the
      reset() function (and potentially forget it).
      This change has fixed the _ctxLocalMV not being reinitialized and
      that was leading to mismatches between full sequence encoding and per
      gop encoding.

 - `raht: m70728 - hierarchical QP for inter`  
    This commit adds encoder parameters to provide hierarchical QP
    for inter frame coding with RAHT.

    - Put hierarchical QP in encoding parameters and then in
      `abh.attr_qp_delta_luma` and `abh.attr_qp_delta_chroma`.

    - Deactivate hierarchical QP for all intra RAHT.

- `raht: m71271 - improve contexts for zero coefficient flag`  
    This commit improves entropy coding contexts for zero coefficient flag.
    Contexts are grouped according to the decomposition characteristics of
    the 8 raht subbands, thus improving coding efficiency and reducing the
    number of contexts.

- `trisoup/fix: m71267 - align decomposition tree an MSOctree`  
    This commit fixes encoder misalignment between decomposition octree and
    the MSOctree built on the current frame for motion search (root node
    could be at a different level with non-power of 2 trisoup node size).
    There is no impact on existing tools, but it will be useful for skip
    mode RDO in next commit.

    The only difference with test conditions is that with non power of two
    trisoup node size, the octree generated by the encoder could contain all
    the points of the point cloud within its first child node.
    With the provided change, this first child node becomes the root node.

    **Note:** The main purpose of mSOctreeCurr with current tools is to preorder
      the points for motion search which has been accelerated by assuming
      the points of the current frame are already partitionned in Morton
      order when spitting a node of the motion field in the RDO process.
      With the skip mode the tree will also have to be aligned with
      occupancy tree to simplify the extraction of nodes information in the
      RDO decision.

- `trisoup: m71267 - skip coding mode`  
    This commit introduces skip coding mode for trisoup.

    The original code provided by the proponent has been significantly
    revised by the software coordinator during the integration.
    Several cleanup have been performed with no observed impact on the
    results (outside of speeding up the execution):
    - Some statistics computed on the current frame for skip mode RDO were
      also computed on the reference frame, in both encoder and decoder, but
      not used. Computation has been removed when not necessary.
    - These encoder specific computation for skip RDO were performed within
      the constructor of MSOctree (motion search octree, also used for
      motion compensation). They have been moved into a dedicated function.
    - Some of the computation made in double have been simplified to use
      integer arithmetic in the RDO process.
    - Skip encoder parameters have been separated from the skip HLS and
      moved outside of the gps to be put in encoder specific structure.
    - Hashtables and code used to supposedly remove duplicated points have
      been removes, has no duplicated point have been observed.
    - Hashtable storing the vertices for skipped blocs has been replaced by
      a small ring buffer, as the vertices are processed in the same order
      as they are created.

    An additional change has been performed to simplify the code, with
    minor impact on the results:
    - Instead of processing the points of the skipped blocks within the
      flush2PointCloud() function the points have been put within the
      renderedBlock before the call to the function, thus processing them as
      other points rendered by the rasterization.
    - It has an impact on the results (mostly noise) because the points
      in renderedBlock are ordered in lexicographic order. Any duplicates
      that could occur during the rendering of a trisoup node are removed
      here.
    - Doing like this allowed to fix the support for the slab blocks,
      which was not provided by the original source code, while avoiding
      additional code duplication within flush2SlabBuff() function.

    **Note:** The coding of thVertexDetermination which was originally an
      encoding parameter is coded in m71267 as a 2-bit integer number with
      no check/limitation being performed on the encoding parameter.
      Also, no study showing if a higher value than 3 could not be
      beneficial has been released. So, it might be better to use a non
      fixed bit size representation here (e.g. Golomb-Rice code) or to add
      checks at the encoder.

### `ges-tm-v8.0`

Tag `ges-tm-v8.0` includes some minor bug fixes for particular situations
not occurring under CTCs.

### `ges-tm-v8.0-rc2`

Tag `ges-tm-v8.0-rc2` includes some optimizations and cleanings of the code,
some bug-fixes and the integration of m69536.

- `enc/raht: m69536 - chroma RDOQ`  
    This commit add specific RDOQ decision for chroma component within the
    RAHT encoder. It is slighlty different from original proposal, following
    the integration phase on top of the simplifications and improvement of
    RAHT.

### `ges-tm-v8.0-rc1`

Tag `ges-tm-v8.0-rc1` includes most of the adoptions made during the meeting
#17 of the WG7 held in Kemer (November 2024), and some additional fixes and
cleanups. It also includes an update of the scripts for output log analysis in
order to extract geometry and attributes runtime for the reporting under
new CTCs.

The contribution m69536 for chroma RDOQ will need to be revisited for a better
integration with the other RAHT improvements and will be provided with
`ges-tm-v8.0-rc2`.

- `entropy: m69533 - fix bypass coding defect`  
    Under some conditions, bypass coding could fail on beeing decoded
    properly.
    When used really early in the bitstream, high bound may be higher than
    0x8000 and bypass coding may discard the significant bit.

    This is fixed by initializing the coding range with $[0, 0.5)$ interval
    instead of $[0, 1)$ (i.e. ~1 bit is lost in the stream).

    **Note:** This fix is integrated to ensure bitstream can be properly
    decoded. Possibility to use a better range interval will be studied
    after the release of `ges-tm-v8.0-rc1`.

- `raht/encoder: fix duplicated points handling`  
    Raht encoder could crash in presence of duplicated points.
    The coefficients buffer size is updated to remove duplicates.

- `raht: m69439 - massive simplification of RAHT`  
    This commit fixes RDO computation formula and removes tools that
    were more compensating the wrong RDO formulation rather than bringing
    something efficient. And they were introducing multiple passes RDO.
    Layer mode selection, intra layer, inter layer and upper layer mode are
    removed.

    The raht prediction thresholds 0 and 1, that were mostly used in tmc13
    for non-dense content, are also removed.

- `raht: m69440 - localization of coefficients coding`  
    This commit provides localization of RAHT coefficients coding/decoding
    within the RAHT transform to remove the two pass coding scheme.

- `raht: m69574 - improved coefficients coding`  
    This commit finishes the localization of RAHT coefficients coding by
    removing the runlength syntax. It brings local coding syntax to
    compensate for the runlenth removal and to provide additional
    compression gains.

- `raht: m69843 - 420 chroma sampling`  
    This commit provides 420-like chroma sampling for RAHT to lower the
    complexity and provide some compression gains.
    The feature could be used for lower complexity profiles, but will be
    kept deactivated in the CTCs.

- `tuc/local-attr: m69453 - using slab blocks for localized attributes`  
    This commit further split slab attributes processing into slab blocks
    to reduce the number of passes within a complete slab. It provides even
    more localized processing for attributes.
    The motion vector field for dual motion is currently signaled at slab
    block level.

- `tuc/local-attr/raht: m69453 - remove RAHT per block`  
    This commit removes RAHT per block which is not needed anymore.

- `tuc/local-attr/raht: m69498 - use boundaries' coeffs for RAHT intra prediction`  
    This commit stores the RAHT coefficients on the external boundary of a
    block, retrieve them for a next block, and use them with RAHT intra
    prediction to perform similarly as if there was no block. It
    corresponds to what was originaly suggested to study within m68790.
    In addition information on boundary for contexts formation of m69574 is
    also stored.

- `tuc/local-attr/raht: m70339 - enable root node intra DC prediction`  
    This commit enables intra DC prediction from the neighboring slab blocks
    boundaries in RAHT.

- `enc: m70412 - MV search improvements`  
    This commit fixes the R/D cost formula by using a mathematically more
    suited Lagrange cost.
    It also provides some fixes and optimizations, as well as new presets.

- `trisoup: m69835 - reconstructing non-closed surfaces`  
    This commit introduces improvements for the representation of non-closed
    surfaces within Triangle Soup framework.

- `trisoup: m69926 - triangle construction improvement`  
    This commit improves visual quality of constructed triangles in TriSoup.

- `raht: m69535 - chroma component prediction harmonization`  
    This commit harmonizes the chroma component prediction
    processes/implementations for:
     - the cross-chroma component prediction (prediction of Cr from Cb), and
     - the cross component prediction (prediction of Cb and Cr from Y).

- `raht: m69931 - enable chroma component predition on top layers`  
    This commit modifies the activation conditions for the chroma component
    prediction so that it can be efficiently activated on upper layers.
    It also adjusts the semantics and signaling conditions for the number
    of layers CCRP is active.

- `tools: update scripts for geom/attr runtime reporting`  
    This commit provide an updated version of the scripts for encoder/decoder
    output log analysis. The update allows extraction of the respective
    geometry and attributes processing runtime, in addition of the already
    extracted overall runtime.

    These scripts may be used to extract and fill the runtime expected in the
    new version of the excel sheet templates for reporting results under CTCs.

### `ges-tm-v7.0`

Tag `ges-tm-v7.0` includes some cleanups of the code as well as some fixes for
cases that could occur outside of the CTCs and mostly when the technology
under consideration for localized attributes processing by slabs was active.

### `ges-tm-v7.0-rc2`

Tag `ges-tm-v7.0-rc2` includes a minor fix for compilation issues, and major
refactoring of RAHT:

- `raht: major refactoring`  
    This commits continues the major refactoring of the RAHT code.
    It continues the work started in `ges-tm-v6.0-rc2`

    Most of the changes that were made in the decoder have been applied to
    the encoder in a similar way, when possible, to fit the RDO nature of the
    encoder.

    Some of the changes are reported bellow:
     - encoder uses similar bottom-up process as already implemented in
       decoder
     - inter predictor computation is moved into bottom-up process (as
       geometry of the inter predictor is now always the same as the
       reconstructed geometry of the current frame)
     - translateLayer and interTree are removed, since now useless
     - bottom-up process is now templatized

    The encoder's RDO process has been optimized by removing most of the
    Lagrangian costs computations made in double, they are now performed
    using integers. The fixed point fpLog2 approximation of log2
    -- which was only applied to probabilities, i.e to the (0; 1] interval,
       for entropy/bits cost computation --
    has been replaced by a linear approximation from the tabulated interval.

    Other RDO/encoder optimizations have been applied, like:
     - storing of pred mode in child nodes,
     - getting the prediction mode,
     - ...

    Various additional code compactions, optimizations and cleanups have
    been applied to either the encoder, the decoder or both.

### `ges-tm-v7.0-rc1`

Tag `ges-tm-v7.0-rc1` includes adoptions made during the meeting #16 of the
WG7 held in Sapporo (July 2024), and some additional fixes and cleanups.

- `trisoup/tidy: refactor to use dedicated node structure`  
    This commit replaces the use of "big" PCCOctree3Node structure within
    trisoup by the use of two new node structures respectively dedicated to
    the encoder and to the decoder. These structures only contains necessary
    information. This is achieved through templatization.

    This commit also includes minor cleanups.

- `attr/tidy: build and encode dual motion field in raster scan order`  
    makes dual motion MV field being encoded in same order as for geometry

- `attr/tidy: avoid generating full depth MSOctree`  
    In case octree geometry is used, it is not necessary to build a full
    depth MSOctree. This commit limits the depth to reach only the
    min_pu_size in case trisoup is not used.

- `attr/tidy: replace PUTrees with MVField`  
    Occupancy tree must be known to properly use PUTrees, this is not
    convenient for reusing geometry motion field with attributes.

    This commits replaces the PUTrees, used for motion field
    representation, with a more versatile and ease to use MVField
    representation.

- `attr/tidy: clipping compensated point cloud to slab boundaries`  
    This commit should slightly improve the prediction for localized
    attributes and accelerate its use in inter RAHT by clipping the
    compensated point cloud to fit within the slab boundaries.

- `octree/tidy: refactor to use common raster scan`  
    In this commit, Octree is refactored to use a common raster scan
    framework with dual motion field which makes the code a bit cleaner.

- `attr: m68327 - allow using mcap without using dual motion`  
    In this commit, HLS and code structure is changed to allow using motion
    compensated attributes projection even if dual motion is not enabled.

    This commit also avoids unneeded sort for motion compensated
    attributes.

- `attr/tidy: do not copy current attributes for mcap`  
    Attributes from reconstructed point cloud do not need to be copied with
    geometry before performing motion compensated attributes projection.
    This commit simplifies the code accordingly.

- `attr: m68327 - always use mcap`  
    As agreed during the meeting, motion compensated attributes projection
    (mcap) will always be used to simplify the codec and keep best coding
    performances. This commit enforce always using mcap and removes HLS
    previously used to enable it.

- `attr/tidy: remove not needed attributes compensation`  
    Now that mcap is always used, it in not needed to compensate the
    attributes when building compensated point cloud for geometry.

- `raht/tidy: simplify decoder translate layer`  
    Since mcap is now always used, num points, positions and weights are
    now the same in the compensated frame and in current frame.

    Decoder is simplified accordingly.

- `raht: m68283 - cross component residual prediction`  
    This commit uses buffer to dynamically maintain correlation information
    between the RAHT luma and RAHT Cb/Cr residues in order to derive a
    prediction of the chroma components residues with a linear model of the
    luma residues.

- `raht: m68397 - improve average prediction`  
    This commit improves inter frame average prediction by also taking into
    account the grand parent mode to determine average prediction weights.

- `geom/trisoup: m68255 - local quality units`  
    This commit introduces local quality units for trisoup.

- `trisoup: m68285 - reduce OBUF memory usage`  
    This commit reduces the memory usage when coding the Trisoup vertex
    position bits using OBUF scheme.

    This is scheme 1 of m68285.

- `trisoup: m68350 - memory reduction over m68285`  
    This commit further reduces memory usage and also takes care of memory
    allocation to further improve runtimes.

- `raht/tidy: remove FixedPoint and VecAttr`  
    FixedPoint was hiding many conversion and shift operations.
    These operations were progressively made visible again. Now, this
    commit totally removes the use of FixedPoint type.

    VecAttr is also removed, as the dynamic allocation nature of the vector
    it uses is not really justified.

- `raht/tidy: uniformization of buffers definition`  
    Some local buffers in raht are still defined as double entry arrays.
    Since single entry array convention is currently used for the main
    buffers, this commit uses the same notations for intermediate buffers.

- `raht/tidy: simplify encoder translate layer`  
    This commit simplifies encoder side as previously made for the decoder,
    since motion compensated projection is now always used.

- `enc: restrict slab thickness`  
    This commit restricts the slab thickness to be a multiple of any
    trisoup node sizes. Other values may provide unexpected results.

- `raht/tidy: use attr_t instead of int`  
    This commit makes using vectors of attr_t (i.e. uint16_t) instead of int
    as input/output of raht.

### `ges-tm-v6.0`

Tag `ges-tm-v6.0` has been released after some cleanups and a few minor bug
fixes (outside of CTCs) have been integrated.

### `ges-tm-v6.0-rc2`

Tag `ges-tm-v6.0-rc2` includes a few bug fixes and a refactoring of the RAHT.

- `raht: fix issue #2 - wrong condition used to access a table`:  
    Encoder was crashing when using inter frame coding with encoder parameter
    --rahtLayerInterModeSelectionEnabled=0, because it was accessing a non
    allocated table element. The table is not supposed to be used.

    A condition enabling the use of the table was based on the wrong variable
    name. This commit fixes the condition.

- `raht: fix issue #4 - RDO could select intra layer mode on root node`:  
    Encoder could select intra layer mode on root node, leading to mismatch
    between encoder and decoder output in some situations.
    This commit fixes this issue by enforcing infinite RD cost for intra
    layer mode on first layer.

- `raht: major refactoring`:  
    This commit provides a major refactoring of the RAHT code:

    - The code is split between the encoder and decoder into two
      separate files:
      - move RAHT.cpp -> RAHTEncoder.cpp
      - add RAHTDecoder.cpp
    - Many optimization are introduced within the decoder
    - Some are reported into the encoder, but the work on the encoder will be
      continued after next meeting cycle

- `enc: fix issue #5 - do not use std::log2 for RDO`:  
    This commit removes dependencies on std::log2() for the encoder when
    performing Rate Distortion Optimizations (RDO). Instead, a fixed point
    approximation of log2 pcc::fpLog2 is introduced so that the results will
    not depend on the platform and its implementation dependent math library.

    This commit thus changes computations for RAHT RDO, as well as for motion
    search, and should solve issue #5.

- `misc: fix issue #3`:  
    Theres is a bug in gcc-8.1/8.2 (fixed in 8.3). This commit avoids the
    compiler complaining for no reason.

### `ges-tm-v6.0-rc1`

Tag `ges-tm-v6.0-rc1` includes adoptions made during the meeting #15 of the
WG7 held in Rennes (April 2024), and some additional fixes and cleanups.

- `raht: tidy - remove unnecessary operations`:  
    Performs some cleanups of the code for RAHT.

- `octree: refactor - restore a single function for enc/dec`:  
    Refactors octree's encoding and decoding functions.
    A single function is restored to support both octree and trisoup.

- `geom/attr: refactor - payload handling structure for single pass decoding`:  
    Refactors the decoding process to allow expected feature of
    being able to perform geometry and attributes decoding together in a
    single pass (with more local processing of the data).

- `geom/attr: tidy - dual motion parameters derivation and HLS`:  
    Provides general cleanups and fixes for the dual motion's
    parameters and HLS.

- `geom/attr: refactor - simplify parameters passing to enc/dec functions`:  
    Simplifies parameters passing to several Encoding/Decoding
    functions.

- `geom/attr: tidy - minor cleanups`:  
    Provides minor code cleanups.

- `geom/attr: refactor - template octree enc/dec with/without trisoup`:  
    Uses templatized versions for the octree encoding/decoding
    functions, to obtain better execution performances.

- `trisoup: m67002 - versatile quantization`:  
    Introduces versatile geometry quantization of TriSoup to provide
    better rate control with quantization parameter.

- `trisoup: misc - optimize OBUF memory use`:  
    Optimizes OBUF memory usage and management to improve
    execution speed.

- `enc/trisoup: m67002 - update motion presets`:  
    Improves presets for inter-frame motion search and
    prediction with trisoup.

- `trisoup: m67015 - support for non-powers of two node sizes`:  
    Provides finer granularity for quality control of the
    TriSoup model by introducing the support for TriSoup node sizes which
    are not powers of two.

- `trisoup: m67539 - improve trisoup surface reconstruction`:  
    Improves trisoup surface reconstruction around vertices.
    It reduces visual artefacts when the surface is passing through the
    corner using one vertex instead of multiple vertices.

- `trisoup: m68006 - simplify face vertex`:  
    Provides several cleanups and simplifications of the
    face vertex for TriSoup.

- `trisoup: m67544 - reduce OBUF memory usage`:  
    Refines the OBUF scheme of Trisoup vertex presence flag.

- `octree: m67540 - improve inter states contextualization`:  
    EE13.60-Test2a which improved inter-frame
    contextual modeling for occupancy coding.

- `octree: m67004 - improve inter prediction`:  
    Combination 1 of m67004 with m67540.

- `enc/raht: m67627 - RDOQ improvement`:  
    The reconstruction distortion is also taken into account
    within RDOQ for RAHT.

- `enc/raht: m67627 - fix RDO with lossless`:  
    Fixes the RDO issue when determining prediction mode under CW
    condition.

- `raht: m67627 - fix prediction mode signaling`:  
    Modifies the software to match with the spec when signalling
    the prediction mode.

- `raht: m67628 - remove redundant signalling and RDO calculations`:  
    Removes unnecessary signaling for the prediction mode of
    RAHT, and avoids corresponding RDO calculations.

- `raht: m67541 - cross chroma component prediction`:  
    Introduces cross chroma component prediction, for intra
    frame coding. The first chroma component is used to predict the second
    chroma component.

- `raht: m67542 - remove unused contexts`:  
    Removes unused contexts for prediction Mode coding.

- `raht: m67406 - refactor div/mod 3 RAHT level/layer`:  
    Provides an example of implementation to avoid using
    divisions by 3 and modulo by 3 when computing levels and layers in RAHT
    transform.

- `raht: m67406 - use fixed point for attributes`:  
    Refactors RAHT to use fixed point representation for the
    intermediate values of the attributes.

- `misc: m67406 - approximate division for FixedPoint representation`:  
    Introduces a new division approximation function to be used
    for divisions removal and uniformization with other existing division
    approximations.

- `raht: m67406 - share a single approximate division`:  
    Removes two LUTs and increases precision when computing
    weighted intra/inter prediction.

- `raht: m67406 - rounding inter integer haar predictor`:  
    Slightly improves the prediction by rounding the integer
    haar predictor.

- `tuc/raht: m67093 - RAHT per block`:  
    Adding possibility of performing RAHT transform by blocks as a technology
    under consideration.

- `tuc/geom/attr: m67808 - localized attributes`:  
    Add localized attributes processing as a technology under consideration.
    Localized attributes improve the design of the codec to let it be more
    friendly to being embedded in a device.

- `enc/trisoup: m67838 - fix slice/tile alignment`:  
    From discussions about m67838, try fixing slice/tile alignment in case
    multiple slices/tiles are used, to remove artefacts on slice boundaries.

- `enc/trisoup: m67601 - improve vertex determination`:  
    Introduces an encoding parameter to control the threshold
    for vertex determination. This threshold is put to 1 by default
    (instead of 0).

- `enc/trisoup: m67017 - new CTCs with 6 rate points`:  
    New Common Test Conditions for TriSoup.

- `enc/trisoup/raht: m67018 - new CTCs for colors with trisoup`:  
    New Common Test Conditions for colors with RAHT coded on top of
    a TriSoup geometry.

### `ges-tm-v5.0`

Tag `ges-tm-v5.0` has been released after some cleanups and improvements and
a few minor bug fixes have been integrated. The `README.md` file has also been
updated to advise building for Release, for runtime comparisons.

### `ges-tm-v5.0-rc1`

Tag `ges-tm-v5.0-rc1` includes adoptions made during the meeting #14 of the
WG7 held online (January 2024), and some additional fixes and cleanups.

- `tidy: motion - remove unused AttributeInterPredParams::motionVector`:  
    Removes `AttributeInterPredParams::motionVector` and
    `tmc3/attr_tools.h:MotionVector` which are not used.

- `geom/tidy: motion - remove unused variables and parameters`:  
    Removes `LPUnumInAxis` variables and `numLPUPerLine` parameter which are not
    used.

- `geom/tidy: octree - remove unused PCCOctree3Node::siblingOccupancy`:  
    Removes `PCCOctree3Node::siblingOccupancy` which is not used.

- `geom/tidy: octree - disable remaining IDCM code`:  
    - Comments out some IDCM related coding/decoding functions
    - Comments out `PCCOctree3Node::numSiblingsPlus1` and
      `PCCOctree3Node::numSiblingsMispredicted` which are only used by IDCM
    - Comments out IDCM coding contexts
    - Removes remaining angular mode related coding contexts

- `geom/trisoup: m65829 - Single pass geometry coding`:  
    - two octrees
    - one pass geo for octree-trisoup
    - context class modified to handle interleaved octree and TriSoup
      geometry coding

- `geom/trisoup: m65831 - centroid residual quantization`:  
    Improves centroid residual quantization.

- `attr/fix: m66895 - point cloud should not be offset`:  
    Point cloud and compensated point cloud should not be offset
    before and after attributes coding with current design of GeS-TM.

- `opti: m66238 - preliminary optimizations of NN search`:  
    This commit optimizes the nearest neighbour search used for motion
    search. It has been used as preliminary improvements in m66238 for
    EE13.60 Test 1.

- `attr: m66238 - dual motion field for color`:  
    - add dual motion field for color
    - add encoder options and hls flags for dual motion
    - update cfg files for CTCs (active with TriSoup)

- `attr: m66238 - motion compensated attributes projection`:  
    - add motion compensated attributes projection from reference frame to
      reconstructed geometry, to replace inter-frame attributes predictor,
    - use depth first approximate NN for motion compensated attributes
      projection,
    - add encoder options and hls flags for motion compensated attributes
      projection and enable it by default when dual motion is active (CTCs),
    - use depth first constructed MSOctree.

- `geom: m66238 - accelerated motion search`:  
    - add depth first NN search for motion search,
    - add depth first approximate NN search to further accelerate encoding,
      at the cost of a decrease in geometry coding performances,
    - add encoder option to enable approximate NN search for motion search.

- `raht: m66286 - RDO for inter mode selection`:  
    This is adopted EE13.60 - Test6, providing rate distorsion optimization
    of inter mode selection for RAHT (also adds signaling of the selected
    mode).

- `raht: fix - use safer approximation of tree depth`:  
    Method used for estimating raht tree depth in m66238 does not
    look safe: the distance between the first an last point in morton
    order is not representative of the actual boundaries of the point
    cloud. Also the function roundLog2 does not look correct.

    This is fixed by only considering the last point in morton order.

- `raht: fix - intra behavior shall not be changed`:  
    Complexity reduction introduced by prediction_threshold1 should be kept
    as it was not discussed to change the behavior of intra frame coding with
    m66286.

- `raht: m66274 - weighted average prediction`:  
    This is adopted EE13.60 - Test 5, after being merged with Test 6. It
    provides a prediction of raht coefficient made by a weighted average
    between inter and intra predictions.

- `raht: fix - out of bound votesum in divisionPredictionVotesLUT`:  
    `votesum` could be 12, when `getNeighborsMode()` is called and the parent
    node has no neighbors.

    This is solved by also considering the mode of the parent node itself
    when computing the values for vote(Inter/Intra[Layer])Weight.

- `misc: add helper functions to FixedPoint`:  
    - Add `+`, `-`, `*` and `/` operators to `FixedPoint`.
    - Add `static FixedPoint::fromVal(int64_t)` to build a `FixedPoint` value
      from its internal `int64_t` representation.

- `raht: m66152 - sample domain prediction`:  
    This is adopted EE13.60 - Test 2, after being merged with Test 5 and
    Test 6. It uses sample domain prediction within a raht layer to avoid
    computing a forward transform in decoding process.

- `raht/enc: m66308 - transform domain distorsion estimation`:  
    Compute distorsion estimation in transform domain to avoid an extra
    inverse transform.

- `raht: m66373 - fix prediction mode determination`:  
    This commit aligns the software with the text description of m62583
    which introduced inter RAHT in GeS-TM. It changes the prediction mode
    determination for RAHT.

### `ges-tm-v4.0`

Tag `ges-tm-v4.0` has been released on top of `ges-tm-v4.0-rc1` since no other
issue had been reported. Gitlab references have been update to reflect the new
address of mpeg gitlab server and the new organization of the projects.

### `ges-tm-v4.0-rc1`

Tag `ges-tm-v4.0-rc1` includes adoptions made during the meeting in
Hannover (October 2023), and some additional cleanups.

- `attr/tidy: m58315 - remove unused quantization variable`:  
    This is same as m58315 adoption, but quantization variables was
    not used at all in GeS-TM since predlift is not part of GeS-TM.

- `attr/raht: m65386 - fix overflow`:  
    Use of int64_t instead of int to avoid potential overflow with high
    bit depth attributes.

- `attr/raht: m65386 - fix attribute QP parameter restriction`:  
    Increases attributes QP range to match G-PCC specification text.

- `attr/raht: m65387 - update encoding QP table`:  
    Increases bitdepth resolution for encoding QP. It has been shown
    on TMC13 that it is far better for high bitdepth attributes, and
    slightly better for 8 bits attributes.

- `attr/raht: m65387 - extend number of coding contexts`:  
    Increase the number of contexts for coding attributes. It has
    been shown on TMC13 that it is much better for high bitdepth
    attributes, and slightly better for 8 bits attributes.

- `attr/raht: m65089 - Test 4.3, intra_mode_level=4`:  
    In RAHT, every frame is divided into binary layers for all-intra
    configuration.

- `attr/tidy: remove unused referencePointCloud`:  
    AttributeInterPredParams::referencePointCloud is not used anymore.
    Removes it.

- `geom/tidy: use const references instead of local copies`:  
    Replaces some local instances/copies by const references in Vec3 and
    PointType member functions.

- `geom/tidy: motion - remove unneeded copies of LPU window`:  
    Removes unnecessary local copies of LPU windows.

- `geom/tidy: motion - optimize memory usage`:  
    Builds indices array of points, before extracting points in well
    dimentionned arrays, to get more optimal use of the memory.

- `geom/tidy: motion - put non compensated nodes in compensatedPointCloud`:  
    Parts of refPointCloud that were not compensated are now copied
    into compensatedPointCloud.
  
    refPointCloud can be removed from trisoup (not necessary anymore).

- `geom/refactor: motion - use PCCPointSet3 instead of vectors`:  
    Uses PCCPointSet3 instead of vectors of coordinates in motion related
    algorithm (search and compensation) to be able to more easily compensate
    with attributes.

- `geom/motion: m64918 - octree based motion search`:  
    Uses an octree based motion search and compensation to remove LPU windows
    and there constrains.
    It allows to use wider motion search and compensation without increasing
    complexity.
  
    The motion window size is also removed from the gps as is no more needed
    for motion compensation in decoder side.
  
    Parameters for octree are modified to keep similar complexity

- `geom/trisoup: m65235 - FaceVertex for Trisoup`:  
    Adds possibility of using vertices on faces in trisoup coding.  
    Also updates CTCs.

- `geom/tidy: m64911 - remove HLS, dead code and sparse algo`:  
    Cleaning GeS-TM, as proposed in m64911.

- `geom/trisoup: m64912 - CTC changes and conditions on centroid activation`:  
    Centroid can now be used with at least 3 points, and it is activated for
    all rate points. CTCs' thikness parameters are also updated.

- `tidy: removal of unused variables and code`:  
    Provides attitional code cleanups:  
    - Removal of unused variables,
    - Removal of code in trisoup:  
      The code had been duplicated in a previous loop iterating on the
      nodes of a slice, for centroid determination with face vertices,
      and the original loopbecame unnecessary.

### `ges-tm-v3.0`

Tag `ges-tm-v3.0` has been released on top of `ges-tm-v3.0-rc1` since no other
issue had been reported.

### `ges-tm-v3.0-rc1`

Tag `ges-tm-v3.0-rc1` includes adoptions made during the meeting in
Geneva (July 2023), a few new features, cleanups and some refactoring.

- `misc: add parenthesis delimiter for arguments parsing`:  
    Arguments/arguments-list can now be delimited by parenthesis to allow
    nested sequencial types.

- `trisoup: m63661 - inter skip mode`:  
    In non-moving parts of a dynamic point cloud, it has been observed that
    the geometry inter prediction based on (zero) motion compensation may
    not be optimal due to lack of invariance of the compensation process.
  
    A new inter coding mode based on colocated vertices and nodes has been
    introduced to mimic a kind of skip mode for point cloud coding.

- `trisoup/ctc: m63660 - align slice BB to zero`:  
    Aligning slice origin to (0, 0, 0) coordinates by default.

- `trisoup: m63660 - hls for skip mode and grid alignment`:  
    Adds high level syntax to align trisoup slices according to trisoup node
    sizes. Also adds restrictions to enforce this alignment with skip mode.

- `refactor: m64005 - refactoring Attribute{Enc,Dec}oder`:  
    This is extracted from monolitic commit made for m64005  in EE repository.
    It refactors Attributes code to get a single templated implementation for
    both 1 dimention attribute luminance, and 3 dimention attributes color.

- `refactor: m64005 - refactoring of RAHT intra`:  
    This is extracted from monolitic commit made for m64005  in EE repository.
    It provides some refactoring to RAHT.

- `restore: divisionless RAHT prediction`:  
    Reverts normative change made during the refactoring of RAHT. Integer
    division was (re-)introduced for inter layer prediction, but had not been
    presented nor discussed with WG7 during the meeting.  
    The original behavior (with tabulated division approximation) has been
    restored.

- `restore: fixed point log2 estimation for RAHT RDOQ`:  
    Reverts non-normative change made during the refactoring of RAHT.
    Fixed point arithmetic for RAHT RDOQ was replaced by floating points
    operations, but had not been presented nor discussed with WG7 during the
    meeting.  
    The original behavior (with fixed point approximation) has been
    restored.

- `restore: original neighbours retrieval for intra RAHT`:  
    Reverts normative change made during the refactoring of RAHT.
    The refactoring of the neighbours retrieval for RAHT was modifying the
    behaviour of the inter layer prediction, but had not been presented nor
    discussed with WG7 during the meeting.
    The original code and behavior has been restored.

- `raht/inter: m64005 - inter RAHT`:  
    This is extracted from monolitic commit made for m64005  in EE repository.
    In provides the inter RAHT discussed with WG7 and adopted.  
  
    Integration notes:  
    - Some modifications affecting all-intra behaviour have been removed,
      as they were not discussed in WG7. It includes the coding of a mode
      for choosing the enabling or not of the inter layer prediction.
    - The syntax element `mode_level` has been moved as an inter-prediction
      dependent syntax element, since it should not affect all-intra.

- `attr/raht: m64218 - inference of prediction modes`:  
    Provides the inference of the prediction modes at certains levels of the
    RAHT decomposition.  
  
    Integration notes:  
    - The syntax element `upper_mode_level` has been moved and made
      dependent on RAHT `enable_inter_prediction`, since the mode
      should be infered only when inter prediction is enabled.

- `attr/raht: m64112 - fix rootLevel value`:  
    This is a fix on the derivation of the rootLevel value within RAHT.
    The rootLevel should be obtained by rounding to the upper integer value
    the division by 3 of the number of RAHT decomposition layers, to get the
    number of dyadic decomposition levels. Because latest one might be
    incomplete.

- `attr/raht: m64118 - integer Haar` 
    Introduces lossless extension to RAHT.  
  
    Integration notes:  
    - Added support for inter-RAHT with m64005:  
      - A simple rounding of the motion compensated predictor has
        been added,
      - The refactoring from m64005 has been modified to re-introduce
        templated use of RAHT/Haar kernels,
      - Kernel support for already normalized weights has been added.
    - Added cfg file for octree-raht-inter with lossless attributes

### `ges-tm-v2.0`

Tag `ges-tm-v2.0` was released on top of `ges-tm-v2.0-rc3` since no other
issue had been reported.

### `ges-tm-v2.0-rc3`

Tag `ges-tm-v2.0-rc3`,
- fixes an out of bound memory access, causing an assertion to occure in debug
  mode. This was unsucessfuly fixed in `ges-tm-v2.0-rc2`.

### `ges-tm-v2.0-rc2`

Tag `ges-tm-v2.0-rc2`,
- fixes uninitialized memory access to an encoding parameters;
- restores the output points being reordered in decoding order within octree
  encoder, which had been accidentaly removed during refactoring.

### `ges-tm-v2.0-rc1`

Tag `ges-tm-v2.0-rc1` includes a few adoptions made during the meeting in
Antalya (April 2023).

- `octree: m62531 - nodes processing in raster scan order`:  
    Use lexicographic order (a.k.a. raster scan order) in octree geometry
    nodes processing to get same node ordering as used in trisoup.

- `entropy: m62547 - probability bounds for dynamic OBUF coders`

- `raht/enc: m63002 - replace log2() with LUT in RDOQ`:  
    Current RDOQ code for RAHT uses the log2() function to estimate the
    rate.

    The proposal simplifies this by using a LUT with a maximum of 16
    entries.

- `trisoup/enc: m62527 - align slices to trisoup node size grid`

- `trisoup: m62526 - Optimal weighted centroid`

- `trisoup/enc: m62981 - centroid quantization offset`

- `trisoup/enc: m62982 - determine more precise centroid`


Adoption of raster-scan ordered nodes processing for octree provides an
unified processing order between octree and trisoup elements.
Tag `ges-tm-v2.0-rc1` includes important refactoring changes on trisoup and
optimizations that were made possible by this unified raster scan ordering,
as well as some code cleanups and simplifications.

The input contribution **m63611: \[GPCC\] Report on new architecture for GeS TM 
v2.0-rc1 and related integrations**
is intended to provide more details on theses modifications.

- `trisoup: m62531 - simplify with raster scan ordered nodes`:  
    simplify trisoup processing by reusing raster scan ordered nodes to
    build edges and vertex informations for trisoup.

- `trisoup: optimize findDominantAxis`

- `trisoup: unit sampling only, simplify accordingly`:  
    GeS-TM is intended to adress solid content.
    Unit sampling shall always be used for triangle projections.

    Encoding is simplified and code is optimized accordingly.

- `trisoup: refactor centroid operations to functions`:  

    Centroid operations are moved to dedicated functions:
    1. determination of centroid and dominant axis
    2. determination of normal vector and bounds
    3. determination of residual/drift
    4. determination of inter predictor for residual/drift

- `trisoup/tidy: cleaning in loop on triangles`

- `trisoup/tidy: better memory management for recPointCloud`:  
    Using preallocation for recPointCloud.

- `trisoup/tidy: better memory management for reconstructed voxels`:  
    Reserve memory for reconstructed voxels of current TriSoup node
    to improve memory management.

- `trisoup/tidy: clean and comment RasterScanTrisoupEdges`

- `trisoup: refactor vertex determination`:  
    Vertex determination is performed in same loop as neighbours and edges
    determination thanks to raster-scan order.

- `trisoup: optimize vertex encoding/decoding`

- `trisoup: refactor vertex inter prediction determination`:  
    Vertex determination for inter prediction is performed in same loop
    as neighbours, edges and intra vertex determination thanks to raster-scan
    order.

    Obsolete determineTrisoupVertices() function is removed.

- `trisoup: refactor vertices encoding/decoding`:  
    Vertex encoding/decoding is performed in same loop as neighbours, edges
    intra and inter vertex determination thanks to raster-scan order.

    Obsolete functions encodeTrisoupVertices() and decodeTrisoupVertices()
    are removed.

- `trisoup/tidy: localize some buffers`:  
    Some buffer are localized inside of function rather than being passed
    as parameters, since they are not used outside anymore.

- `trisoup: refactor triangles rendering`:  
    Put rendering of triangles of a node in specific function.
    Apply it in the same loop as the other trisoup operations,
    to get a single raster scan ordered nodes traversal.

    Remove obsolete function decodeTrisoupCommon().

- `trisoup/tidy: localize some arrays`:  
    Some array are now only used locally for raster scan traversal.

- `trisoup: optimize memory usage`:  
    Use queue instead of vector to get more localized memory, with
    edgePattern, xForedgeOfVertex, and TriSoupVerticesPred.

- `trisoup: optimize rendering/ray tracing`:  
  - use int64 values to represent in 1D 3 dimensional coordinates to
    accelerate duplicate points removal.
  - use one optimized function for each ray tracing directions.
  - optimize/simplify ray tracing.

- `trisoup: optimize duplicates removal and buffer allocation`

- `trisoup/tidy: various cleanups and optimizations`


On trisoup, all the divisions at decoder have been removed and replaced
by fixed point inverse multiplications (normative).

- `trisoup: remove divisions`:  
    remove all divisions at decoder.


Additional cleanups, simplifications and code optimization have been made,
mainly on inter-frame motion search and on octree coding.

- `tidy: minor cleanups`

- `geometry: optimize LPU inter search window generation`:  
  - avoid using unnecessary unordered_map -> decoder twice faster.
  - use BB for points inside the loop.
  - offset the points to slice origin during the loop (avoid multiple passes).

- `geometry: reduce size of LPU inter search window`

- `octree/tidy: remove IDCM variables in node + moved planar container for QTBT up`:  
  !! this may impact the code in case there are duplicated point removing turned on.

- `octree/tidy: optimizations and cleanups`

- `trisoup/tidy: more reasonable fifo size with trisoup`


Other changes includes the disabling of non-cubic nodes that is currently not
supported properly, and a bug fix for and issue that could happen outside of
CTCs.

- `ctc: disable non-cubic nodes`

- `trisoup: fix - determination of centroid predictor`:  
    When inter prediction is not used for a node, the centroid predictor
    could be wrongly estimated.


Documentation files for `ges-tm-v2.0-rc1` are also updated.

- `doc: update documentation files`


### `ges-tm-v1.0`

Tag `ges-tm-v1.0` adds a few fixes.

- `hls/entropy: add backward compatibility flag for m62220`:  
  Applies a fix for backward compatibility issue made in TMC13 for
  entropy coding related contribution m62220.

- `fix: read access outside of an array`:  
  The encoder could read outside of an array.
  Even if the value was not used, it could cause issues when debugging.
  This is fixed by accessing the array only when necessary, and so,
  when the index is inside of the array.

- `raht: fix RDOQ overflows`:  
  In rare cases overflows could occur during computation of Rate/Distortion,
  if distortion was really huge (ex: DC coefficient).
  This is fixed by limiting the computation to small coefficients only.


Additional cleanups were also made, mainly removing unused tools previously
commented out.

- `QTBT is not yet supported, avoid using it`:  
  QTBT is not yet supported by inter frame coding. When inter frame is
  used the tool cannot be activated.
- `removing global motion`
- `removing predictive geometry`
- `removing angular`
- `removing bytewise coder`
- `removing planar`  
  Note: some planar things are kept for QTBT handling.
- `removing predlift`


And the documentation files for `ges-tm-v1.0-rc1` have been added.

- `doc: update documentation files`


### `ges-tm-v1.0-rc1`

Tag `ges-tm-v1.0-rc1` adds all the latest adoptions integrated in the
release 21.0 of TMC13 and related/applicable to dense geometry content,
RAHT attributes and entropy coding.

- `trisoup/m61577: non-cubic nodes`:  
  Add possibility to use non cubic nodes on boundaries of the slices in
  order to avoid reconstruction artefacts and discontinuities between
  neighbouring slices.
  Node cubic nodes are determined based on a bounding box provided for the
  trisoup volume within the slice.

  It can be enabled separately for lowest coordinates in the slice and
  highest ones.

  Encoder decision is also added to determine the bounding box parameters.

- `trisoup/m61561: provide max num reconstructed points to encoder`:  
  Add possibility to provide separately to the encoder a maximum number of
  points reconstructed using trisoup before having to downsample the
  reconstruction.

  This allows to better control the number of input points per slice such
  that better coherence is obtained in the overall reconstructed content.

- `geom/m61583: reduce dynamic OBUF memory footprint`:  
  Reduce the amount of memory being used by dynamic OBUF.

- `raht/m61151: increase buffer precision to remove rounding ops`:  
  The precision is increased in the memory buffer used in RAHT to reconstruct
  attributes at each layer. Now, values are directly stored with fixed
  point precision, avoiding rounding operations during the transform.

- `entropy/m62220: simpler bypass bit with arithmetic coder`:  
  Simplify the operations for arithmetic coding in case of bypassed bits.

- `trisoup/m61982: rasterization with unitary sampling`:  
  Optimization to replace ray tracing by rasterization when no subsampling
  has to be used (trisoup_sampling_value_minus1 is equal to zero), and when
  integer precision rendering shall be used (trisoup_fine_ray_tracing_flag
  is equal to zero).

- `trisoup/m61982: voxelization of vertices on edges`:  
  Replace previous voxelization of trisoup vertices, which was not aligned
  on the edges, by a voxelization aligned on the edges.

  This voxelisation is not needed anymore when no sub-sampling is used:
  the vertices' points are already included in the rendered triangles.

- `trisoup/m61982: use thickness and adaptive halo with rasterization`:  
  Improve trisoup rendering by using thickness and adaptive halo with
  rasterization.


This release includes some extra optimizations.

- `ctc: consider level limit of 5 million points`:  
  For better efficiency, the GeS-TM currently considers hypothetical profile
  level limit of 5 million points per slice so that every frame is coded in
  an single slice.

- `trisoup: optimization of ray tracing`:  
  The ray tracing is performed along a single direction, and reconstructed
  points are restricted to fit inside the node (this could occur because
  of a lack of precision) by clipping the coordinates.

- `trisoup: unique points are selected at block level`:  
  Since points are restricted to be reconstructed locally to a node,
  the selection of the unique points can be (and is) performed at a block
  level rather than at a global point cloud level.

  Sort has ~N*log(N) complexity. Working locally reduces complexity and
  is better for memory access.

Some additionnal code simplifications and optimizations have been made,
as well as some changes in the motion search presets in order to improve
the execution speed:
- `trisoup: optimize determineTriSoupVertices function`,
- `trisoup: code reorganization for faster execution`,
- `trisoup/ctc: change motion search preset`.

Finally `ges-tm-v1.0-rc1` also includes additional code cleanups, minor bug
fixes and modifications to support the common test conditions which have been
agreed during the meeting.

### `ges-tm-v0.1`

Tag `ges-tm-v0.1` corresponds to the software being output in branch
`mpeg140/ee13.60/m61562_ETM_inter_dense` for the 140th MPEG meeting as
experimental model for dynamic dense content coding, including the latest
tools and improvement in `release-v20.0-rc1` adopted in G-PCC for dense
content, and part of earlier exploratory tools studied in EE13.2 for
inter-frame coding in `mpeg-pcc-em13`.


Bug reporting
-------------
Bugs should be reported on the issue tracker set up at
<https://git.mpeg.expert/MPEG/3dgh/g-pcc/software/tm/mpeg-pcc-ges-tm/-/issues>.

