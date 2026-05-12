#!/bin/bash
#
# Generate a configuration tree in $PWD from YAML files in the same
# directory.

set -e
shopt -s nullglob

script_dir="$(dirname $0)"
geom="trisoup"
attr="raht"
pred="intra"
local="full"
src_cfg_dir=""

while (( $# )); do
	case $1 in
	--octree) geom="octree" ;;
	--trisoup) geom="trisoup" ;;
	--raht) attr="raht" ;;
	--intra) pred="intra" ;;
	--inter) pred="inter" ;;
	--local) local="local" ;;
	--all) all=1 ;;
	--cfgdir=*) src_cfg_dir="${1#--cfgdir=}/" ;;
	--) shift; break ;;
	--help|*)
		echo -e "usage:\n $0\n" \
			"    [[--octree|--trisoup] [--raht] [--intra|--inter] [--local] | --all]\n" \
			"    [--cfgdir=<dir>]"
		exit 1
	esac
	shift;
done

extra_args=("$@")

##
# NB: it is important that the configs in each config set are
# capable of being merged together by gen-cfg.pl.  Ie, no two
# configs may have different definitions of one category.
cfg_octree_raht_full=(
	octree-raht-ctc-lossless-geom-lossy-attrs.yaml
	octree-raht-ctc-lossy-geom-lossy-attrs.yaml
	octree-raht-ctc-lossless-geom-lossless-attrs.yaml
)

cfg_trisoup_raht_full=(
	trisoup-raht-ctc-lossy-geom-lossy-attrs.yaml
)

cfg_octree_raht_local=(
	octree-raht-local-ctc-lossless-geom-lossy-attrs.yaml
	octree-raht-local-ctc-lossy-geom-lossy-attrs.yaml
	octree-raht-local-ctc-lossless-geom-lossless-attrs.yaml
)

cfg_trisoup_raht_local=(
	trisoup-raht-local-ctc-lossy-geom-lossy-attrs.yaml
)

do_one_cfgset() {
	local geom=$1
	local attr=$2
	local pred=$3
	local local=$4

	old_src_cfg_dir=${src_cfg_dir}
	seq_src_cfg_dir=${src_cfg_dir}
	if [[ $pred != "inter" ]]
	then
		outdir="$geom-$attr/"
	else
		outdir="$geom-$attr-inter/"
		src_cfg_dir=${src_cfg_dir}inter/
	fi
	if [[ $local == "local" ]]
	then
		outdir="${outdir%/}-local/"
	fi
	mkdir -p "$outdir"

	cfgset="cfg_${geom}_${attr}_${local}[@]"

	for f in ${!cfgset}
	do
		echo "${src_cfg_dir}$f -> $outdir" ...
	done

	# NB: specifying extra_args at the end does not affect option
	# processing since gen-cfg.pl is flexible in argument positions
	$script_dir/gen-cfg.pl \
		--prefix="$outdir" --no-skip-sequences-without-src \
		"${!cfgset/#/${src_cfg_dir}}" \
		"${seq_src_cfg_dir}sequences-cat2.yaml" \
		"${extra_args[@]}"

	rm -f "$outdir/config-merged.yaml"
	src_cfg_dir=${old_src_cfg_dir}
}

if [[ "$all" != "1" ]]
then
	do_one_cfgset "$geom" "$attr" "$pred" "$local"
else
	do_one_cfgset "octree" "raht" "intra" "full"
	do_one_cfgset "trisoup" "raht" "intra" "full"
	do_one_cfgset "octree" "raht" "inter" "full"
	do_one_cfgset "octree" "raht" "inter" "local"
	do_one_cfgset "trisoup" "raht" "inter" "full"
	do_one_cfgset "trisoup" "raht" "inter" "local"
fi
