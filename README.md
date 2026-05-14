# Subjective Quality Towards Simple yet Effective Geometry Processing for Solid G-PCC Application

This repository provides the source code and datasets used in the paper
**Subjective Quality Towards Simple yet Effective Geometry Processing for
Solid G-PCC Application**.

## Repository Structure

The repository is built around the MPEG G-PCC test model and includes two
geometry processing methods evaluated in the paper.

| Path | Description |
| --- | --- |
| `ges-tm-v12rc4/` | The original GPCC/GES-TM v12rc4 source code used as the baseline. |
| `ges-tm-v12rc4-HighFidelityCoordinateScale/` | The implementation of **Method A** in the paper, named High Fidelity Coordinate Scale. |
| `ges-tm-v12rc4-FourLoopDecoupledTriangleReconstruction/` | The implementation of **Method B** in the paper, named Four-Loop Decoupled Triangle Reconstruction for Non-Closed Nodes. |

## Datasets

The point cloud that records large-scale and high-precision electric power system in the wild is available at:

- `https://pan.baidu.com/s/1001yds9Rs4fcOdP3Uqb3xA?pwd=km6m`

The experiments also use MPEG point cloud datasets. The MPEG dataset website
is added here:

- MPEG point cloud datasets: `https://content.mpeg.expert/data/MPEG-I/`

## Building

Each code project contains its own README file with build and usage
instructions:

- Baseline GPCC/GES-TM v12rc4:
  [`ges-tm-v12rc4/README.md`](ges-tm-v12rc4/README.md)
- Method A, High Fidelity Coordinate Scale:
  [`ges-tm-v12rc4-HighFidelityCoordinateScale/README.md`](ges-tm-v12rc4-HighFidelityCoordinateScale/README.md)
- Method B, Four-Loop Decoupled Triangle Reconstruction:
  [`ges-tm-v12rc4-FourLoopDecoupledTriangleReconstruction/README.md`](ges-tm-v12rc4-FourLoopDecoupledTriangleReconstruction/README.md)

All three projects are based on CMake. A typical build workflow is:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Please refer to the README file in each project directory for platform-specific
commands, release build options, and usage examples.

## Notes

- `ges-tm-v12rc4/` should be treated as the reference GPCC/GES-TM v12rc4 code.
- `HighFidelityCoordinateScale` and `FourLoopDecoupledTriangleReconstruction`
  are modified versions of the baseline code for the two methods proposed in
  the paper.