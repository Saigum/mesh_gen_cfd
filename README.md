# Coronary Mesh & Simulation Setup Pipeline

This repository provides a pipeline to segment coronary vessels from medical imaging data, prepare 3D surface meshes, and generate an OpenFOAM simulation setup. It uses Python scripts along with OpenFOAM utilities to automate the workflow.

---

## 📁 Repository Structure

```
├── environment_.yml            # Conda environment specification
├── gen.py                      # Main pipeline script
├── generate_foam_setup.py      # Module to write OpenFOAM dictionaries and run mesh generation
└── README.md                   # Project instructions and usage
```

---

## 🛠 Prerequisites

* **Operating System**: Linux (Ubuntu recommended)
* **Python**: 3.8 or later
* **Conda**: Anaconda or Miniconda to manage environments
* **OpenFOAM**: Version 2412 installed and accessible (e.g., via `/usr/lib/openfoam/openfoam2412/etc/bashrc`)
* **VTK / VMTK**: For mesh and skeleton processing
* **TotalSegmentator**: Python package for automated medical image segmentation

---

## 🔧 Installation and Environment Setup

1. **Clone the repository**:

   ```bash
   git clone <your-repo-url>
   cd <your-repo-folder>
   ```

2. **Create and activate the Conda environment**:

   ```bash
   conda env create -f environment_.yml
   conda activate coronary-pipeline
   ```

3. **Source OpenFOAM environment**:

   ```bash
   source /usr/lib/openfoam/openfoam2412/etc/bashrc
   ```

> **Tip:** You can add the above line to your `~/.bashrc` or `~/.zshrc` to load OpenFOAM automatically.

---

## 🚀 Usage

All processing is handled by the `gen.py` script, which invokes segmentation, mesh extraction, visualization, and OpenFOAM setup.

### Command-line Arguments

* `--label_nii_path` (required): Path to the labeled NIfTI file containing coronary segmentation (e.g., `.nii.gz`).
* `--img_nii_path` (required): Path to the original medical image NIfTI file for aorta segmentation.
* `--stl_output_path` (optional, default: `coronary_surface.stl`): Filename for the extracted coronary STL.
* `--prefix` (required): Output directory prefix for storing pipeline results and OpenFOAM case.

### Example

```bash
python gen.py \
  --label_nii_path data/coronary_label.nii.gz \
  --img_nii_path data/patient_scan.nii.gz \
  --stl_output_path coronary_surface.stl \
  --prefix output/run1
```

This will:

1. Segment the coronary vessels (`gen.py` → `output_stl`).
2. Identify vessel terminal points and skeleton.
3. Segment the aorta and intersect masks.
4. Visualize combined segmentation in HTML (`output/run1/coronary_aorta_segmentation.html`).
5. Copy STL to `output/run1/mesh/constant/triSurface/`.
6. Generate OpenFOAM dictionaries and run `blockMesh`, `surfaceFeatureExtract`, `snappyHexMesh`, and `foamToVTK`.
7. Clone mesh directory into a working OpenFOAM case at `output/run1/foam`.
8. Save terminal points in `output/run1/foam/constant/terminal_points.txt`.

---

## ⚙️ Module: `generate_foam_setup.py`

This module exports:

* `mesh_generation(...)`:

  * Writes OpenFOAM dictionary files under `<prefix>/mesh/system/`:

    * `blockMeshDict`, `createPatchDict`, `controlDict`, `fvSchemes`, `fvSolution`, `meshQualityDict`, `snappyHexMeshDict`, `surfaceFeatureExtractDict`
  * Invokes OpenFOAM commands (`blockMesh`, `surfaceFeatureExtract`, `snappyHexMesh -overwrite`, `foamToVTK`).

You normally do not call this script directly; it is invoked by `gen.py` with the appropriate arguments.

---

## 🐞 Troubleshooting

* **File not found errors**:

  * Check that your NIfTI paths (`--label_nii_path`, `--img_nii_path`) are correct.
* **OpenFOAM utilities missing**:

  * Ensure OpenFOAM 2412 is installed and that you have sourced its `bashrc`.
* **TotalSegmentator failures**:

  * Verify that the `totalsegmentator` package is installed in your Conda environment.

---


