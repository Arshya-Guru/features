# EpLink Cortical Feature Extraction Pipeline

A [Snakebids](https://snakebids.readthedocs.io) workflow that extracts vertex-wise cortical surface features from [DeepPrep](https://github.com/pBFSLab/DeepPrep) (FreeSurfer-style) derivatives of T1-weighted MRI data.

Supports multiple FreeSurfer-style derivative layouts: **DeepPrep**, **FreeSurfer**, **FastSurfer**, and **custom** path templates.

## Feature Groups (48+ features)

### 1. Core Morphometrics (4)
Direct from FreeSurfer outputs, converted to GIFTI:
- **Cortical thickness** — distance between white and pial surfaces
- **Sulcal depth** — depth of cortical sulci
- **Mean curvature** — average surface curvature
- **White-gray percent contrast** — T1w intensity contrast at the gray/white boundary

### 2. Curvature-Derived (3)
Computed from principal curvatures (k1, k2) via `wb_command -surface-curvature`:
- **Gaussian curvature** — k1 x k2
- **Curvedness** — sqrt((k1^2 + k2^2)/2)
- **Shape index** — (2/pi) arctan((k1+k2)/(k1-k2))

### 3. Surface Geometry (4)
Derived from comparing pial and white surface meshes:
- **Pial surface area** — vertex-wise area on pial surface
- **White surface area** — vertex-wise area on white surface
- **Area ratio** — pial/white (cortical folding complexity)
- **Surface normal angle** — arccos(dot(pial_normal, white_normal))

### 4. Expanded Features (4)
- **Local gyrification index (LGI)** — smoothed area ratio approximation
- **Fractal dimension** — local fractal dimension from multi-scale surface analysis
- **Cortical complexity** — surface roughness from local distance variance
- **Fundal depth** — sulcal depth at local maxima (fundal vertices)

### 5. Spatial Gradients (11)
`wb_command -metric-gradient` applied to morphometric features:
- **First-order gradients (8):** thickness, sulc, curv, wgpct, Gaussian curvature, curvedness, shape index, area ratio
- **Second-order gradients (3):** gradient of gradient for thickness, sulc, curv

### 6. T1w Depth Profiling (22)
Sampled from `brain.mgz` along white-to-pial normals:
- **9 depth intensity samples** (10%-90% cortical depth)
- **2 surface intensities** (white, pial)
- **3 WM subsurface intensities** (1-3mm below white surface)
- **6 profile summaries** (slope, curvature, variance, skewness, superficial/deep ratio, peak depth)
- **WM slope** (intensity gradient in white matter)
- **Intensity gradient** (pial - white T1w contrast)

### 7. Myelin-Proxy Features (5, optional)
Computed if T2w data is available:
- **T1w/T2w ratio** at white, midthickness, pial surfaces
- **Myelin depth profile slope**
- **Myelin depth profile variance**

### 8. Atlas-Based Parcellation (optional)
Given FreeSurfer annotation files (e.g., Desikan-Killiany, Destrieux):
- Regional **mean/std/median** for each feature
- Output as CSV per subject/hemisphere/atlas

## Requirements

- Python >= 3.10
- [Snakemake](https://snakemake.github.io) >= 8.0
- [Snakebids](https://snakebids.readthedocs.io) >= 0.14
- [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) (for `mris_convert`, `mri_convert`)
- [Connectome Workbench](https://www.humanconnectome.org/software/workbench-command) (for `wb_command`)
- Python packages: `nibabel`, `numpy`, `scipy`

## Supported Derivative Layouts

| Layout | Path pattern |
|--------|-------------|
| `deepprep` | `{dir}/sub-{subject}/Recon/sub-{subject}/surf/` |
| `freesurfer` | `{dir}/sub-{subject}/surf/` |
| `fastsurfer` | `{dir}/sub-{subject}/surf/` |
| `custom` | User-specified template |

## Usage

### As a BIDS App

```bash
python run.py /path/to/bids_dataset /path/to/output participant \
    --deepprep_dir /path/to/derivatives/deepprep \
    --skip-bids-validation \
    --cores 8
```

### With a different layout

```bash
python run.py /path/to/bids /path/to/output participant \
    --deepprep_dir /path/to/freesurfer_output \
    --freesurfer_layout freesurfer \
    --cores 4
```

### With atlas parcellation

```bash
python run.py /path/to/bids /path/to/output participant \
    --deepprep_dir /path/to/derivatives \
    --atlas aparc aparc.a2009s \
    --skip-bids-validation \
    --cores 4
```

### Single subject

```bash
python run.py /path/to/bids /path/to/output participant \
    --participant-label 001 \
    --deepprep_dir /path/to/derivatives/deepprep \
    --skip-bids-validation \
    --cores 4
```

## Output Structure

```
output/
├── sub-001/
│   └── anat/
│       ├── sub-001_hemi-L_features.shape.gii     # Combined multi-feature GIFTI
│       ├── sub-001_hemi-R_features.shape.gii
│       ├── sub-001_hemi-L_white.surf.gii          # Converted surfaces
│       ├── sub-001_hemi-L_pial.surf.gii
│       ├── sub-001_hemi-L_thickness.shape.gii     # Individual feature maps
│       ├── sub-001_hemi-L_gaussiancurvature.shape.gii
│       ├── sub-001_hemi-L_lgi.shape.gii
│       ├── sub-001_hemi-L_fractaldim.shape.gii
│       ├── sub-001_hemi-L_atlas-aparc_parcellatedfeatures.csv
│       └── ...
```
