# EpLink Cortical Feature Extraction Pipeline

A [Snakebids](https://snakebids.readthedocs.io) workflow that extracts 34 vertex-wise cortical surface features from [DeepPrep](https://github.com/pBFSLab/DeepPrep) (FreeSurfer-style) derivatives of T1-weighted MRI data.

## Feature Groups (34 total)

### 1. Curvature-Derived (3)
Computed from principal curvatures (k₁, k₂) via `wb_command -surface-curvature`:
- **Gaussian curvature** — k₁ × k₂ (dome vs saddle vs flat)
- **Curvedness** — √((k₁² + k₂²)/2) (how "bent" the surface is)
- **Shape index** — (2/π) arctan((k₁+k₂)/(k₁−k₂)) (normalized shape descriptor)

### 2. Surface Geometry (4)
Derived from comparing pial and white surface meshes:
- **Pial surface area** — vertex-wise area on pial surface
- **White surface area** — vertex-wise area on white surface
- **Area ratio** — pial/white (cortical folding complexity)
- **Surface normal angle** — arccos(dot(pial_normal, white_normal))

### 3. Spatial Gradients (4)
`wb_command -metric-gradient` applied to base morphometric features:
- **Gradient of thickness, sulcal depth, curvature, intensity gradient**

### 4. T1w Depth Profiling (23)
Sampled from `brain.mgz` along white→pial normals:
- **9 depth intensity samples** (10%–90% cortical depth)
- **6 profile summary statistics** (slope, curvature, variance, skewness, superficial/deep ratio, peak depth)
- **4 WM subsurface features** (intensity at 1–3mm below white + slope)
- **Intensity gradient** (pial − white T1w contrast)

## Requirements

- Python ≥ 3.10
- [Snakemake](https://snakemake.github.io) ≥ 8.0
- [Snakebids](https://snakebids.readthedocs.io) ≥ 0.14
- [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) (for `mris_convert`, `mri_convert`)
- [Connectome Workbench](https://www.humanconnectome.org/software/workbench-command) (for `wb_command`)
- Python packages: `nibabel`, `numpy`, `scipy`

## Expected DeepPrep Directory Structure

The pipeline expects FreeSurfer-style outputs from DeepPrep:

```
deepprep_derivatives/
├── sub-001/
│   ├── surf/
│   │   ├── lh.white          # White matter surface
│   │   ├── rh.white
│   │   ├── lh.pial           # Pial surface
│   │   ├── rh.pial
│   │   ├── lh.thickness      # Cortical thickness
│   │   ├── lh.sulc           # Sulcal depth
│   │   ├── lh.curv           # Mean curvature
│   │   └── lh.w-g.pct.mgh   # White-gray percent contrast
│   └── mri/
│       └── brain.mgz         # Bias-corrected T1w volume
├── sub-002/
│   └── ...
```

> **Note**: If your DeepPrep output uses different filenames (e.g., the gradient
> overlay is named differently), edit the `convert_gradient` rule in
> `workflow/rules/convert_surfaces.smk`.

## Usage

### As a BIDS App

```bash
python run.py /path/to/bids_dataset /path/to/output participant \
    --deepprep_dir /path/to/derivatives/deepprep \
    --cores 8
```

### Workflow Mode (for development)

1. Edit `config/snakebids.yml` with your paths:
   ```yaml
   bids_dir: '/path/to/bids'
   output_dir: '/path/to/output'
   deepprep_dir: '/path/to/derivatives/deepprep'
   ```

2. First run (generates config):
   ```bash
   python run.py /path/to/bids . participant --deepprep_dir /path/to/deepprep
   ```

3. Subsequent runs (direct Snakemake):
   ```bash
   snakemake --snakefile workflow/Snakefile --configfile config/snakebids.yml --cores 8
   ```

### Single Subject

```bash
python run.py /path/to/bids /path/to/output participant \
    --participant-label 001 \
    --deepprep_dir /path/to/derivatives/deepprep \
    --cores 4
```

## Output Structure

```
output/
├── sub-001/
│   └── anat/
│       ├── sub-001_hemi-L_features.shape.gii     # Combined 34-feature GIFTI
│       ├── sub-001_hemi-R_features.shape.gii
│       ├── sub-001_hemi-L_white.surf.gii          # Converted surfaces
│       ├── sub-001_hemi-L_pial.surf.gii
│       ├── sub-001_hemi-L_thickness.shape.gii     # Individual feature maps
│       ├── sub-001_hemi-L_gaussiancurvature.shape.gii
│       └── ...
```

The combined `features.shape.gii` file contains all 34 features as named data arrays in a single GIFTI, ready for downstream analysis or atlas-based parcellation.

## Pipeline DAG

```
T1w (BIDS) ──> enumerate subjects
                     │
DeepPrep outputs ────┤
    │                │
    ├─ surfaces ─────┼── convert to GIFTI ──┬── curvature-derived (3)
    │   (white,pial) │                      ├── surface geometry (4)
    │                │                      ├── spatial gradients (4)
    ├─ morphometry ──┼── convert to GIFTI ──┘
    │   (thick,sulc, │
    │    curv,grad)  │
    │                │
    └─ brain.mgz ────┼── convert to NIfTI ──── depth profiling (23)
                     │
                     └── combine all ──> features.shape.gii
```
# features
