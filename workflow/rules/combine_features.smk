"""Rule for combining all individual feature maps into a single GIFTI file.

Gathers all per-feature .shape.gii outputs and merges them into a single
multi-array GIFTI file for convenient downstream analysis.
"""
import os


def get_all_feature_inputs(wildcards):
    """Collect all individual feature GIFTI files for a subject/hemisphere.

    Returns a dict of named inputs for the combine rule.
    """
    subject = wildcards.subject
    hemi = wildcards.hemi

    def _bids(**kw):
        return bids(root=root, subject=subject, hemi=hemi, datatype="anat", **kw)

    features = {}

    # --- Core morphometrics (3) ---
    for morph in ["thickness", "sulc", "curv"]:
        features[morph] = _bids(suffix=morph, extension=".shape.gii")

    # --- White-gray percent contrast (1) ---
    features["wgpct"] = _bids(suffix="wgpct", extension=".shape.gii")

    # --- Curvature-derived (3) ---
    for feat in ["gaussiancurvature", "curvedness", "shapeindex"]:
        features[feat] = _bids(suffix=feat, extension=".shape.gii")

    # --- Surface geometry (4) ---
    features["pialarea"] = _bids(suffix="pialarea", extension=".shape.gii")
    features["whitearea"] = _bids(suffix="whitearea", extension=".shape.gii")
    features["arearatio"] = _bids(suffix="arearatio", extension=".shape.gii")
    features["normalangle"] = _bids(suffix="normalangle", extension=".shape.gii")

    # --- First-order spatial gradients (8) ---
    for metric in ["thickness", "sulc", "curv", "wgpct",
                   "gaussiancurvature", "curvedness", "shapeindex", "arearatio"]:
        features[f"{metric}gradient"] = _bids(
            suffix=f"{metric}gradient", extension=".shape.gii"
        )

    # --- Second-order gradients (3) ---
    for metric in ["thickness", "sulc", "curv"]:
        features[f"{metric}gradient2"] = _bids(
            suffix=f"{metric}gradient2", extension=".shape.gii"
        )

    # --- T1w depth samples (9) ---
    for d in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
        features[f"T1wdepth{d}"] = _bids(suffix=f"T1wdepth{d}", extension=".shape.gii")

    # --- T1w surface samples (2) ---
    features["T1wwhite"] = _bids(suffix="T1wwhite", extension=".shape.gii")
    features["T1wpial"] = _bids(suffix="T1wpial", extension=".shape.gii")

    # --- WM subsurface intensities (3) ---
    for d in [1, 2, 3]:
        features[f"T1wwm{d}mm"] = _bids(suffix=f"T1wwm{d}mm", extension=".shape.gii")

    # --- Depth profile summaries (8) ---
    for stat in ["T1wprofileslope", "T1wprofilecurvature", "T1wprofilevariance",
                 "T1wprofileskewness", "T1wsdratio", "T1wpeakdepth",
                 "T1wwmslope", "T1wintgradient"]:
        features[stat] = _bids(suffix=stat, extension=".shape.gii")

    # --- Expanded features (4) ---
    for feat in ["lgi", "fractaldim", "corticalcomplexity", "fundaldepth"]:
        features[feat] = _bids(suffix=feat, extension=".shape.gii")

    # --- Myelin features (only if T2 data is available) ---
    # These are optional; the combine script handles missing inputs gracefully
    mri_dir = get_mri_dir(wildcards)
    has_t2 = any(
        os.path.exists(os.path.join(mri_dir, name))
        for name in ["T2.mgz", "T2w.mgz", "T2.nii.gz"]
    )
    if has_t2:
        for feat in ["myelinwhite", "myelinmid", "myelinpial",
                     "myelinslope", "myelinvariance"]:
            features[feat] = _bids(suffix=feat, extension=".shape.gii")

    return features


rule combine_features:
    """Combine all individual feature maps into one multi-array GIFTI."""
    input:
        unpack(get_all_feature_inputs),
    output:
        combined=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="features",
            extension=".shape.gii",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="combinefeatures",
            extension=".log",
        ),
    script:
        "../scripts/combine_features.py"
