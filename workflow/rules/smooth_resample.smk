"""Rules for smoothing and resampling cortical features to standard meshes.

Smooths each feature on the white surface (Gaussian kernel), then resamples
to standard mesh densities (fsaverage, fsaverage5, fsaverage6) using
barycentric interpolation via wb_command.
"""


# --- Feature list (all 48 features from the combined GIFTI) ---
SMOOTH_FEATURES = [
    # Core morphometrics
    "thickness", "sulc", "curv", "wgpct",
    # Curvature-derived
    "gaussiancurvature", "curvedness", "shapeindex",
    # Surface geometry
    "pialarea", "whitearea", "arearatio", "normalangle",
    # First-order gradients
    "thicknessgradient", "sulcgradient", "curvgradient", "wgpctgradient",
    "gaussiancurvaturegradient", "curvednessgradient", "shapeindexgradient",
    "arearatiogradient",
    # Second-order gradients
    "thicknessgradient2", "sulcgradient2", "curvgradient2",
    # T1w depth samples
    "T1wdepth10", "T1wdepth20", "T1wdepth30", "T1wdepth40", "T1wdepth50",
    "T1wdepth60", "T1wdepth70", "T1wdepth80", "T1wdepth90",
    # T1w surface samples
    "T1wwhite", "T1wpial",
    # WM subsurface
    "T1wwm1mm", "T1wwm2mm", "T1wwm3mm",
    # Depth profile summaries
    "T1wprofileslope", "T1wprofilecurvature", "T1wprofilevariance",
    "T1wprofileskewness", "T1wsdratio", "T1wpeakdepth",
    "T1wwmslope", "T1wintgradient",
    # Expanded features
    "lgi", "fractaldim", "corticalcomplexity", "fundaldepth",
]

# Wildcard constraint regex for feature names
_smooth_feat_re = "|".join(SMOOTH_FEATURES)

# Config values
smooth_sigma = int(config.get("smooth_sigma", 10))
smooth_desc = f"smooth{smooth_sigma}mm"
template_dir = config.get("template_dir", None)
resolutions = config.get("resolutions", ["fsaverage", "fsaverage5", "fsaverage6"])
_resolutions_re = "|".join(resolutions)


# --- Smoothing ---
rule smooth_feature:
    """Smooth a metric map on the white surface using a Gaussian kernel."""
    input:
        metric=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{featname}",
            extension=".shape.gii",
        ),
        surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="white",
            extension=".surf.gii",
        ),
    output:
        smoothed=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            desc=smooth_desc,
            suffix="{featname}",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        featname=_smooth_feat_re,
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            desc=smooth_desc,
            suffix="{featname}",
            extension=".smooth.log",
        ),
    shell:
        """
        wb_command -metric-smoothing {input.surf} {input.metric} \
            {smooth_sigma} {output.smoothed} &> {log}
        """


# --- Resampling ---
rule resample_feature:
    """Resample a smoothed metric to a standard mesh resolution."""
    input:
        smoothed=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            desc=smooth_desc,
            suffix="{featname}",
            extension=".shape.gii",
        ),
        native_sphere=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="sphere",
            extension=".surf.gii",
        ),
        template_sphere=lambda wc: (
            f"{template_dir}/{wc.resolution}/surf/"
            f"{hemi_map[wc.hemi]}.sphere.surf.gii"
        ),
    output:
        resampled=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            space="{resolution}",
            datatype="anat",
            desc=smooth_desc,
            suffix="{featname}",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        featname=_smooth_feat_re,
        resolution=_resolutions_re,
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            space="{resolution}",
            datatype="anat",
            desc=smooth_desc,
            suffix="{featname}",
            extension=".resample.log",
        ),
    shell:
        """
        wb_command -metric-resample {input.smoothed} \
            {input.native_sphere} {input.template_sphere} \
            BARYCENTRIC {output.resampled} &> {log}
        """


# --- Combine resampled features ---
def get_all_resampled_inputs(wildcards):
    """Collect all resampled feature files for a subject/hemi/resolution."""
    features = {}
    for feat in SMOOTH_FEATURES:
        features[feat] = bids(
            root=root,
            subject=wildcards.subject,
            hemi=wildcards.hemi,
            space=wildcards.resolution,
            datatype="anat",
            desc=smooth_desc,
            suffix=feat,
            extension=".shape.gii",
        )
    return features


rule combine_resampled:
    """Combine all resampled features into one multi-array GIFTI per resolution."""
    input:
        unpack(get_all_resampled_inputs),
    output:
        combined=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            space="{resolution}",
            datatype="anat",
            desc=smooth_desc,
            suffix="features",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        resolution=_resolutions_re,
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            space="{resolution}",
            datatype="anat",
            desc=smooth_desc,
            suffix="combineresampled",
            extension=".log",
        ),
    script:
        "../scripts/combine_features.py"
