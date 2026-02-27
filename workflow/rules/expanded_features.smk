"""Rules for expanded cortical features.

Computes:
- Local gyrification index (LGI) approximated from area ratio
- Fractal dimensionality / cortical complexity from pial surface
- Sulcal fundal depth (local minima of sulcal depth)
- Sulcal depth gradient direction (vector magnitude of gradient components)
"""


rule compute_expanded_features:
    """Compute expanded cortical features: LGI, fractal dimension, sulcal morphometry."""
    input:
        pial_surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="pial",
            extension=".surf.gii",
        ),
        white_surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="white",
            extension=".surf.gii",
        ),
        pial_area=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="pialarea",
            extension=".shape.gii",
        ),
        white_area=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="whitearea",
            extension=".shape.gii",
        ),
        sulc=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="sulc",
            extension=".shape.gii",
        ),
    output:
        lgi=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="lgi",
            extension=".shape.gii",
        ),
        fractal_dim=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="fractaldim",
            extension=".shape.gii",
        ),
        cortical_complexity=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="corticalcomplexity",
            extension=".shape.gii",
        ),
        fundal_depth=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="fundaldepth",
            extension=".shape.gii",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="expandedfeatures",
            extension=".log",
        ),
    script:
        "../scripts/compute_expanded_features.py"
