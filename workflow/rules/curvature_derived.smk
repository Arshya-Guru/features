"""Rules for computing curvature-derived features.

Computes from principal curvatures (k1, k2) via wb_command:
- Gaussian curvature (k1 * k2)
- Curvedness (sqrt((k1^2 + k2^2) / 2))
- Shape index ((2/pi) * arctan((k1+k2) / (k1-k2)))
"""


rule compute_principal_curvatures:
    """Compute principal curvatures (k1, k2) from a surface using wb_command."""
    input:
        surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{surfname}",
            extension=".surf.gii",
        ),
    output:
        k1=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{surfname}k1",
            extension=".shape.gii",
        ),
        k2=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{surfname}k2",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        surfname="white|pial",
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{surfname}k1k2",
            extension=".log",
        ),
    shell:
        """
        # wb_command -surface-curvature outputs optional -mean, -k1, -k2
        # At least one output must be specified
        wb_command -surface-curvature {input.surf} \
            -k1 {output.k1} \
            -k2 {output.k2} \
            &> {log}
        """


rule compute_curvature_features:
    """Compute Gaussian curvature, curvedness, and shape index from k1/k2."""
    input:
        k1=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="whitek1",
            extension=".shape.gii",
        ),
        k2=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="whitek2",
            extension=".shape.gii",
        ),
    output:
        gaussian=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="gaussiancurvature",
            extension=".shape.gii",
        ),
        curvedness=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="curvedness",
            extension=".shape.gii",
        ),
        shapeindex=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="shapeindex",
            extension=".shape.gii",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="curvaturefeatures",
            extension=".log",
        ),
    script:
        "../scripts/compute_curvature_features.py"
