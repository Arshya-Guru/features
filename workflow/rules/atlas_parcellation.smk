"""Rule for atlas-based parcellation of vertex-wise features.

Given a FreeSurfer annotation (e.g., Desikan-Killiany, Destrieux, Glasser),
computes regional mean/std/median for each feature and outputs as CSV.
"""


rule atlas_parcellation:
    """Compute regional statistics for each feature within atlas parcels."""
    input:
        features=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="features",
            extension=".shape.gii",
        ),
        annot=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            atlas="{atlas}",
            datatype="anat",
            suffix="dseg",
            extension=".label.gii",
        ),
    output:
        csv=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            atlas="{atlas}",
            datatype="anat",
            suffix="parcellatedfeatures",
            extension=".csv",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            atlas="{atlas}",
            datatype="anat",
            suffix="parcellation",
            extension=".log",
        ),
    script:
        "../scripts/atlas_parcellation.py"
