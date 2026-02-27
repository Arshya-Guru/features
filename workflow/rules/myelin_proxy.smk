"""Rules for T1w/T2w myelin-proxy features (optional, requires T2w data).

If T2.mgz is available in the subject's mri/ directory, computes:
- T1w/T2w ratio at white, midthickness, and pial surfaces
- T1w/T2w ratio at depth fractions (10%-90%)
- T1w/T2w depth profile slope and variance
"""
import os


def t2_available(wildcards):
    """Check if T2 volume exists for this subject."""
    mri_dir = get_mri_dir(wildcards)
    for t2_name in ["T2.mgz", "T2w.mgz", "T2.nii.gz"]:
        if os.path.exists(os.path.join(mri_dir, t2_name)):
            return True
    return False


rule sample_t2_at_surface:
    """Sample T2w volume at a specific surface (white, pial, or midthickness)."""
    input:
        volume=bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="T2brain",
            extension=".nii.gz",
        ),
        surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{surfname}",
            extension=".surf.gii",
        ),
    output:
        sampled=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="T2w{surfname}",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        surfname="white|pial|midthickness",
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="T2w{surfname}",
            extension=".log",
        ),
    shell:
        """
        wb_command -volume-to-surface-mapping {input.volume} {input.surf} \
            {output.sampled} \
            -trilinear \
            &> {log}
        """


rule sample_t2_at_depth:
    """Sample T2w volume at cortical depth fractions using ribbon-constrained mapping."""
    input:
        volume=bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="T2brain",
            extension=".nii.gz",
        ),
        white=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="white",
            extension=".surf.gii",
        ),
        pial=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="pial",
            extension=".surf.gii",
        ),
        mid=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="midthickness",
            extension=".surf.gii",
        ),
    output:
        sampled=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="T2wdepth{depth}",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        depth="[0-9]+",
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="T2wdepth{depth}",
            extension=".log",
        ),
    shell:
        """
        wb_command -volume-to-surface-mapping {input.volume} {input.mid} \
            {output.sampled} \
            -ribbon-constrained {input.white} {input.pial} \
            &> {log}
        """


rule compute_myelin_features:
    """Compute T1w/T2w myelin-proxy features at multiple depths and surfaces."""
    input:
        # T1w surface samples
        t1w_white=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                       suffix="T1wwhite", extension=".shape.gii"),
        t1w_pial=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                      suffix="T1wpial", extension=".shape.gii"),
        # T2w surface samples
        t2w_white=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                       suffix="T2wwhite", extension=".shape.gii"),
        t2w_pial=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                      suffix="T2wpial", extension=".shape.gii"),
        t2w_mid=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T2wmidthickness", extension=".shape.gii"),
        # T1w and T2w depth samples for profile analysis
        t1w_depth50=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                         suffix="T1wdepth50", extension=".shape.gii"),
        t2w_depth50=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                         suffix="T2wdepth50", extension=".shape.gii"),
    output:
        ratio_white=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                         suffix="myelinwhite", extension=".shape.gii"),
        ratio_mid=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                       suffix="myelinmid", extension=".shape.gii"),
        ratio_pial=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                        suffix="myelinpial", extension=".shape.gii"),
        ratio_slope=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                         suffix="myelinslope", extension=".shape.gii"),
        ratio_variance=bids(root=root, subject="{subject}", hemi="{hemi}",
                            datatype="anat", suffix="myelinvariance",
                            extension=".shape.gii"),
    log:
        bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
             suffix="myelinfeatures", extension=".log"),
    script:
        "../scripts/compute_myelin_features.py"
