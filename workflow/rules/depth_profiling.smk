"""Rules for T1w depth profiling along cortical normals.

Computes:
- 9 depth intensity samples (10%-90% of cortical depth)
- 6 profile summary statistics (slope, curvature, variance, skewness,
  superficial/deep ratio, peak depth)
- 3 WM subsurface intensities (1mm, 2mm, 3mm below white surface) + WM slope
- Intensity gradient (pial - white T1w contrast)
"""


rule compute_midthickness:
    """Compute the midthickness surface (average of white and pial coordinates)."""
    input:
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
    output:
        mid=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="midthickness",
            extension=".surf.gii",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="midthickness",
            extension=".log",
        ),
    shell:
        """
        wb_command -surface-average {output.mid} \
            -surf {input.white} \
            -surf {input.pial} \
            &> {log}
        """


rule create_depth_surface:
    """Create an interpolated surface at a given cortical depth fraction.

    Computes coords = (1-frac)*white + frac*pial via wb_command arithmetic.
    depth=10 means 10% from white toward pial, etc.
    """
    input:
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
    output:
        surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="depth{depth}",
            extension=".surf.gii",
        ),
    wildcard_constraints:
        depth="[0-9]+",
    params:
        frac=lambda wc: float(wc.depth) / 100.0,
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="depth{depth}surf",
            extension=".log",
        ),
    script:
        "../scripts/create_depth_surface.py"


rule sample_volume_at_depth:
    """Sample T1w volume at a specific cortical depth fraction.

    Uses trilinear interpolation on the depth-specific surface.
    """
    input:
        volume=bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="brain",
            extension=".nii.gz",
        ),
        surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="depth{depth}",
            extension=".surf.gii",
        ),
    output:
        sampled=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="T1wdepth{depth}",
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
            suffix="T1wdepth{depth}",
            extension=".log",
        ),
    shell:
        """
        wb_command -volume-to-surface-mapping {input.volume} {input.surf} \
            {output.sampled} \
            -trilinear \
            &> {log}
        """


rule sample_volume_at_surface:
    """Sample T1w volume directly at a surface (white or pial)."""
    input:
        volume=bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="brain",
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
            suffix="T1w{surfname}",
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
            suffix="T1w{surfname}",
            extension=".log",
        ),
    shell:
        """
        wb_command -volume-to-surface-mapping {input.volume} {input.surf} \
            {output.sampled} \
            -trilinear \
            &> {log}
        """


rule create_wm_subsurface:
    """Create a surface at a fixed distance below the white matter surface.

    Offsets white surface inward along vertex normals.
    """
    input:
        white=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="white",
            extension=".surf.gii",
        ),
    output:
        surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="wm{depth}mm",
            extension=".surf.gii",
        ),
    wildcard_constraints:
        depth="[1-3]",
    params:
        depth_mm=lambda wc: float(wc.depth),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="wm{depth}mmsurf",
            extension=".log",
        ),
    script:
        "../scripts/create_wm_subsurface.py"


rule sample_wm_subsurface:
    """Sample T1w intensity at a subsurface below white matter."""
    input:
        volume=bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="brain",
            extension=".nii.gz",
        ),
        surf=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="wm{depth}mm",
            extension=".surf.gii",
        ),
    output:
        sampled=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="T1wwm{depth}mm",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        depth="[1-3]",
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="T1wwm{depth}mm",
            extension=".log",
        ),
    shell:
        """
        wb_command -volume-to-surface-mapping {input.volume} {input.surf} \
            {output.sampled} \
            -trilinear \
            &> {log}
        """


rule compute_depth_profiles:
    """Compute depth profile summary statistics from sampled T1w intensities.

    Takes 9 depth samples (10-90%), white and pial surface intensities,
    3 WM subsurface samples, and computes summary features.
    """
    input:
        # 9 depth samples (10% through 90%)
        depth10=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth10", extension=".shape.gii"),
        depth20=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth20", extension=".shape.gii"),
        depth30=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth30", extension=".shape.gii"),
        depth40=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth40", extension=".shape.gii"),
        depth50=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth50", extension=".shape.gii"),
        depth60=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth60", extension=".shape.gii"),
        depth70=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth70", extension=".shape.gii"),
        depth80=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth80", extension=".shape.gii"),
        depth90=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                     suffix="T1wdepth90", extension=".shape.gii"),
        # Surface samples
        t1w_white=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                       suffix="T1wwhite", extension=".shape.gii"),
        t1w_pial=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                      suffix="T1wpial", extension=".shape.gii"),
        # WM subsurface samples
        wm1=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                 suffix="T1wwm1mm", extension=".shape.gii"),
        wm2=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                 suffix="T1wwm2mm", extension=".shape.gii"),
        wm3=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                 suffix="T1wwm3mm", extension=".shape.gii"),
    output:
        slope=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                   suffix="T1wprofileslope", extension=".shape.gii"),
        curvature=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                       suffix="T1wprofilecurvature", extension=".shape.gii"),
        variance=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                      suffix="T1wprofilevariance", extension=".shape.gii"),
        skewness=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                      suffix="T1wprofileskewness", extension=".shape.gii"),
        sd_ratio=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                      suffix="T1wsdratio", extension=".shape.gii"),
        peak_depth=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                        suffix="T1wpeakdepth", extension=".shape.gii"),
        wm_slope=bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
                      suffix="T1wwmslope", extension=".shape.gii"),
        intensity_gradient=bids(root=root, subject="{subject}", hemi="{hemi}",
                                datatype="anat", suffix="T1wintgradient",
                                extension=".shape.gii"),
    log:
        bids(root=root, subject="{subject}", hemi="{hemi}", datatype="anat",
             suffix="depthprofiles", extension=".log"),
    script:
        "../scripts/compute_depth_profiles.py"
