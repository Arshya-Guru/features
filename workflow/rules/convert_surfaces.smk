"""Rules for converting FreeSurfer surfaces and overlays to GIFTI format.

Converts:
- White and pial surfaces (.white, .pial -> .surf.gii)
- Morphometric overlays (.thickness, .sulc, .curv -> .shape.gii)
- White-gray intensity contrast (.w-g.pct.mgh -> .shape.gii)
- brain.mgz -> brain.nii.gz (for depth profiling)
"""


# --- Surface conversion (FreeSurfer -> GIFTI) ---
rule convert_surface:
    """Convert a FreeSurfer surface file to GIFTI using mris_convert."""
    input:
        surf=lambda wc: f"{get_surf_dir(wc)}/{hemi_map[wc.hemi]}.{wc.surfname}",
    output:
        gii=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{surfname}",
            extension=".surf.gii",
        ),
    wildcard_constraints:
        surfname="white|pial",
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{surfname}",
            extension=".convert_surface.log",
        ),
    shell:
        "mris_convert {input.surf} {output.gii} &> {log}"


# --- Morphometric overlay conversion (FreeSurfer curv format -> GIFTI) ---
rule convert_morph:
    """Convert a FreeSurfer morphometric overlay to GIFTI shape.

    Uses wb_command -metric-convert after mris_convert to get shape.gii format.
    """
    input:
        overlay=lambda wc: f"{get_surf_dir(wc)}/{hemi_map[wc.hemi]}.{wc.morphname}",
        surf=lambda wc: f"{get_surf_dir(wc)}/{hemi_map[wc.hemi]}.white",
    output:
        gii=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{morphname}",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        morphname="thickness|sulc|curv",
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{morphname}",
            extension=".convert_morph.log",
        ),
    shell:
        """
        # Convert FreeSurfer curv -> GIFTI via mris_convert
        mris_convert -c {input.overlay} {input.surf} {output.gii} &> {log}
        # FreeSurfer 8 prepends hemisphere prefix to output filename - fix it
        _hemi=$(basename {input.surf} | cut -d. -f1)
        _dir=$(dirname {output.gii})
        _base=$(basename {output.gii})
        if [ ! -f "{output.gii}" ] && [ -f "$_dir/$_hemi.$_base" ]; then
            mv "$_dir/$_hemi.$_base" "{output.gii}" >> {log} 2>&1
        fi
        """


# --- White-gray intensity gradient conversion ---
rule convert_gradient:
    """Convert white-gray percent contrast overlay (.w-g.pct.mgh) to GIFTI."""
    input:
        mgh=lambda wc: f"{get_surf_dir(wc)}/{hemi_map[wc.hemi]}.w-g.pct.mgh",
        surf=lambda wc: f"{get_surf_dir(wc)}/{hemi_map[wc.hemi]}.white",
    output:
        gii=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="wgpct",
            extension=".shape.gii",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="wgpct",
            extension=".convert_gradient.log",
        ),
    shell:
        """
        # w-g.pct.mgh is a surface overlay in MGH format
        mris_convert -c {input.mgh} {input.surf} {output.gii} &> {log}
        # FreeSurfer 8 prepends hemisphere prefix to output filename - fix it
        _hemi=$(basename {input.surf} | cut -d. -f1)
        _dir=$(dirname {output.gii})
        _base=$(basename {output.gii})
        if [ ! -f "{output.gii}" ] && [ -f "$_dir/$_hemi.$_base" ]; then
            mv "$_dir/$_hemi.$_base" "{output.gii}" >> {log} 2>&1
        fi
        """


# --- Volume conversion (brain.mgz -> NIfTI for depth profiling) ---
rule convert_brain_volume:
    """Convert brain.mgz to NIfTI for volume sampling along surface normals."""
    input:
        mgz=lambda wc: f"{get_mri_dir(wc)}/brain.mgz",
    output:
        nii=bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="brain",
            extension=".nii.gz",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="brain",
            extension=".convert_volume.log",
        ),
    shell:
        "mri_convert {input.mgz} {output.nii} &> {log}"


# --- Optional T2 volume conversion ---
def get_t2_path(wildcards):
    """Find the T2 volume path for this subject. Returns first match found."""
    import os
    mri_dir = get_mri_dir(wildcards)
    for t2_name in ["T2.mgz", "T2w.mgz", "T2.nii.gz"]:
        path = os.path.join(mri_dir, t2_name)
        if os.path.exists(path):
            return path
    # Return a path that snakemake will report as missing (but this rule
    # should only be triggered when T2 exists, guarded by combine_features)
    return os.path.join(mri_dir, "T2.mgz")


rule convert_t2_volume:
    """Convert T2.mgz to NIfTI (if available) for myelin-proxy features."""
    input:
        mgz=lambda wc: get_t2_path(wc),
    output:
        nii=bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="T2brain",
            extension=".nii.gz",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            datatype="anat",
            suffix="T2brain",
            extension=".convert_volume.log",
        ),
    shell:
        "mri_convert {input.mgz} {output.nii} &> {log}"


# --- FreeSurfer annotation conversion (for atlas parcellation) ---
rule convert_annotation:
    """Convert a FreeSurfer annotation file to GIFTI label format."""
    input:
        annot=lambda wc: f"{get_annot_dir(wc)}/{hemi_map[wc.hemi]}.{wc.atlas}.annot",
        surf=lambda wc: f"{get_surf_dir(wc)}/{hemi_map[wc.hemi]}.white",
    output:
        gii=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            atlas="{atlas}",
            datatype="anat",
            suffix="dseg",
            extension=".label.gii",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            atlas="{atlas}",
            datatype="anat",
            suffix="dseg",
            extension=".convert_annot.log",
        ),
    shell:
        """
        mris_convert --annot {input.annot} {input.surf} {output.gii} &> {log}
        # FreeSurfer 8 prepends hemisphere prefix to output filename - fix it
        _hemi=$(basename {input.surf} | cut -d. -f1)
        _dir=$(dirname {output.gii})
        _base=$(basename {output.gii})
        if [ ! -f "{output.gii}" ] && [ -f "$_dir/$_hemi.$_base" ]; then
            mv "$_dir/$_hemi.$_base" "{output.gii}" >> {log} 2>&1
        fi
        """
