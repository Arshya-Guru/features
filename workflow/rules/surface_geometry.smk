"""Rules for computing surface geometry features.

Computes:
- Pial surface vertex area
- White surface vertex area
- Area ratio (pial/white) — proxy for cortical folding complexity
- Surface normal angle between pial and white surfaces
"""


rule compute_vertex_area:
    """Compute per-vertex area on a surface using wb_command."""
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
        area=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{surfname}area",
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
            suffix="{surfname}area",
            extension=".log",
        ),
    shell:
        "wb_command -surface-vertex-areas {input.surf} {output.area} &> {log}"


rule compute_surface_geometry:
    """Compute area ratio and normal angle from pial and white surfaces."""
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
    output:
        area_ratio=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="arearatio",
            extension=".shape.gii",
        ),
        normal_angle=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="normalangle",
            extension=".shape.gii",
        ),
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="surfgeom",
            extension=".log",
        ),
    script:
        "../scripts/compute_surface_geometry.py"
