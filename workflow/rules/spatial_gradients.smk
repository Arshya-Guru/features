"""Rules for computing spatial gradients of metric maps.

Computes surface gradients using wb_command -metric-gradient:
- Gradient of thickness
- Gradient of sulcal depth
- Gradient of curvature
- Gradient of intensity gradient (w-g percent contrast)
- Gradients of curvature-derived features (Gaussian, curvedness, shape index)
- Gradients of area ratio
- Second-order gradients (gradient of gradient) for key features
"""


rule compute_metric_gradient:
    """Compute the spatial gradient of a metric on the surface."""
    input:
        metric=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{metricname}",
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
        gradient=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{metricname}gradient",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        # All metrics for which we compute first-order gradients
        metricname="thickness|sulc|curv|wgpct|gaussiancurvature|curvedness|shapeindex|arearatio",
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{metricname}gradient",
            extension=".log",
        ),
    shell:
        """
        wb_command -metric-gradient {input.surf} {input.metric} {output.gradient} \
            &> {log}
        """


rule compute_second_order_gradient:
    """Compute second-order gradient (gradient of gradient) for key features."""
    input:
        metric=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{metricname}gradient",
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
        gradient2=bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{metricname}gradient2",
            extension=".shape.gii",
        ),
    wildcard_constraints:
        # Second-order gradients for the most informative features
        metricname="thickness|sulc|curv",
    log:
        bids(
            root=root,
            subject="{subject}",
            hemi="{hemi}",
            datatype="anat",
            suffix="{metricname}gradient2",
            extension=".log",
        ),
    shell:
        """
        wb_command -metric-gradient {input.surf} {input.metric} {output.gradient2} \
            &> {log}
        """
