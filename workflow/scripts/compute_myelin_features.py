"""Compute T1w/T2w myelin-proxy features.

Features:
- T1w/T2w ratio at white surface
- T1w/T2w ratio at midthickness surface
- T1w/T2w ratio at pial surface
- T1w/T2w depth profile slope (white-to-pial gradient of myelin ratio)
- T1w/T2w depth profile variance

The T1w/T2w ratio is a well-established proxy for cortical myelin content.
Higher ratios indicate greater myelin density.

References:
- Glasser & Van Essen (2011) Mapping human cortical areas in vivo based on
  myelin content as revealed by T1- and T2-weighted MRI
"""
import logging
import sys

import nibabel as nib
import numpy as np

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)
log = logging.getLogger(__name__)


def save_metric(data, out_path, name):
    """Save a 1D array as a GIFTI shape file."""
    darray = nib.gifti.GiftiDataArray(
        data=data.astype(np.float32),
        intent="NIFTI_INTENT_SHAPE",
        datatype="NIFTI_TYPE_FLOAT32",
        meta=nib.gifti.GiftiMetaData.from_dict({"Name": name}),
    )
    img = nib.gifti.GiftiImage(darrays=[darray])
    nib.save(img, out_path)


try:
    # Load T1w samples
    t1w_white = nib.load(snakemake.input.t1w_white).darrays[0].data.astype(np.float64)
    t1w_pial = nib.load(snakemake.input.t1w_pial).darrays[0].data.astype(np.float64)

    # Load T2w samples
    t2w_white = nib.load(snakemake.input.t2w_white).darrays[0].data.astype(np.float64)
    t2w_pial = nib.load(snakemake.input.t2w_pial).darrays[0].data.astype(np.float64)
    t2w_mid = nib.load(snakemake.input.t2w_mid).darrays[0].data.astype(np.float64)

    # Load midthickness T1w (use depth50 as proxy)
    t1w_mid = nib.load(snakemake.input.t1w_depth50).darrays[0].data.astype(np.float64)
    t2w_depth50 = nib.load(snakemake.input.t2w_depth50).darrays[0].data.astype(np.float64)

    n_vertices = len(t1w_white)
    log.info(f"Computing myelin features for {n_vertices} vertices")

    # Compute T1w/T2w ratio, with protection against division by zero
    def safe_ratio(t1w, t2w, name):
        ratio = np.zeros(n_vertices, dtype=np.float64)
        valid = np.abs(t2w) > 1e-10
        ratio[valid] = t1w[valid] / t2w[valid]
        ratio[~valid] = 0.0
        bad = ~np.isfinite(ratio)
        if np.any(bad):
            log.warning(f"{name}: {np.sum(bad)} non-finite values replaced with 0")
            ratio[bad] = 0.0
        return ratio

    ratio_white = safe_ratio(t1w_white, t2w_white, "ratio_white")
    ratio_pial = safe_ratio(t1w_pial, t2w_pial, "ratio_pial")
    ratio_mid = safe_ratio(t1w_mid, t2w_mid, "ratio_mid")

    log.info(f"Myelin ratio (white): mean={ratio_white.mean():.3f}")
    log.info(f"Myelin ratio (mid): mean={ratio_mid.mean():.3f}")
    log.info(f"Myelin ratio (pial): mean={ratio_pial.mean():.3f}")

    # Slope of myelin ratio from white to pial (3 points: white, mid, pial)
    # Depth positions: 0 (white), 0.5 (mid), 1.0 (pial)
    depths = np.array([0.0, 0.5, 1.0])
    ratio_stack = np.column_stack([ratio_white, ratio_mid, ratio_pial])

    myelin_slope = np.zeros(n_vertices, dtype=np.float64)
    myelin_variance = np.zeros(n_vertices, dtype=np.float64)

    for v in range(n_vertices):
        profile = ratio_stack[v, :]
        if np.all(np.isfinite(profile)):
            coeffs = np.polyfit(depths, profile, 1)
            myelin_slope[v] = coeffs[0]
            myelin_variance[v] = np.var(profile)

    # Clean
    for name, arr in [("myelin_slope", myelin_slope),
                      ("myelin_variance", myelin_variance)]:
        bad = ~np.isfinite(arr)
        if np.any(bad):
            log.warning(f"{name}: {np.sum(bad)} non-finite values replaced with 0")
            arr[bad] = 0.0

    # Save outputs
    save_metric(ratio_white, snakemake.output.ratio_white, "myelin_ratio_white")
    save_metric(ratio_mid, snakemake.output.ratio_mid, "myelin_ratio_midthickness")
    save_metric(ratio_pial, snakemake.output.ratio_pial, "myelin_ratio_pial")
    save_metric(myelin_slope, snakemake.output.ratio_slope, "myelin_slope")
    save_metric(myelin_variance, snakemake.output.ratio_variance, "myelin_variance")

    log.info("Myelin features computed successfully")

except Exception as e:
    log.exception(f"Error computing myelin features: {e}")
    sys.exit(1)
