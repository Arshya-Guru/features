"""Compute depth profile summary statistics from sampled T1w intensities.

Takes 9 cortical depth samples (10-90%), white and pial surface intensities,
and 3 WM subsurface samples to compute:

Profile summaries (6):
- slope: linear slope of intensity across cortical depth
- curvature: quadratic term of polynomial fit to depth profile
- variance: variance of intensity across depths
- skewness: skewness of intensity distribution across depths
- superficial/deep ratio: mean(60-90%) / mean(10-40%)
- peak depth: depth fraction with maximum intensity

WM features (2):
- WM slope: linear intensity gradient in white matter (1-3mm below surface)
- Intensity gradient: pial T1w - white T1w contrast
"""
import logging
import sys

import nibabel as nib
import numpy as np
from scipy import stats as sp_stats

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
    # Load depth samples into a (n_vertices x 9) array
    depth_names = [
        "depth10", "depth20", "depth30", "depth40", "depth50",
        "depth60", "depth70", "depth80", "depth90",
    ]
    depth_fracs = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

    depth_data = []
    for dname in depth_names:
        img = nib.load(getattr(snakemake.input, dname))
        depth_data.append(img.darrays[0].data.astype(np.float64))
    depth_array = np.column_stack(depth_data)  # (n_vertices, 9)

    n_vertices = depth_array.shape[0]
    log.info(f"Loaded depth profiles for {n_vertices} vertices")

    # Load surface intensities
    t1w_white = nib.load(snakemake.input.t1w_white).darrays[0].data.astype(np.float64)
    t1w_pial = nib.load(snakemake.input.t1w_pial).darrays[0].data.astype(np.float64)

    # Load WM subsurface samples
    wm1 = nib.load(snakemake.input.wm1).darrays[0].data.astype(np.float64)
    wm2 = nib.load(snakemake.input.wm2).darrays[0].data.astype(np.float64)
    wm3 = nib.load(snakemake.input.wm3).darrays[0].data.astype(np.float64)

    # --- Profile slope: linear fit across depths ---
    # Use polyfit degree 1 for each vertex
    slope = np.zeros(n_vertices, dtype=np.float64)
    curvature_prof = np.zeros(n_vertices, dtype=np.float64)

    for v in range(n_vertices):
        profile = depth_array[v, :]
        if np.all(np.isfinite(profile)):
            # Quadratic fit: a*x^2 + b*x + c
            coeffs = np.polyfit(depth_fracs, profile, 2)
            curvature_prof[v] = coeffs[0]  # quadratic term
            slope[v] = coeffs[1]  # linear term
        else:
            slope[v] = 0.0
            curvature_prof[v] = 0.0

    log.info("Computed slope and curvature")

    # --- Profile variance ---
    variance = np.nanvar(depth_array, axis=1)
    variance[~np.isfinite(variance)] = 0.0

    # --- Profile skewness ---
    # Use scipy skew with nan_policy
    with np.errstate(all="ignore"):
        skewness = sp_stats.skew(depth_array, axis=1, nan_policy="omit")
    skewness[~np.isfinite(skewness)] = 0.0

    # --- Superficial/deep ratio ---
    # Superficial = mean of depths 60-90%, deep = mean of depths 10-40%
    deep_mean = np.nanmean(depth_array[:, :4], axis=1)      # 10-40%
    superficial_mean = np.nanmean(depth_array[:, 5:], axis=1)  # 60-90%
    sd_ratio = np.ones(n_vertices, dtype=np.float64)
    valid = np.abs(deep_mean) > 1e-10
    sd_ratio[valid] = superficial_mean[valid] / deep_mean[valid]
    sd_ratio[~np.isfinite(sd_ratio)] = 1.0

    # --- Peak depth: depth fraction with maximum intensity ---
    peak_idx = np.nanargmax(depth_array, axis=1)
    peak_depth = depth_fracs[peak_idx]

    # --- WM slope: linear intensity gradient in WM (1-3mm below surface) ---
    wm_stack = np.column_stack([wm1, wm2, wm3])  # (n_vertices, 3)
    wm_depths = np.array([1.0, 2.0, 3.0])
    wm_slope = np.zeros(n_vertices, dtype=np.float64)
    for v in range(n_vertices):
        wm_profile = wm_stack[v, :]
        if np.all(np.isfinite(wm_profile)):
            coeffs = np.polyfit(wm_depths, wm_profile, 1)
            wm_slope[v] = coeffs[0]

    # --- Intensity gradient: pial - white T1w ---
    intensity_gradient = t1w_pial - t1w_white
    intensity_gradient[~np.isfinite(intensity_gradient)] = 0.0

    # Clean all outputs
    for name, arr in [("slope", slope), ("curvature", curvature_prof),
                      ("variance", variance), ("skewness", skewness),
                      ("sd_ratio", sd_ratio), ("peak_depth", peak_depth),
                      ("wm_slope", wm_slope), ("intensity_gradient", intensity_gradient)]:
        bad = ~np.isfinite(arr)
        if np.any(bad):
            log.warning(f"{name}: {np.sum(bad)} non-finite values replaced with 0")
            arr[bad] = 0.0

    # Save all outputs
    save_metric(slope, snakemake.output.slope, "profile_slope")
    save_metric(curvature_prof, snakemake.output.curvature, "profile_curvature")
    save_metric(variance, snakemake.output.variance, "profile_variance")
    save_metric(skewness, snakemake.output.skewness, "profile_skewness")
    save_metric(sd_ratio, snakemake.output.sd_ratio, "superficial_deep_ratio")
    save_metric(peak_depth, snakemake.output.peak_depth, "peak_depth")
    save_metric(wm_slope, snakemake.output.wm_slope, "wm_slope")
    save_metric(intensity_gradient, snakemake.output.intensity_gradient, "intensity_gradient")

    log.info("Depth profile features computed successfully")

except Exception as e:
    log.exception(f"Error computing depth profiles: {e}")
    sys.exit(1)
