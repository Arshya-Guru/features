"""Compute curvature-derived features from principal curvatures (k1, k2).

Features:
- Gaussian curvature: k1 * k2 (positive = dome/cup, negative = saddle)
- Curvedness: sqrt((k1^2 + k2^2) / 2) (overall bending magnitude)
- Shape index: (2/pi) * arctan((k1+k2)/(k1-k2)) (normalized shape descriptor [-1,1])

References:
- Koenderink & van Doorn (1992) Surface shape and curvature scales
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

try:
    # Load principal curvatures
    k1_img = nib.load(snakemake.input.k1)
    k2_img = nib.load(snakemake.input.k2)

    k1 = k1_img.darrays[0].data.astype(np.float64)
    k2 = k2_img.darrays[0].data.astype(np.float64)
    n_vertices = len(k1)
    log.info(f"Loaded principal curvatures for {n_vertices} vertices")

    # Gaussian curvature = k1 * k2
    gaussian = k1 * k2

    # Curvedness = sqrt((k1^2 + k2^2) / 2)
    curvedness = np.sqrt((k1**2 + k2**2) / 2.0)

    # Shape index = (2/pi) * arctan((k1+k2) / (k1-k2))
    # Handle division by zero where k1 == k2 (umbilical points)
    denom = k1 - k2
    shape_index = np.zeros_like(k1)
    valid = np.abs(denom) > 1e-10
    shape_index[valid] = (2.0 / np.pi) * np.arctan2(k1[valid] + k2[valid], denom[valid])

    # Replace NaN/Inf with 0
    for arr_name, arr in [("gaussian", gaussian), ("curvedness", curvedness),
                          ("shape_index", shape_index)]:
        bad = ~np.isfinite(arr)
        n_bad = np.sum(bad)
        if n_bad > 0:
            log.warning(f"{arr_name}: {n_bad} non-finite vertices replaced with 0")
            arr[bad] = 0.0

    # Save as GIFTI shape files
    def save_metric(data, out_path, name):
        darray = nib.gifti.GiftiDataArray(
            data=data.astype(np.float32),
            intent="NIFTI_INTENT_SHAPE",
            datatype="NIFTI_TYPE_FLOAT32",
            meta=nib.gifti.GiftiMetaData.from_dict({"Name": name}),
        )
        img = nib.gifti.GiftiImage(darrays=[darray])
        nib.save(img, out_path)
        log.info(f"Saved {name} to {out_path}")

    save_metric(gaussian, snakemake.output.gaussian, "gaussian_curvature")
    save_metric(curvedness, snakemake.output.curvedness, "curvedness")
    save_metric(shape_index, snakemake.output.shapeindex, "shape_index")

    log.info("Curvature features computed successfully")

except Exception as e:
    log.exception(f"Error computing curvature features: {e}")
    sys.exit(1)
