"""Compute principal curvatures (k1, k2) from a surface.

Uses wb_command -surface-curvature to get mean (H) and Gaussian (K) curvature,
then derives principal curvatures:
    k1 = H + sqrt(H^2 - K)
    k2 = H - sqrt(H^2 - K)

This works with wb_command 1.4.x which supports -mean and -gauss but not -k1/-k2.
"""
import logging
import subprocess
import sys
import tempfile

import nibabel as nib
import numpy as np

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)
log = logging.getLogger(__name__)

try:
    surf_path = snakemake.input.surf
    k1_path = snakemake.output.k1
    k2_path = snakemake.output.k2

    # Compute mean and Gaussian curvature via wb_command
    with tempfile.NamedTemporaryFile(suffix=".shape.gii") as mean_f, \
         tempfile.NamedTemporaryFile(suffix=".shape.gii") as gauss_f:

        mean_path = mean_f.name
        gauss_path = gauss_f.name

        cmd = (
            f"wb_command -surface-curvature {surf_path}"
            f" -mean {mean_path}"
            f" -gauss {gauss_path}"
        )
        log.info(f"Running: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            log.error(f"wb_command failed: {result.stderr}")
            sys.exit(1)

        # Load mean and Gaussian curvature
        H = nib.load(mean_path).darrays[0].data.astype(np.float64)
        K = nib.load(gauss_path).darrays[0].data.astype(np.float64)

    log.info(f"Loaded curvatures for {len(H)} vertices")

    # Derive principal curvatures from H and K
    # k1 = H + sqrt(H^2 - K), k2 = H - sqrt(H^2 - K)
    discriminant = H**2 - K
    # Clamp negative discriminant (numerical noise) to zero
    discriminant = np.maximum(discriminant, 0.0)
    sqrt_disc = np.sqrt(discriminant)

    k1 = H + sqrt_disc
    k2 = H - sqrt_disc

    # Replace any non-finite values
    for name, arr in [("k1", k1), ("k2", k2)]:
        bad = ~np.isfinite(arr)
        n_bad = np.sum(bad)
        if n_bad > 0:
            log.warning(f"{name}: {n_bad} non-finite vertices replaced with 0")
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

    save_metric(k1, k1_path, "k1")
    save_metric(k2, k2_path, "k2")
    log.info("Principal curvatures computed successfully")

except Exception as e:
    log.exception(f"Error computing principal curvatures: {e}")
    sys.exit(1)
