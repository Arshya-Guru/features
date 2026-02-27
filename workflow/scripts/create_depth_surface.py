"""Create an interpolated surface at a given cortical depth fraction.

Interpolates between white and pial surfaces:
    coords = (1 - frac) * white_coords + frac * pial_coords

where frac = depth_percentage / 100 (0 = white, 1 = pial).
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
    frac = snakemake.params.frac
    log.info(f"Creating depth surface at fraction {frac}")

    # Load white and pial surfaces
    white_img = nib.load(snakemake.input.white)
    pial_img = nib.load(snakemake.input.pial)

    white_coords = white_img.darrays[0].data.astype(np.float64)
    pial_coords = pial_img.darrays[0].data.astype(np.float64)
    faces = white_img.darrays[1].data  # topology is the same

    n_vertices = len(white_coords)
    log.info(f"Interpolating {n_vertices} vertices at depth {frac:.2f}")

    # Linear interpolation
    depth_coords = ((1.0 - frac) * white_coords + frac * pial_coords).astype(np.float32)

    # Create GIFTI surface
    coord_array = nib.gifti.GiftiDataArray(
        data=depth_coords,
        intent="NIFTI_INTENT_POINTSET",
        datatype="NIFTI_TYPE_FLOAT32",
    )
    face_array = nib.gifti.GiftiDataArray(
        data=faces,
        intent="NIFTI_INTENT_TRIANGLE",
        datatype="NIFTI_TYPE_INT32",
    )
    img = nib.gifti.GiftiImage(darrays=[coord_array, face_array])
    nib.save(img, snakemake.output.surf)

    log.info(f"Saved depth surface to {snakemake.output.surf}")

except Exception as e:
    log.exception(f"Error creating depth surface: {e}")
    sys.exit(1)
