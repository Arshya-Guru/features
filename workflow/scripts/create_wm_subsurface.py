"""Create a surface offset below the white matter surface along vertex normals.

Computes inward-pointing normals from the white surface and offsets each
vertex by a specified distance (in mm) into the white matter.
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


def compute_vertex_normals(coords, faces):
    """Compute outward-pointing vertex normals from mesh."""
    v0 = coords[faces[:, 0]]
    v1 = coords[faces[:, 1]]
    v2 = coords[faces[:, 2]]

    face_normals = np.cross(v1 - v0, v2 - v0)

    vertex_normals = np.zeros_like(coords)
    for i in range(3):
        np.add.at(vertex_normals, faces[:, i], face_normals)

    norms = np.linalg.norm(vertex_normals, axis=1, keepdims=True)
    norms[norms < 1e-10] = 1.0
    vertex_normals /= norms

    return vertex_normals


try:
    depth_mm = snakemake.params.depth_mm
    log.info(f"Creating WM subsurface at {depth_mm}mm below white surface")

    white_img = nib.load(snakemake.input.white)
    white_coords = white_img.darrays[0].data.astype(np.float64)
    faces = white_img.darrays[1].data

    n_vertices = len(white_coords)

    # Compute outward normals, then offset inward (subtract)
    normals = compute_vertex_normals(white_coords, faces)
    offset_coords = (white_coords - depth_mm * normals).astype(np.float32)

    # Create GIFTI surface
    coord_array = nib.gifti.GiftiDataArray(
        data=offset_coords,
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

    log.info(f"Saved WM subsurface ({depth_mm}mm) to {snakemake.output.surf}")

except Exception as e:
    log.exception(f"Error creating WM subsurface: {e}")
    sys.exit(1)
