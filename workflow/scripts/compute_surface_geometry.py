"""Compute surface geometry features from pial and white surfaces.

Features:
- Area ratio: pial_area / white_area (folding complexity proxy)
- Normal angle: angle between pial and white surface normals at each vertex
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
    """Compute per-vertex normals from mesh coordinates and faces.

    Uses area-weighted face normals averaged at each vertex.
    """
    v0 = coords[faces[:, 0]]
    v1 = coords[faces[:, 1]]
    v2 = coords[faces[:, 2]]

    # Face normals via cross product
    face_normals = np.cross(v1 - v0, v2 - v0)

    # Accumulate face normals at each vertex
    vertex_normals = np.zeros_like(coords)
    for i in range(3):
        np.add.at(vertex_normals, faces[:, i], face_normals)

    # Normalize
    norms = np.linalg.norm(vertex_normals, axis=1, keepdims=True)
    norms[norms < 1e-10] = 1.0  # avoid division by zero
    vertex_normals /= norms

    return vertex_normals


try:
    # Load surfaces
    pial_img = nib.load(snakemake.input.pial_surf)
    white_img = nib.load(snakemake.input.white_surf)

    pial_coords = pial_img.darrays[0].data
    pial_faces = pial_img.darrays[1].data
    white_coords = white_img.darrays[0].data
    white_faces = white_img.darrays[1].data

    n_vertices = len(pial_coords)
    log.info(f"Loaded surfaces with {n_vertices} vertices")

    # Load pre-computed vertex areas
    pial_area_img = nib.load(snakemake.input.pial_area)
    white_area_img = nib.load(snakemake.input.white_area)
    pial_area = pial_area_img.darrays[0].data.astype(np.float64)
    white_area = white_area_img.darrays[0].data.astype(np.float64)

    # Area ratio (pial / white) — proxy for cortical folding complexity
    area_ratio = np.ones(n_vertices, dtype=np.float64)
    valid = white_area > 1e-10
    area_ratio[valid] = pial_area[valid] / white_area[valid]
    area_ratio[~valid] = 1.0
    log.info(f"Area ratio: mean={area_ratio.mean():.3f}, std={area_ratio.std():.3f}")

    # Surface normal angle — angle between pial and white normals
    pial_normals = compute_vertex_normals(pial_coords, pial_faces)
    white_normals = compute_vertex_normals(white_coords, white_faces)

    # Dot product, clipped to [-1, 1] for numerical stability
    dot = np.sum(pial_normals * white_normals, axis=1)
    dot = np.clip(dot, -1.0, 1.0)
    normal_angle = np.arccos(dot)  # in radians
    log.info(f"Normal angle: mean={np.degrees(normal_angle.mean()):.2f} deg")

    # Replace NaN/Inf
    for name, arr in [("area_ratio", area_ratio), ("normal_angle", normal_angle)]:
        bad = ~np.isfinite(arr)
        n_bad = np.sum(bad)
        if n_bad > 0:
            log.warning(f"{name}: {n_bad} non-finite vertices replaced with 0")
            arr[bad] = 0.0

    # Save outputs
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

    save_metric(area_ratio, snakemake.output.area_ratio, "area_ratio")
    save_metric(normal_angle, snakemake.output.normal_angle, "normal_angle")

    log.info("Surface geometry features computed successfully")

except Exception as e:
    log.exception(f"Error computing surface geometry: {e}")
    sys.exit(1)
