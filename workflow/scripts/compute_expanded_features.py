"""Compute expanded cortical features.

Features:
- Local gyrification index (LGI): smoothed area ratio at multiple scales,
  approximated by computing pial/white area ratio with neighborhood averaging
- Fractal dimensionality: local fractal dimension estimated from surface mesh
  at multiple scales using the box-counting approach on local patches
- Cortical complexity: surface roughness metric from pial surface
- Fundal depth: identifies sulcal fundi as local minima of sulcal depth

References:
- Schaer et al. (2012) local gyrification index
- Yotter et al. (2011) fractal dimension of cortical surfaces
"""
import logging
import sys

import nibabel as nib
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

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


def build_adjacency(faces, n_vertices):
    """Build sparse adjacency matrix from triangle faces."""
    rows = np.concatenate([faces[:, 0], faces[:, 1], faces[:, 2],
                           faces[:, 1], faces[:, 2], faces[:, 0]])
    cols = np.concatenate([faces[:, 1], faces[:, 2], faces[:, 0],
                           faces[:, 0], faces[:, 1], faces[:, 2]])
    data = np.ones(len(rows), dtype=np.float32)
    adj = csr_matrix((data, (rows, cols)), shape=(n_vertices, n_vertices))
    return adj


def get_k_ring_neighbors(adj, vertex, k):
    """Get all vertices within k edges of the given vertex."""
    visited = {vertex}
    frontier = {vertex}
    for _ in range(k):
        new_frontier = set()
        for v in frontier:
            neighbors = adj[v].indices
            for n in neighbors:
                if n not in visited:
                    visited.add(n)
                    new_frontier.add(n)
        frontier = new_frontier
    return np.array(list(visited))


def compute_local_area_ratio_smoothed(pial_area, white_area, adj, n_vertices, k=5):
    """Approximate LGI as smoothed pial/white area ratio over k-ring neighborhoods.

    For each vertex, sum pial area and white area in the k-ring neighborhood
    and compute the ratio. This approximates the local gyrification index
    without requiring the full outer hull computation.
    """
    lgi = np.ones(n_vertices, dtype=np.float64)

    # Build k-ring smoothing by repeated matrix multiplication
    smooth_pial = pial_area.copy()
    smooth_white = white_area.copy()
    adj_binary = (adj > 0).astype(np.float64)

    for _ in range(k):
        # Add neighbor values (including self)
        smooth_pial = adj_binary.dot(smooth_pial) + smooth_pial
        smooth_white = adj_binary.dot(smooth_white) + smooth_white

    valid = smooth_white > 1e-10
    lgi[valid] = smooth_pial[valid] / smooth_white[valid]
    lgi[~valid] = 1.0

    return lgi


def compute_fractal_dimension(coords, faces, adj, n_vertices, scales=(3, 5, 8)):
    """Estimate local fractal dimension from surface mesh at multiple scales.

    Uses a simplified approach: for each vertex, measure the surface area
    within neighborhoods of increasing radius (k-ring) and estimate the
    fractal dimension from the log-log scaling of area vs radius.
    """
    fractal_dim = np.full(n_vertices, 2.0, dtype=np.float64)

    # Precompute face areas
    v0 = coords[faces[:, 0]]
    v1 = coords[faces[:, 1]]
    v2 = coords[faces[:, 2]]
    face_areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)

    # Map faces to vertices (which faces touch each vertex)
    vertex_faces = [[] for _ in range(n_vertices)]
    for fi, face in enumerate(faces):
        for vi in face:
            vertex_faces[vi].append(fi)

    # For computational efficiency, sample a subset of vertices
    # and interpolate for the rest
    adj_binary = (adj > 0).astype(np.float64)

    log_areas = []
    log_radii = np.log(np.array(scales, dtype=np.float64))

    for scale in scales:
        # Compute neighborhood area at this scale by repeated smoothing
        # Accumulate face areas at each vertex
        vertex_area = np.zeros(n_vertices, dtype=np.float64)
        for vi in range(n_vertices):
            vertex_area[vi] = sum(face_areas[fi] / 3.0 for fi in vertex_faces[vi])

        # Smooth to k-ring neighborhood
        smoothed = vertex_area.copy()
        for _ in range(scale):
            smoothed = adj_binary.dot(smoothed) + smoothed

        log_areas.append(np.log(np.maximum(smoothed, 1e-20)))

    # Fit log(area) vs log(radius) for each vertex
    log_areas = np.column_stack(log_areas)  # (n_vertices, n_scales)
    for v in range(n_vertices):
        y = log_areas[v, :]
        if np.all(np.isfinite(y)):
            coeffs = np.polyfit(log_radii, y, 1)
            fractal_dim[v] = coeffs[0]  # slope = fractal dimension estimate

    return fractal_dim


def compute_cortical_complexity(pial_coords, pial_faces, adj, n_vertices, k=3):
    """Compute cortical complexity as local surface roughness.

    Measures how much the local pial surface deviates from a smooth sphere,
    computed as the ratio of actual surface area to the area of the convex
    hull of the local neighborhood.

    Simplified: variance of vertex-to-centroid distances in local patches.
    """
    complexity = np.zeros(n_vertices, dtype=np.float64)
    adj_binary = (adj > 0).astype(np.float64)

    # For each vertex, compute the variance of distances to neighbors
    for v in range(n_vertices):
        neighbors = adj[v].indices
        if len(neighbors) < 3:
            continue

        # Get k-ring neighborhood
        visited = {v}
        frontier = set(neighbors)
        visited.update(frontier)
        for _ in range(k - 1):
            new_frontier = set()
            for n in frontier:
                for nn in adj[n].indices:
                    if nn not in visited:
                        visited.add(nn)
                        new_frontier.add(nn)
            frontier = new_frontier

        nbr_idx = np.array(list(visited))
        local_coords = pial_coords[nbr_idx]

        # Centroid
        centroid = local_coords.mean(axis=0)
        dists = np.linalg.norm(local_coords - centroid, axis=1)

        # Complexity = coefficient of variation of distances
        if dists.mean() > 1e-10:
            complexity[v] = dists.std() / dists.mean()

    return complexity


def compute_fundal_depth(sulc, adj, n_vertices):
    """Identify sulcal fundi as local minima of sulcal depth.

    A vertex is a fundal point if its sulcal depth is greater than
    all its immediate neighbors (sulc is inverted: more positive = deeper sulcus).
    Returns the fundal depth value for fundal vertices, 0 elsewhere.
    """
    fundal = np.zeros(n_vertices, dtype=np.float64)

    for v in range(n_vertices):
        neighbors = adj[v].indices
        if len(neighbors) == 0:
            continue
        # In FreeSurfer convention, more positive sulc = deeper sulcus
        v_depth = sulc[v]
        neighbor_depths = sulc[neighbors]

        # Check if this vertex is a local maximum of sulcal depth
        if v_depth > np.max(neighbor_depths):
            fundal[v] = v_depth

    return fundal


try:
    # Load inputs
    pial_img = nib.load(snakemake.input.pial_surf)
    white_img = nib.load(snakemake.input.white_surf)

    pial_coords = pial_img.darrays[0].data.astype(np.float64)
    pial_faces = pial_img.darrays[1].data
    white_coords = white_img.darrays[0].data.astype(np.float64)

    pial_area = nib.load(snakemake.input.pial_area).darrays[0].data.astype(np.float64)
    white_area = nib.load(snakemake.input.white_area).darrays[0].data.astype(np.float64)
    sulc = nib.load(snakemake.input.sulc).darrays[0].data.astype(np.float64)

    n_vertices = len(pial_coords)
    log.info(f"Loaded data for {n_vertices} vertices")

    # Build adjacency matrix from pial surface
    adj = build_adjacency(pial_faces, n_vertices)
    log.info("Built adjacency matrix")

    # Local gyrification index
    lgi = compute_local_area_ratio_smoothed(pial_area, white_area, adj, n_vertices, k=5)
    log.info(f"LGI: mean={lgi.mean():.3f}, std={lgi.std():.3f}")

    # Fractal dimension
    fractal_dim = compute_fractal_dimension(pial_coords, pial_faces, adj, n_vertices)
    log.info(f"Fractal dim: mean={fractal_dim.mean():.3f}, std={fractal_dim.std():.3f}")

    # Cortical complexity
    complexity = compute_cortical_complexity(pial_coords, pial_faces, adj, n_vertices, k=2)
    log.info(f"Complexity: mean={complexity.mean():.4f}, std={complexity.std():.4f}")

    # Fundal depth
    fundal = compute_fundal_depth(sulc, adj, n_vertices)
    n_fundal = np.sum(fundal > 0)
    log.info(f"Fundal depth: {n_fundal} fundal vertices identified")

    # Clean non-finite values
    for name, arr in [("lgi", lgi), ("fractal_dim", fractal_dim),
                      ("complexity", complexity), ("fundal", fundal)]:
        bad = ~np.isfinite(arr)
        if np.any(bad):
            log.warning(f"{name}: {np.sum(bad)} non-finite values replaced with 0")
            arr[bad] = 0.0

    # Save outputs
    save_metric(lgi, snakemake.output.lgi, "local_gyrification_index")
    save_metric(fractal_dim, snakemake.output.fractal_dim, "fractal_dimension")
    save_metric(complexity, snakemake.output.cortical_complexity, "cortical_complexity")
    save_metric(fundal, snakemake.output.fundal_depth, "fundal_depth")

    log.info("Expanded features computed successfully")

except Exception as e:
    log.exception(f"Error computing expanded features: {e}")
    sys.exit(1)
