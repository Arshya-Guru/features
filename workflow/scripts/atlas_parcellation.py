"""Compute atlas-based regional statistics for each cortical feature.

Given a FreeSurfer annotation (converted to GIFTI label format) and the
combined feature GIFTI, computes mean/std/median for each feature within
each atlas parcel. Outputs a CSV file with rows as regions and columns
as feature statistics.
"""
import logging
import sys

import nibabel as nib
import numpy as np
import csv

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
)
log = logging.getLogger(__name__)

try:
    # Load combined features
    features_img = nib.load(snakemake.input.features)
    n_features = len(features_img.darrays)
    n_vertices = len(features_img.darrays[0].data)
    log.info(f"Loaded {n_features} features for {n_vertices} vertices")

    # Get feature names from metadata
    feature_names = []
    feature_data = []
    for darray in features_img.darrays:
        name = darray.meta.get("Name", f"feature_{len(feature_names)}")
        feature_names.append(name)
        feature_data.append(darray.data.astype(np.float64))

    feature_data = np.column_stack(feature_data)  # (n_vertices, n_features)

    # Load atlas annotation (GIFTI label format)
    annot_img = nib.load(snakemake.input.annot)
    labels = annot_img.darrays[0].data.astype(np.int32)

    # Get label table (mapping from label index to region name)
    label_table = annot_img.labeltable
    label_dict = {}
    if label_table is not None and label_table.labels:
        for label in label_table.labels:
            label_dict[label.key] = label.label
    else:
        # Fall back to unique label values
        for lbl in np.unique(labels):
            label_dict[int(lbl)] = f"region_{lbl}"

    log.info(f"Atlas has {len(label_dict)} regions")

    # Compute regional statistics
    rows = []
    header = ["region_index", "region_name", "n_vertices"]
    for fname in feature_names:
        header.extend([f"{fname}_mean", f"{fname}_std", f"{fname}_median"])

    unique_labels = sorted(np.unique(labels))
    for lbl in unique_labels:
        if lbl == 0:
            # Skip the "unknown" / "medial wall" label (index 0)
            continue

        region_name = label_dict.get(lbl, f"region_{lbl}")
        mask = labels == lbl
        n_region = int(np.sum(mask))

        if n_region == 0:
            continue

        row = [lbl, region_name, n_region]
        for fi in range(n_features):
            vals = feature_data[mask, fi]
            valid = np.isfinite(vals)
            if np.sum(valid) > 0:
                v = vals[valid]
                row.extend([float(np.mean(v)), float(np.std(v)), float(np.median(v))])
            else:
                row.extend([0.0, 0.0, 0.0])

        rows.append(row)

    # Write CSV
    with open(snakemake.output.csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

    log.info(f"Wrote parcellated features to {snakemake.output.csv} "
             f"({len(rows)} regions x {n_features} features)")

except Exception as e:
    log.exception(f"Error in atlas parcellation: {e}")
    sys.exit(1)
