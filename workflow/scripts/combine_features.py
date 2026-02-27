"""Combine all individual feature maps into a single multi-array GIFTI file.

Reads each named input from snakemake.input, extracts the first data array,
and combines them into a single GIFTI with named arrays.
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
    darrays = []
    feature_names = []

    # snakemake.input is a named dict-like object
    # Iterate over all named inputs
    for name in sorted(snakemake.input.keys()):
        path = snakemake.input[name]
        if not isinstance(path, str):
            # skip non-string entries
            continue

        try:
            img = nib.load(path)
            data = img.darrays[0].data.astype(np.float32)

            # Replace NaN/Inf with 0
            bad = ~np.isfinite(data)
            if np.any(bad):
                log.warning(f"{name}: {np.sum(bad)} non-finite vertices replaced with 0")
                data[bad] = 0.0

            darray = nib.gifti.GiftiDataArray(
                data=data,
                intent="NIFTI_INTENT_SHAPE",
                datatype="NIFTI_TYPE_FLOAT32",
                meta=nib.gifti.GiftiMetaData.from_dict({"Name": name}),
            )
            darrays.append(darray)
            feature_names.append(name)
            log.info(f"Added feature: {name} ({len(data)} vertices, "
                     f"mean={data.mean():.4f}, std={data.std():.4f})")

        except Exception as e:
            log.warning(f"Failed to load {name} from {path}: {e}")
            continue

    if not darrays:
        log.error("No features were loaded!")
        sys.exit(1)

    # Create combined GIFTI
    combined = nib.gifti.GiftiImage(darrays=darrays)
    nib.save(combined, snakemake.output.combined)

    log.info(f"Combined {len(darrays)} features into {snakemake.output.combined}")
    log.info(f"Features: {', '.join(feature_names)}")

except Exception as e:
    log.exception(f"Error combining features: {e}")
    sys.exit(1)
