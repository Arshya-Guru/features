#!/usr/bin/env python3
"""EpLink Cortical Feature Extraction Pipeline.

A snakebids BIDS App that extracts vertex-wise cortical surface features
from FreeSurfer-style derivatives (DeepPrep, FreeSurfer, FastSurfer).
"""
from pathlib import Path

from snakebids import bidsapp, plugins

app = bidsapp.app(
    [
        # snakemake_dir = project root (where config/ and workflow/ live)
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.BidsValidator(),
    ]
)


def get_parser():
    """Expose parser for sphinx doc generation."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
