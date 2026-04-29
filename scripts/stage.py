"""Shared stage() helper for Python (scanpy) modules.

Mirrors scripts/stage.R: emit start, run block, emit end with elapsed_s,
then force gc.collect() so the next stage's RSS baseline is clean.

Usage:
    from stage import stage, stage_timings_df

    with stage("load"):
        ad = sc.read_h5ad(path)
"""
from __future__ import annotations

import gc
import time
from contextlib import contextmanager

from omnibench_logger import emit   # canonical schema


_timings: dict[str, float] = {}


@contextmanager
def stage(name: str):
    emit(name, "start")
    t0 = time.monotonic()
    try:
        yield
    finally:
        elapsed = time.monotonic() - t0
        emit(name, "end", attrs={"elapsed_s": elapsed})
        _timings[name] = elapsed
        gc.collect()


def stage_timings_df(dataset: str | None = None):
    import pandas as pd
    return pd.DataFrame(
        {"stage": list(_timings), "elapsed_s": list(_timings.values()),
         "dataset": dataset}
    )
