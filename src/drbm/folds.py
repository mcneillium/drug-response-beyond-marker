"""Cross-validation folds, made once and saved with a fingerprint.

Five-fold cross-validation repeated five times. Folds are stratified by tissue
so each has a similar mix of cancers. Saving them means every model in every
experiment is scored on exactly the same held-out cell lines.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

N_FOLDS = 5
N_REPEATS = 5


def make_folds(ids: pd.Index, strata: pd.Series) -> list[list[list[str]]]:
    """folds[repeat][fold] = the held-out cell line ids. Repeat r uses seed r."""
    ids = pd.Index(ids)
    labels = strata.loc[ids].to_numpy()
    out = []
    for repeat in range(N_REPEATS):
        splitter = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=repeat)
        out.append([sorted(ids[test].tolist()) for _, test in splitter.split(np.zeros(len(ids)), labels)])
    return out


def _fingerprint(content: dict) -> str:
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def save_folds(all_folds: dict[str, list], path: Path) -> None:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"{path} already exists. Folds are fixed.")
    content = {"version": 1, "n_folds": N_FOLDS, "n_repeats": N_REPEATS, "folds": all_folds}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({**content, "sha256": _fingerprint(content)}))


def load_folds(path: Path) -> dict[str, list]:
    saved = json.loads(Path(path).read_text())
    content = {key: value for key, value in saved.items() if key != "sha256"}
    if _fingerprint(content) != saved.get("sha256"):
        raise ValueError(f"{path} does not match its fingerprint; it has been altered.")
    return saved["folds"]


def iterate(folds_for_drug: list, ids: pd.Index):
    """Yield (repeat, fold, train ids, test ids)."""
    ids = pd.Index(ids)
    for repeat, repeat_folds in enumerate(folds_for_drug):
        for fold, test in enumerate(repeat_folds):
            test_index = pd.Index(test)
            yield repeat, fold, ids.difference(test_index), test_index
