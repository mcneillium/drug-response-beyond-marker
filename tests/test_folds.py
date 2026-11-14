import json

import pandas as pd
import pytest

from drbm import folds


def test_every_cell_line_is_held_out_exactly_once_per_repeat(toy):
    made = folds.make_folds(toy["ids"], toy["meta"]["tissue_group"])
    assert len(made) == folds.N_REPEATS
    for repeat in made:
        assert len(repeat) == folds.N_FOLDS
        held_out = [i for fold in repeat for i in fold]
        assert sorted(held_out) == sorted(toy["ids"])


def test_folds_are_stratified_by_tissue(toy):
    made = folds.make_folds(toy["ids"], toy["meta"]["tissue_group"])
    overall = (toy["meta"]["tissue_group"] == "Lung").mean()
    for fold in made[0]:
        assert abs((toy["meta"].loc[fold, "tissue_group"] == "Lung").mean() - overall) < 0.05


def test_folds_are_reproducible_and_repeats_differ(toy):
    first = folds.make_folds(toy["ids"], toy["meta"]["tissue_group"])
    second = folds.make_folds(toy["ids"], toy["meta"]["tissue_group"])
    assert first == second
    assert first[0] != first[1]


def test_iterate_gives_disjoint_train_and_test(toy):
    made = folds.make_folds(toy["ids"], toy["meta"]["tissue_group"])
    count = 0
    for _, _, train, test in folds.iterate(made, toy["ids"]):
        assert len(train.intersection(test)) == 0
        assert len(train) + len(test) == len(toy["ids"])
        count += 1
    assert count == folds.N_REPEATS * folds.N_FOLDS


def test_saved_folds_cannot_be_overwritten_or_altered(toy, tmp_path):
    path = tmp_path / "folds.json"
    made = {"drug": folds.make_folds(toy["ids"], toy["meta"]["tissue_group"])}
    folds.save_folds(made, path)
    assert folds.load_folds(path) == made
    with pytest.raises(FileExistsError):
        folds.save_folds(made, path)
    saved = json.loads(path.read_text())
    saved["folds"]["drug"][0][0][0] = "TAMPERED"
    path.write_text(json.dumps(saved))
    with pytest.raises(ValueError):
        folds.load_folds(path)


def test_the_real_folds_file_matches_its_fingerprint():
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "splits" / "folds_v1.json"
    if not path.exists():
        pytest.skip("folds not made yet")
    loaded = folds.load_folds(path)
    for drug_folds in loaded.values():
        for repeat in drug_folds:
            ids = [i for fold in repeat for i in fold]
            assert len(ids) == len(set(ids))
