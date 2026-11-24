"""The bundled template hands over finished inputs.

Up to 1.2.11, `workflows/pdbmdauto-pipeline.json` kept the values it had been tested
with. Every step's case name was `e2e_4z8j` and step 1 wrote under `pdbmdauto-e2e-full`,
so a workflow made from it for another protein still ran as "e2e_4z8j". Every step was
stored as "unconfigured", though nothing in it needs setting. The steps stood in a
staircase 100 px apart, so their cards overlapped. And two steps were both "GMX MD
Relaxation", so a run that failed at one did not say which. Salpa 0.34 also reads the
inputs a template declares (`template_info.inputs`) and fills them in from what a person
asks for when they start.

This checks the template as it ships:

- the four copies a template keeps of each setting agree: the graph's `options`, which
  is what runs, and on the step's canvas model `options.config`, which a load reads,
  `config` and `options.option_types[<setting>].value`;
- no setting's value holds "e2e" or "test";
- every step is stored as configured;
- `inputs` declares step 1's PDB ID, of kind `pdb_id`;
- the case name is empty on every step, so step 1 names the case after the PDB ID and
  each later step takes the case from the step before it, and step 1 writes under
  `pdbmdauto/`;
- no two cards' boxes overlap, each reckoned at 240 x 80 px, and the cards read in run
  order, left to right and then row by row;
- nothing on the canvas is stored as selected;
- every step's name is its own.

Run from the package root: pixi run test
"""

import json
import re
from pathlib import Path

import pytest

TEMPLATE_PATH = (
    Path(__file__).resolve().parents[1] / "workflows" / "pdbmdauto-pipeline.json"
)
TEMPLATE = json.loads(TEMPLATE_PATH.read_text(encoding="utf-8"))

# A card's box on the canvas, as an estimate: an unchosen card is at most 240 px wide
# and about 64 px high
CARD_WIDTH = 240
CARD_HEIGHT = 80
# Words that mark a value as left from a test run
TEST_LIKE = re.compile(r"e2e|test", re.IGNORECASE)
# Salpa's rule for a PDB ID: four characters, the first a digit
PDB_ID = re.compile(r"[1-9][A-Za-z0-9]{3}")


def _steps():
    return TEMPLATE["graph"]["nodes"]


def _models():
    """The canvas models, by step id."""
    layers = [
        layer
        for layer in TEMPLATE["react"]["layers"]
        if layer.get("type") == "diagram-nodes"
    ]
    assert len(layers) == 1, "the template has one layer of steps"
    return layers[0]["models"]


def _run_order():
    """The steps in the order a run takes them: each after every step linked into it."""
    steps = {step["id"]: step for step in _steps()}
    into = {step_id: 0 for step_id in steps}
    out = {step_id: [] for step_id in steps}
    for link in TEMPLATE["graph"]["links"]:
        out[link["source"]].append(link["target"])
        into[link["target"]] += 1
    free = [step_id for step_id in steps if into[step_id] == 0]
    order = []
    while free:
        step_id = free.pop(0)
        order.append(steps[step_id])
        for after in out[step_id]:
            into[after] -= 1
            if into[after] == 0:
                free.append(after)
    assert len(order) == len(steps), "every step is reached from the first"
    return order


def _copies(step):
    """The four copies of each of a step's settings: {setting: {copy: value}}."""
    model = _models()[step["id"]]
    copies = {}
    for name, values in (
        ("graph options", step["options"]),
        ("options.config", model["options"]["config"]),
        ("config", model["config"]),
    ):
        for setting, value in values.items():
            copies.setdefault(setting, {})[name] = value
    for setting, spec in model["options"]["option_types"].items():
        copies.setdefault(setting, {})["option_types value"] = spec.get("value")
    return copies


STEP_IDS = [step["id"] for step in _steps()]
STEP_NAMES = {step["id"]: step["name"] for step in _steps()}


def _step(step_id):
    return next(step for step in _steps() if step["id"] == step_id)


def test_the_check_finds_what_it_checks():
    """Without this, a template read as empty would pass everything below."""
    assert len(_steps()) == 11
    assert set(_models()) == set(
        STEP_IDS
    ), "each step has its canvas model, and each model its step"
    assert len(TEMPLATE["graph"]["links"]) == 10
    assert sum(len(_copies(step)) for step in _steps()) >= 70


@pytest.mark.parametrize("step_id", STEP_IDS, ids=lambda step_id: STEP_NAMES[step_id])
def test_the_four_copies_of_every_setting_agree(step_id):
    for setting, copies in _copies(_step(step_id)).items():
        assert (
            len(copies) == 4
        ), f"{setting} is kept in {sorted(copies)}, not in all four copies"
        assert (
            len({json.dumps(value) for value in copies.values()}) == 1
        ), f"{setting}'s copies disagree: {copies}"


def test_no_setting_holds_a_value_left_from_a_test_run():
    found = [
        f"{STEP_NAMES[step['id']]}: {setting} ({where}) = {value!r}"
        for step in _steps()
        for setting, copies in _copies(step).items()
        for where, value in copies.items()
        if isinstance(value, str) and TEST_LIKE.search(value)
    ]
    assert not found, "\n".join(found)


@pytest.mark.parametrize("step_id", STEP_IDS, ids=lambda step_id: STEP_NAMES[step_id])
def test_every_step_is_stored_as_configured(step_id):
    assert _step(step_id)["status"] == "configured"
    assert _models()[step_id]["options"]["status"] == "configured"


def test_the_template_declares_step_1s_pdb_id_as_its_input():
    first = _run_order()[0]
    assert first["node_key"] == "PdbFastaBiopython"
    assert (
        first["options"]["input_mode"] == "pdb_id"
    ), "step 1 fetches its structure by PDB ID"
    assert TEMPLATE["template_info"].get("inputs") == [
        {"step": first["id"], "option": "pdb_id", "kind": "pdb_id", "label": "PDB ID"}
    ]
    spec = _models()[first["id"]]["options"]["option_types"]["pdb_id"]
    # Salpa fills in only a setting that holds text
    assert spec["type"] == "string"
    assert spec["label"] == "PDB ID"
    # The template's example, which a workflow keeps when nothing else is given
    assert PDB_ID.fullmatch(first["options"]["pdb_id"]), first["options"]["pdb_id"]


def test_the_case_and_its_folder_follow_the_pdb_id():
    named = {
        STEP_NAMES[step["id"]]: step["options"]["case_name"]
        for step in _steps()
        if step["options"]["case_name"]
    }
    assert (
        not named
    ), f"a case name set on a step stops it following the PDB ID: {named}"
    assert _run_order()[0]["options"]["output_dir"] == "rel:pdbmdauto"


def _boxes():
    models = _models()
    boxes = {}
    for step_id in STEP_IDS:
        x, y = models[step_id]["x"], models[step_id]["y"]
        assert isinstance(x, (int, float)) and isinstance(
            y, (int, float)
        ), f"{STEP_NAMES[step_id]} has no place"
        boxes[step_id] = (x, y, x + CARD_WIDTH, y + CARD_HEIGHT)
    return boxes


def test_no_two_cards_overlap():
    boxes = _boxes()
    overlapping = [
        f"{STEP_NAMES[a]} {boxes[a][:2]} and {STEP_NAMES[b]} {boxes[b][:2]}"
        for i, a in enumerate(STEP_IDS)
        for b in STEP_IDS[i + 1 :]
        if boxes[a][0] < boxes[b][2]
        and boxes[b][0] < boxes[a][2]
        and boxes[a][1] < boxes[b][3]
        and boxes[b][1] < boxes[a][3]
    ]
    assert not overlapping, "\n".join(overlapping)


def test_the_cards_read_in_run_order():
    """Left to right, then row by row: the steps sorted by row, and then by their place
    in the row, are in the run's order."""
    models = _models()
    by_place = sorted(
        STEP_IDS, key=lambda step_id: (models[step_id]["y"], models[step_id]["x"])
    )
    assert by_place == [step["id"] for step in _run_order()]


def test_nothing_on_the_canvas_is_stored_as_selected():
    """A card the template stores as selected is drawn selected when the workflow opens
    in an app that reads the flag: neither of a card's two flags (its model's and its
    options'), nor a link's."""
    selected = [
        f"{STEP_NAMES[step_id]}: {where}"
        for step_id, model in _models().items()
        for where, flag in (
            ("selected", model.get("selected")),
            ("options.selected", (model.get("options") or {}).get("selected")),
        )
        if flag
    ]
    links = [
        layer
        for layer in TEMPLATE["react"]["layers"]
        if layer["type"] == "diagram-links"
    ]
    assert len(links) == 1, "the template has one layer of links"
    selected += [
        f"the link {link_id}"
        for link_id, link in links[0]["models"].items()
        if link.get("selected")
    ]
    assert not selected, "\n".join(selected)


def test_every_step_has_its_own_name():
    for step in _steps():
        assert (
            _models()[step["id"]]["options"]["name"] == step["name"]
        ), "the canvas shows the name that runs"
    names = [step["name"] for step in _steps()]
    assert len(set(names)) == len(names), sorted(
        name for name in names if names.count(name) > 1
    )
