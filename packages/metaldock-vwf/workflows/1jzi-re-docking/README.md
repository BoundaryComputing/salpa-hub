# 1JZI Re docking — importable workflow

`workflow.bcflow` is an importable BoCoFlow workflow for the six-node metaldock
pipeline on the 1JZI Re-complex case. Node configuration, inputs, and expected
results are documented in [`../1jzi_re_demo.md`](../1jzi_re_demo.md).

[`slide.html`](slide.html) is a self-contained 16:9 SVG slide illustrating the
pipeline (open in any browser).

## How this file was produced

It is a **real BoCoFlow export** (`exportedBy: BoCoFlow BF2`) — the workflow is
assembled, wired, and configured *through the live app UI*, then exported. It is
**not** hand-generated JSON.

This matters: a hand-emitted `react` section imports but its edges never paint —
the links target the flow-control port instead of the data `in-0` port, so each
link's last point resolves to `(0,0)` and `VPLinkModel.getSVGPath` suppresses the
line. A genuine app export wires `out-0 → in-0` with render-compatible geometry,
so the edges show on the canvas.

Regenerate (rebuild in the UI + re-export) with the Playwright spec — this is the
always-correct generator:

```bash
cd operations
# BF2 stack must be running (admin UI :18001, orchestrator API :18000) with the
# metaldock-vwf package installed. Optionally point the QM node at ORCA:
METALDOCK_E2E_ORCA=/abs/orca_6_1_1_dir \
  npx playwright test export-showcase.spec.ts
```

It creates "MetalDock 1JZI Re Pipeline" in the app, configures all six nodes
(inputs from `collect/MetalDock/...` example data, `rel:` output dirs under the
`operations/e2e_run` working path), saves, then exports and writes this file —
asserting along the way that the export has 6 nodes, 5 edges, all wired
`out-0 → in-0`.

The `abs:` input refs and working path baked into the committed file point at
this machine — adjust them (or re-run the spec with different inputs) for another
environment.

## Verification status

Round-trip verified live: importing this file into the BoCoFlow GUI renders all
six nodes **and the five connecting edges** (confirmed via Playwright —
`[data-linkid]` link elements present on the canvas).
