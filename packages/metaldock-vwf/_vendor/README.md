# _vendor/

A vendored snapshot of `metaldock_modules` (from the repo's `src/metaldock_modules`)
so each node is self-contained when installed into a BoCoFlow node registry. The
node resolver (`_ensure_metaldock_modules` in every `node.py`) tries, in order:
`import metaldock_modules` → `METALDOCK_SRC` env → bundled `scripts/` → this
`_vendor/` → the repo `src/` (dev fallback).

Re-sync after changing the source modules:

    rsync -a --exclude=__pycache__ --exclude='*.pyc' \
        ../../../src/metaldock_modules ./
