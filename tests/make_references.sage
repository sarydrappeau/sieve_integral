r"""
Regenerate ``tests/reference_values.json``. From the repository root::

    sage tests/make_references.sage                 # every example
    sage tests/make_references.sage 'stadlmann/*'   # only those, the rest kept

Each entry holds the ball the example returns at its own default precision, as
exact midpoint and radius strings, together with that precision. The file also
records where the values come from: commit, branch, date, Sage and LattE.

These are regression references, not certified values: they say what the code
gave, and the tests check that a later value still overlaps them. Regenerating
is deliberate, and the commit that does it must say why each entry moved.
"""
import fnmatch
import json
import sys
from time import time

load("tests/harness.sage")

patterns = sys.argv[1:]
names = [name for name in example_names()
         if not patterns or any(fnmatch.fnmatch(name, p) for p in patterns)]

try:
    previous = read_references()
except FileNotFoundError:
    previous = {"examples": {}}

entries = dict(previous["examples"])
for name in names:
    start = time()
    value = example_value(name)
    entry = ball_to_json(value)
    entry["precision"] = str(default_precision(name))
    old = entries.get(name)
    entries[name] = entry
    moved = "" if old is None else (
        "  (unchanged)" if old["mid"] == entry["mid"] and old["rad"] == entry["rad"]
        else f"  (was {old['value']})")
    print(f"{name:20s} {time() - start:7.2f}s  {entry['value']}{moved}", flush=True)

record = {"provenance": provenance(),
          "note": ("Regression references: what the code gave on this commit, "
                   "at each example's own precision. Not certified values."),
          "examples": {name: entries[name] for name in sorted(entries)}}

with open(REFERENCES, "w") as f:
    f.write('{"provenance": ' + json.dumps(record["provenance"], sort_keys=True)
            + ',\n "note": ' + json.dumps(record["note"])
            + ',\n "examples": {\n')
    f.write(",\n".join(f'  {json.dumps(name)}: {json.dumps(record["examples"][name], sort_keys=True)}'
                       for name in record["examples"]))
    f.write("\n }}\n")
print(f"\n{len(record['examples'])} references written to {REFERENCES}")
