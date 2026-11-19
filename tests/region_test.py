#!/usr/bin/env python3
"""Live test: cut out every DP2 object in a small random region.

Picks a random point inside a known DP2 field, finds all dp2.Object rows
within --search-radius of it, then fetches a cutout for each one with a single
batch run (`lsst-cutout -t -o DIR BAND < ids`), or with one run per object
(--sequential, the old baseline), and checks that each output is FITS.
Prints per-stage timing statistics. Needs RSP_TOKEN. Exits 1 on any failure.
"""

import argparse
import importlib.machinery
import importlib.util
import math
import os
import random
import re
import statistics
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

TOOL = Path(__file__).resolve().parent.parent / "lsst-cutout"

# Centres of fields with DP2 deep coadds (from dp2.lsst.io tutorials).
FIELDS = {"ECDFS": (53.13, -28.10)}
FIELD_SPREAD = 0.2  # deg; random offset from the field centre

TIMING_LINE = re.compile(r"^lsst-cutout: (\S+)\s+([\d.]+)s")


def load_tool():
    loader = importlib.machinery.SourceFileLoader("lsst_cutout", str(TOOL))
    spec = importlib.util.spec_from_loader("lsst_cutout", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def collect(stderr, stages):
    """Add the -t timing lines in an lsst-cutout stderr to stages."""
    for line in stderr.splitlines():
        m = TIMING_LINE.match(line)
        if m:
            stages[m[1]].append(float(m[2]))


def random_point(rng):
    name, (ra0, dec0) = rng.choice(sorted(FIELDS.items()))
    r = FIELD_SPREAD * math.sqrt(rng.random())  # uniform over a disc
    theta = rng.uniform(0, 2 * math.pi)
    dec = dec0 + r * math.sin(theta)
    ra = ra0 + r * math.cos(theta) / math.cos(math.radians(dec0))
    return name, ra % 360, dec


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--band", default="r", choices="ugrizy")
    p.add_argument("--search-radius", type=float, default=0.01,
                   help="object search radius in deg (default: 0.01)")
    p.add_argument("--cutout-radius", type=float, default=0.01,
                   help="radius of each cutout in deg (default: 0.01)")
    p.add_argument("--seed", type=int, help="random seed (default: random)")
    p.add_argument("--max", type=int, help="cut out at most this many objects")
    p.add_argument("-j", "--jobs", type=int, default=1,
                   help="parallel cutout requests in batch mode (default: 1)")
    p.add_argument("--sequential", action="store_true",
                   help="one lsst-cutout run per object instead of a batch")
    p.add_argument("-o", "--outdir", help="where to keep the FITS files "
                   "(default: a new temporary directory)")
    a = p.parse_args()

    token = os.environ.get("RSP_TOKEN", "").strip()
    if not token:
        sys.exit("RSP_TOKEN not set")
    seed = a.seed if a.seed is not None else random.randrange(2**32)
    field, ra, dec = random_point(random.Random(seed))
    outdir = Path(a.outdir or tempfile.mkdtemp(prefix="lsst-cutout-test-"))
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"seed {seed}: {field} field, centre {ra:.6f} {dec:.6f}, "
          f"search radius {a.search_radius} deg, band {a.band}")

    tool = load_tool()
    t0 = time.monotonic()
    rows = tool.tap(
        "SELECT objectId, coord_ra, coord_dec FROM dp2.Object WHERE "
        "CONTAINS(POINT('ICRS', coord_ra, coord_dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {a.search_radius})) = 1",
        token, ("objectId", "coord_ra", "coord_dec"), "search")
    print(f"found {len(rows)} objects in {time.monotonic() - t0:.2f}s")
    if not rows:
        sys.exit("no objects found; try another --seed")
    rows = rows[:a.max] if a.max else rows

    stages, failures = defaultdict(list), []
    ids = [int(row["objectId"]) for row in rows]
    size = tool.radius_label(a.cutout_radius)
    t0 = time.monotonic()
    if a.sequential:
        for i, oid in enumerate(ids, 1):
            cmd = [sys.executable, str(TOOL), "-t", "-r", str(a.cutout_radius),
                   a.band, str(oid)]
            (outdir / str(oid)).mkdir(exist_ok=True)
            with open(outdir / str(oid) / f"{a.band}_{size}.fits", "wb") as f:
                proc = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE,
                                      text=True)
            collect(proc.stderr, stages)
            print(f"[{i}/{len(ids)}] {oid} exit {proc.returncode}")
    else:
        cmd = [sys.executable, str(TOOL), "-t", "--clobber",
               "-r", str(a.cutout_radius), "-j", str(a.jobs), "-o", str(outdir), a.band]
        proc = subprocess.run(cmd, input="".join(f"{i}\n" for i in ids),
                              capture_output=True, text=True)
        collect(proc.stderr, stages)
        for line in proc.stderr.splitlines():
            if not TIMING_LINE.match(line):
                print(line)
        print(f"batch exit {proc.returncode}, "
              f"{len(proc.stdout.splitlines())} files reported")
    wall = time.monotonic() - t0

    for oid in ids:
        out = outdir / str(oid) / f"{a.band}_{size}.fits"
        if not out.exists():
            failures.append((oid, "no output file"))
            continue
        with open(out, "rb") as f:
            if f.read(9) != b"SIMPLE  =":
                failures.append((oid, f"not FITS ({out.stat().st_size} bytes)"))

    print(f"\ncutouts in {outdir}")
    print(f"wall time {wall:.2f}s for {len(ids)} objects "
          f"({wall / len(ids):.2f}s per object)")
    print(f"{'stage':<9} {'n':>4} {'mean':>7} {'median':>7} {'max':>7} {'sum':>8}")
    for stage in ("object", "obscore", "datalink", "cutout", "total"):
        v = stages.get(stage)
        if v:
            print(f"{stage:<9} {len(v):4d} {statistics.mean(v):7.2f} "
                  f"{statistics.median(v):7.2f} {max(v):7.2f} {sum(v):8.2f}")
    if failures:
        print(f"\n{len(failures)} of {len(rows)} failed:")
        for oid, reason in failures:
            print(f"  {oid}: {reason}")
        sys.exit(1)
    print(f"\nall {len(ids)} cutouts ok")


if __name__ == "__main__":
    main()
