#!/usr/bin/env python3
"""Visual test: compare service cutouts with cutouts of a downloaded patch.

Downloads the full deep_coadd image of one tract/patch in BAND (default: the
patch containing the ECDFS centre), fetches cutouts of N bright objects of that
patch with `lsst-cutout -o`, cuts the same pixels out of the downloaded image,
and compares them. Writes into OUTDIR:

  patch_<tract>_<patch>_<band>.fits   the full image (reused if present)
  cutouts/<objectId>/<band>_<radius>.fits  the service cutouts
  overview.png                        the patch with the cutout boxes
  cutouts.png                         per object: local | service | difference

Needs RSP_TOKEN, astropy, numpy and matplotlib. Exits 1 if any cutout differs
from the downloaded image.
"""

import argparse
import importlib.machinery
import importlib.util
import io
import os
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits, votable  # noqa: E402
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval  # noqa: E402
from astropy.wcs import WCS  # noqa: E402

TOOL = Path(__file__).resolve().parent.parent / "lsst-cutout"
ECDFS = (53.13, -28.10)


def load_tool():
    loader = importlib.machinery.SourceFileLoader("lsst_cutout", str(TOOL))
    spec = importlib.util.spec_from_loader("lsst_cutout", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def image_hdu(path):
    """(data, WCS) of the first HDU with a 2-d image and a celestial WCS."""
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim == 2:
                wcs = WCS(hdu.header)
                if wcs.has_celestial:
                    return np.asarray(hdu.data, dtype=float), wcs
    sys.exit(f"{path}: no 2-d image HDU with a celestial WCS")


def compare(data, wcs, path):
    """Compare a cutout file with the same pixels of the full image.

    Returns (x0, y0, local, cut, diff, max_abs_diff, ok), or an error string.
    """
    cut, cwcs = image_hdu(path)
    # sky position of the cutout's first pixel -> pixel in the full image
    x0, y0 = (int(round(float(v))) for v in
              wcs.world_to_pixel(cwcs.pixel_to_world(0, 0)))
    h, w = cut.shape
    local = data[max(y0, 0):y0 + h, max(x0, 0):x0 + w]
    if x0 < 0 or y0 < 0 or local.shape != cut.shape:
        return f"cutout at x={x0} y={y0} runs off the patch image"
    both = np.isfinite(local) & np.isfinite(cut)
    diff = np.where(both, cut - local, np.nan)
    maxdiff = float(np.nanmax(np.abs(diff))) if both.any() else float("nan")
    ok = bool(both.any()) and maxdiff == 0 and \
        np.array_equal(np.isfinite(local), np.isfinite(cut))
    return x0, y0, local, cut, diff, maxdiff, ok


def find_patch(tool, a, token):
    """(tract, patch, DataLink URL) of the image to download."""
    where = (f"lsst_tract = {a.tract} AND lsst_patch = {a.patch}"
             if a.tract is not None else
             f"CONTAINS(POINT('ICRS', {ECDFS[0]}, {ECDFS[1]}), s_region) = 1")
    rows = tool.tap(
        "SELECT lsst_tract, lsst_patch, access_url, s_ra, s_dec "
        f"FROM ivoa.ObsCore WHERE {where} "
        f"AND obs_collection = '{tool.COLLECTION}' "
        f"AND dataproduct_subtype = '{tool.SUBTYPE}' "
        f"AND lsst_band = '{a.band}'", token,
        ("lsst_tract", "lsst_patch", "access_url", "s_ra", "s_dec"), "obscore")
    if not rows:
        sys.exit("no matching deep_coadd image")
    r = min(rows, key=lambda r: tool.separation(*ECDFS, float(r["s_ra"]),
                                                float(r["s_dec"])))
    return int(r["lsst_tract"]), int(r["lsst_patch"]), r["access_url"]


def download(tool, dl_url, path, token):
    """Save the full image (DataLink '#this' link) to path."""
    links = votable.parse(io.BytesIO(tool.fetch(dl_url, token, stage="datalink"))
                          ).get_first_table().to_table()
    this = [row for row in links if row["semantics"] == "#this"]
    if not this:
        sys.exit(f"no #this link in DataLink response from {dl_url}")
    url = str(this[0]["access_url"])
    req = urllib.request.Request(url)
    if urllib.parse.urlsplit(url).netloc == urllib.parse.urlsplit(tool.API).netloc:
        req.add_header("Authorization", f"Bearer {token}")  # not to signed URLs
    tmp = path.with_suffix(".part")
    with tool.OPENER.open(req) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 20)
    tmp.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--band", default="r", choices="ugrizy")
    p.add_argument("--tract", type=int, help="with --patch: image to use")
    p.add_argument("--patch", type=int)
    p.add_argument("-n", type=int, default=6, help="objects (default: 6)")
    p.add_argument("--mag", type=float, default=22.0,
                   help="only objects brighter than this cModel mag (default: 22)")
    p.add_argument("--cutout-radius", type=float, default=0.002,
                   help="cutout radius in deg (default: 0.002)")
    p.add_argument("-o", "--outdir", default="visual-test",
                   help="output directory (default: ./visual-test)")
    a = p.parse_args()
    if (a.tract is None) != (a.patch is None):
        p.error("give both --tract and --patch, or neither")
    token = os.environ.get("RSP_TOKEN", "").strip() or sys.exit("RSP_TOKEN not set")
    tool = load_tool()
    out = Path(a.outdir)
    (out / "cutouts").mkdir(parents=True, exist_ok=True)

    # 1. the patch image
    tract, patch, dl_url = find_patch(tool, a, token)
    image = out / f"patch_{tract}_{patch}_{a.band}.fits"
    if image.exists():
        print(f"using existing {image}")
    else:
        print(f"downloading tract {tract} patch {patch} band {a.band} ...")
        download(tool, dl_url, image, token)
    print(f"{image}: {image.stat().st_size / 1e6:.1f} MB")
    data, wcs = image_hdu(image)

    # 2. bright objects measured on this patch, and their service cutouts
    rows = tool.tap(
        f"SELECT TOP {a.n} objectId, coord_ra, coord_dec FROM {tool.OBJECT_TABLE} "
        f"WHERE tract = {tract} AND patch = {patch} "
        f"AND {a.band}_cModelMag < {a.mag}", token,
        ("objectId", "coord_ra", "coord_dec"), "object")
    if not rows:
        sys.exit("no objects found; try a fainter --mag")
    ids = [int(r["objectId"]) for r in rows]
    proc = subprocess.run(
        [sys.executable, str(TOOL), "--clobber", "-r", str(a.cutout_radius),
         "-o", str(out / "cutouts"), a.band],
        input="".join(f"{i}\n" for i in ids), capture_output=True, text=True)
    sys.stderr.write(proc.stderr)

    # 3. the same pixels from the downloaded image
    results, bad = [], 0
    for oid in ids:
        size = tool.radius_label(a.cutout_radius)
        path = out / "cutouts" / str(oid) / f"{a.band}_{size}.fits"
        if not path.exists():
            print(f"{oid}: no service cutout")
            bad += 1
            continue
        r = compare(data, wcs, path)
        if isinstance(r, str):
            print(f"{oid}: {r}")
            bad += 1
            continue
        x0, y0, local, cut, diff, maxdiff, ok = r
        bad += not ok
        print(f"{oid}: {cut.shape[1]}x{cut.shape[0]} px at x={x0} y={y0}, "
              f"max |diff| {maxdiff:g} {'ok' if ok else 'MISMATCH'}")
        results.append((oid, x0, y0, local, cut, diff))

    # 4. pictures
    norm = ImageNormalize(data, interval=ZScaleInterval(), stretch=AsinhStretch())
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(data, origin="lower", cmap="gray", norm=norm)
    for oid, x0, y0, _, cut, _ in results:
        h, w = cut.shape
        ax.add_patch(plt.Rectangle((x0, y0), w, h, fill=False, color="red"))
        ax.text(x0, y0 + h, str(oid), color="red", fontsize=6, va="bottom")
    ax.set_title(f"tract {tract} patch {patch} band {a.band}")
    fig.savefig(out / "overview.png", dpi=150, bbox_inches="tight")

    if results:
        fig, axes = plt.subplots(len(results), 3, squeeze=False,
                                 figsize=(7.5, 2.5 * len(results)))
        for row, (oid, _, _, local, cut, diff) in zip(axes, results):
            n = ImageNormalize(local, interval=ZScaleInterval(),
                               stretch=AsinhStretch())
            for axis, img, title in zip(row, (local, cut, diff),
                                        ("downloaded patch", "service", "difference")):
                axis.imshow(img, origin="lower", cmap="gray",
                            norm=n if title != "difference" else None)
                axis.set_title(f"{title}\n{oid}" if title == "downloaded patch"
                               else title, fontsize=8)
                axis.set_xticks([])
                axis.set_yticks([])
        fig.tight_layout()
        fig.savefig(out / "cutouts.png", dpi=120)

    print(f"\nwrote {out / 'overview.png'} and {out / 'cutouts.png'}")
    if bad:
        sys.exit(f"{bad} of {len(ids)} cutouts missing or different")
    print(f"all {len(ids)} cutouts match the downloaded image")


if __name__ == "__main__":
    main()
