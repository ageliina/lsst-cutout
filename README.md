# lsst-cutout

Fetch image cutouts of LSST DP2 deep coadds from the
[Rubin Science Platform](https://data.lsst.cloud) as FITS files.

It is a single Python 3 script with no dependencies outside the standard library.

## Disclaimer

All credits go to Claude Code 2.1.327. Prompted by AV.

## Setup

Create an access token with the `read:tap` and `read:image` scopes at
<https://data.lsst.cloud/settings/tokens> and export it:

```sh
export RSP_TOKEN='gt-…'
```

## Usage

One target, FITS written to stdout:

```sh
lsst-cutout r 775563066640193057 > obj.fits      # by objectId
lsst-cutout r 53.13 -28.10      > pos.fits       # by RA Dec (ICRS degrees)
lsst-cutout -r 0.005 g 53.13 -28.10 > small.fits # cutout radius in degrees
```

Many targets: read from stdin, one `OBJECTID` or `RA DEC` per line, and write to a
directory. Blank lines and `#` comments are ignored, and commas count as spaces.

```sh
lsst-cutout -o cutouts/ r < targets.txt
```

`BAND` can name several bands, such as `gr` or `ugrizy`, and each band is fetched
separately. Each target gets its own directory, and each cutout is written to
`DIR/<name>/<band>_<radius>.fits`, where `<name>` is the objectId or
`<ra>_<dec>`. `<radius>` is the cutout radius in the
largest unit that keeps the number at least 1, e.g. `1.4deg`, `4.2arcmin` or
`36arcsec` (the default, 0.01°). Each path is printed on stdout as soon as the file
has been written. A target that fails is reported on stderr and the others carry
on.

| Option | Meaning |
| --- | --- |
| `BAND` | one or more of `ugrizy`, e.g. `r`, `gr` or `ugrizy` (required). Writing to stdout takes exactly one band |
| `-r RADIUS` | cutout radius in degrees (default 0.01). The cutout is a square with side 2 × RADIUS |
| `-o DIR` | batch mode: read targets from stdin and write the cutouts to DIR |
| `-j N` | number of cutouts fetched in parallel in batch mode (default 1; the server allows about 35 requests per minute) |
| `-p` | batch mode: also draw a gallery of the `BAND` cutouts for each target (see below) |
| `--clobber` | batch mode: overwrite existing output files (see below) |
| `-t` | print per-stage timings to stderr |

### Gallery (`-p`)

```sh
lsst-cutout -p -o cutouts/ ugrizy < targets.txt   # 6-band gallery
lsst-cutout -p -o cutouts/ gr < targets.txt       # g and r only
```

With `-p`, `DIR/<name>/<bands>_<radius>.png` shows the cutouts of the `BAND`
bands side by side.
Every panel has a 0.1′ scale bar. All panels share the colour scale of band `i`
if it is included, otherwise of the reddest band. A marker shows the object, and
the title gives `objectId`, `i_ra`, `i_dec` and the PSF magnitudes of the bands,
all from `dp2.Object`, plus g−r and i−z when both bands of a colour are included. For RA
Dec targets the title shows only the position. A band with no image is left out
of the gallery and reported on stderr. Each panel shows the object ± the cutout
radius. A cutout can come back smaller than requested when the object is near the
edge of its patch image, because the service only returns pixels inside that
image. The panel then stays centred on the object and shows the missing area in
grey; the FITS file keeps exactly what the service returned.
`-p` needs astropy, numpy and matplotlib. Without it, the tool uses only the
standard library.

### Existing files

By default, existing output files are kept, and only what is missing is fetched
or drawn. A target whose files all exist needs no network requests, and a missing
gallery is drawn from the FITS files already on disk. Kept files are still listed
on stdout, so the output is always the complete list of files, and a note on
stderr says how many were kept. `--clobber` fetches everything again and
overwrites it.

Exit status: 0 if every requested cutout succeeded, 1 if any failed, 2 for a usage error.

## How it works

1. **objectId targets:** one TAP query to `dp2.Object` (`objectId IN (…)`) gets the
   position, tract and patch of every object, plus the gallery columns with `-p`. One `ivoa.ObsCore` query then gets the
   deep coadds of those tracts, and each object gets the image of the patch it was
   measured on.
2. **RA Dec targets:** an `ivoa.ObsCore` query per target finds the images that
   contain the position, and the one whose centre is nearest is used.
3. The image's DataLink record gives the cutout service and its image ID. Records
   are cached, so targets on the same image share one record.
4. One SODA `cutout-sync` request per target. The service only accepts one cutout
   per request. The requests can run in parallel with `-j`, but the server
   limits clients to about 35 requests per minute, so the default is one at a time.

As a rough guide, a batch of 10 objectIds takes about 40 s: about 5 s for each of
the two TAP queries, plus about 3 s per cutout (about 18 s with `-j 4`).

## Limitations

- The tool only covers DP2 (`LSST.DP2`, `dp2.Object`) and deep coadds
  (`lsst.deep_coadd`). These are constants at the top of the script.
- Cutouts use the `cutout-sync` service, which returns the pixels, the WCS and the
  primary header. The mask, variance and PSF planes are not included.

## Tests

`tests/region_test.py` is a live test that needs `RSP_TOKEN`. It picks a random
point in the ECDFS field, finds every object within 0.01° of it, cuts them all out,
checks that each file is FITS, and prints a timing summary.

```sh
tests/region_test.py --max 10 --seed 42               # batch mode
tests/region_test.py --max 10 --seed 42 --sequential  # one run per object, for comparison
```

`tests/visual_test.py` is a visual check. It needs astropy, numpy and matplotlib.
It downloads the full deep coadd of one tract/patch (by default, the patch
containing the ECDFS centre), fetches service cutouts of a few bright objects on
that patch, and compares each one pixel by pixel with the same region of the
downloaded image. It writes `overview.png` (the patch with the cutout boxes drawn
on it) and `cutouts.png` (downloaded, service and difference images side by side).
The downloaded patch is kept and reused on later runs.

```sh
tests/visual_test.py                                   # writes into ./visual-test/
tests/visual_test.py --band g --tract 5063 --patch 34 -n 9
```
