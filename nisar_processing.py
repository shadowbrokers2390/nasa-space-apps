# These are the cells run in Google Colab on 2026-11-14 for NASA Space Apps 2026, in order.
# AI acknowledgment: this code was written with Claude (Anthropic), an AI assistant.
# Paste each step into its own Colab cell. get_ipython().system(...) lines are Colab shell commands.

# NISAR ground movement at Portuguese Bend - NASA Space Apps 2026, "Dancing with the SARs"
# Code written with help from Claude (Anthropic), an AI assistant, during the event.

# Step 1. What this does: installs NASA's data-access library and logs in to Earthdata.
get_ipython().system('pip -q install earthaccess')
import earthaccess
earthaccess.login()

# Step 2. What this does: sets up the site and finds every NISAR scene that covers the landslide.
import time, json
from datetime import datetime
import numpy as np, h5py
from pyproj import Transformer
from scipy import ndimage
from google.colab import files

W, S, E, N = -118.45, 33.70, -118.30, 33.80      # site box around Palos Verdes (lon/lat)
SLIDE = (-118.365, 33.742)                        # approximate center of the Portuguese Bend landslide
STEP = 0.0008                                     # output grid spacing in degrees, about 80 m
COH_MIN = 0.3                                     # hide pixels with coherence (quality) below this
ncols, nrows = round((E - W) / STEP), round((N - S) / STEP)
LON, LAT = np.meshgrid(W + STEP * (np.arange(ncols) + 0.5), N - STEP * (np.arange(nrows) + 0.5))

allg = earthaccess.search_data(short_name="NISAR_L2_GUNW_PROVISIONAL_V1", point=SLIDE, count=400)
name = lambda g: g["umm"]["GranuleUR"]
# Two viewing directions. Each key is _track_direction_frame_ inside the scene name.
TRACKS = {"asc": "_034_A_019_", "desc": "_013_D_071_"}
for t, key in TRACKS.items():
    print(t, len([g for g in allg if key in name(g)]), "date pairs")

# Step 3. What this does: turns ONE date pair into centimeters of line-of-sight movement on the site grid.
# It reads only the site's pixels, hides low-quality pixels, converts radar phase to cm,
# and measures against the typical (median) pixel, because InSAR is a relative measurement.
UNW = "science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/"

def process(g):
    h = h5py.File(earthaccess.open([g])[0], "r")
    wl = 299792458 / float(h["science/LSAR/GUNW/grids/frequencyA/centerFrequency"][()])   # wavelength = c / frequency
    xx, yy = h[UNW + "xCoordinates"][:], h[UNW + "yCoordinates"][:]
    tr = Transformer.from_crs(4326, int(h[UNW + "projection"][()]), always_xy=True)
    X, Y = tr.transform(LON, LAT)
    r = np.rint((Y - yy[0]) / (yy[1] - yy[0])).astype(int)
    c = np.rint((X - xx[0]) / (xx[1] - xx[0])).astype(int)
    ok = (r >= 0) & (r < len(yy)) & (c >= 0) & (c < len(xx))
    a0, a1, b0, b1 = r[ok].min(), r[ok].max() + 1, c[ok].min(), c[ok].max() + 1
    ph = h[UNW + "unwrappedPhase"][a0:a1, b0:b1]
    co = h[UNW + "coherenceMagnitude"][a0:a1, b0:b1]
    cc = h[UNW + "connectedComponents"][a0:a1, b0:b1]
    look = h["science/LSAR/identification/lookDirection"][()].decode()
    h.close()
    good = np.isfinite(ph) & (cc > 0) & (co >= COH_MIN)
    los = np.where(good, -wl / (4 * np.pi) * ph * 100, np.nan)   # cm, positive = toward the satellite
    los -= np.nanmedian(los)
    out = np.full(LON.shape, np.nan); out[ok] = los[r[ok] - a0, c[ok] - b0]
    return out, look, wl

# Step 4. What this does: runs Step 3 for every date pair in both viewing directions (about 12 minutes).
TS = {}
for t, key in TRACKS.items():
    TS[t] = []
    for g in sorted([g for g in allg if key in name(g)], key=name):
        p = name(g).split("_")
        los, look, wavelength = process(g)
        TS[t].append({"d1": p[11][:8], "d2": p[13][:8], "los": los, "look": look})

# Step 5. What this does: defines the slide area from all the data, then saves everything the web map needs into ts.js.
import warnings; warnings.simplefilter("ignore", RuntimeWarning)   # ocean pixels are all-empty, which is expected

# Slide area = largest connected area where the typical (median) ascending date pair moved more than 1 cm away.
typical = np.nanmedian(np.array([p["los"] for p in TS["asc"]]), axis=0)
lab, n = ndimage.label(typical <= -1.0)
slide_mask = lab == (1 + int(np.argmax(np.bincount(lab.ravel())[1:])))

def pack(p):
    los = p["los"]; inside_ = los[slide_mask]; outside_ = los[~slide_mask]
    has_in, has_out = np.isfinite(inside_).any(), np.isfinite(outside_).any()
    day = lambda d: datetime.strptime(d, "%Y%m%d")
    return {"d1": p["d1"], "d2": p["d2"], "days": (day(p["d2"]) - day(p["d1"])).days,
            "mm": [int(round(v * 10)) if np.isfinite(v) else None for v in los.ravel()],       # millimeters
            "med_cm": round(float(np.nanmedian(inside_)), 2) if has_in else 0,                 # slide area median
            "mask_good_pct": int(round(100 * np.isfinite(inside_).mean())),                    # how much of it had data
            "noise_cm": [round(float(v), 1) for v in np.nanpercentile(outside_, [5, 95])] if has_out else [0, 0]}

INFO = {"asc": {"pass": "Ascending", "track": 34, "frame": 19, "from": "east", "toward": "west"},
        "desc": {"pass": "Descending", "track": 13, "frame": 71, "from": "west", "toward": "east"}}
TSOUT = {
    "bounds": [[N - nrows * STEP, W], [N, W + ncols * STEP]], "step_deg": STEP, "nrows": nrows, "ncols": ncols,
    "mask": [int(v) for v in np.flatnonzero(slide_mask)],
    "meta": {"product": "NISAR L2 GUNW (collection NISAR_L2_GUNW_PROVISIONAL_V1)", "status": "Provisional",
             "wavelength_m": round(wavelength, 4), "coherence_min": COH_MIN, "color_scale_cm": 4.0, "pixel_m": 80,
             "sign": "positive = toward satellite, negative = away",
             "reference": "each date pair is relative to the median of its good pixels in the map area",
             "mask_rule": "the largest connected area where the typical (median) ascending date pair moved more than 1 cm away from the satellite",
             "generated": "2026-11-14, NASA Space Apps Challenge. Processing code written with AI assistance (Claude, Anthropic); values come from NISAR data."},
    "tracks": [dict(INFO[t], id=t, look=TS[t][0]["look"], pairs=[pack(p) for p in TS[t]]) for t in ["asc", "desc"]],
}
txt = "// NISAR line-of-sight displacement for Portuguese Bend. Processing code written with AI assistance (Claude, Anthropic).\nconst TS = " + json.dumps(TSOUT, allow_nan=False)
open("ts.js", "w").write(txt)
for t in TSOUT["tracks"]:
    ok = [p for p in t["pairs"] if p["mask_good_pct"] >= 50]
    print(t["id"], "pairs used", len(ok), "of", len(t["pairs"]), "| per pair cm:", min(p["med_cm"] for p in ok), "to", max(p["med_cm"] for p in ok), "| total cm:", round(sum(p["med_cm"] for p in ok), 1))
print("slide area cells:", int(slide_mask.sum()), "| ts.js bytes:", len(txt), "| look:", {p["look"] for t in TS for p in TS[t]})

# Step 6. What this does: downloads the data file for the web map.
files.download("ts.js")

# Step 7. What this does: repeats Steps 2-4 for two more places where the ground is known to sink (about 20 minutes).
# These areas are much bigger than the landslide, so the output grid is coarser (about 330 m).
# shortcut: each grid cell takes the nearest 80 m radar pixel instead of averaging; average if the maps look speckled.
def run_site(point, box, step, tracks):
    global LON, LAT                      # process() reads the current site's grid
    w, s, e, n = box
    nc_, nr_ = round((e - w) / step), round((n - s) / step)
    LON, LAT = np.meshgrid(w + step * (np.arange(nc_) + 0.5), n - step * (np.arange(nr_) + 0.5))
    found = earthaccess.search_data(short_name="NISAR_L2_GUNW_PROVISIONAL_V1", point=point, count=400)
    series = {}
    for t, key in tracks.items():
        series[t] = []
        for g in sorted([g for g in found if key in name(g)], key=name):
            p = name(g).split("_")
            try:
                los, look, wl = process(g)
            except Exception as err:     # a scene that cannot be read is skipped, and listed
                SKIPPED.append(name(g)[:70] + " " + type(err).__name__); continue
            series[t].append({"d1": p[11][:8], "d2": p[13][:8], "los": los, "look": look})
    return series

SITES = {
    "mexico":   dict(point=(-99.05, 19.42),  box=(-99.30, 19.20, -98.85, 19.60),   step=0.003, tracks={"asc": "_019_A_011_", "desc": "_113_D_079_"}),
    "corcoran": dict(point=(-119.56, 36.10), box=(-119.85, 35.90, -119.30, 36.35), step=0.003, tracks={"asc": "_034_A_020_", "desc": "_114_D_070_"}),
}
SKIPPED = []
SERIES = {"pv": TS}                      # Portuguese Bend, already processed in Step 4
for sid, cfg in SITES.items():
    SERIES[sid] = run_site(**cfg)

# Step 8. What this does: checks the quality of each new site before anything is shown on the map.
# For the biggest sinking patches it prints where they are, how consistent they are from pair to pair,
# and whether the other viewing direction agrees (sinking should look similar from both sides).
print("skipped scenes:", SKIPPED)
for sid in ["mexico", "corcoran"]:
    w, s, e, n = SITES[sid]["box"]; step = SITES[sid]["step"]
    typ = {t: np.nanmedian(np.array([p["los"] for p in pairs]), axis=0) for t, pairs in SERIES[sid].items()}
    for t, pairs in SERIES[sid].items():
        spread = np.round(np.median([np.nanpercentile(p["los"], [5, 95]) for p in pairs if np.isfinite(p["los"]).any()], axis=0), 1)
        print(sid, t, "| pairs", len(pairs), pairs[0]["d1"], "to", pairs[-1]["d2"], "| typical spread of one pair, cm:", spread)
        lab, k = ndimage.label(typ[t] <= -0.5)
        sizes = np.bincount(lab.ravel())[1:]
        for i in np.argsort(sizes)[::-1][:3]:
            m = lab == i + 1
            rr, cc = np.where(m)
            meds = [np.nanmedian(p["los"][m]) for p in pairs if np.isfinite(p["los"][m]).mean() > 0.5]
            other = [o for o in typ if o != t][0]
            print("   patch cells", int(m.sum()), "| center lon/lat", round(w + (cc.mean() + 0.5) * step, 3), round(n - (rr.mean() + 0.5) * step, 3),
                  "| typical cm", round(float(np.nanmedian(typ[t][m])), 2), "| pairs moving away:", sum(1 for v in meds if v < 0), "of", len(meds),
                  "| total cm", round(float(np.sum(meds)), 1), "| other direction typical cm", round(float(np.nanmedian(typ[other][m])), 2))

# Step 9. What this does: saves one data file per place for the web map (replaces the single ts.js from Step 5).
# The Central Valley site is NOT exported: Step 8 showed no sinking patch that was consistent from pair to pair
# or that both viewing directions agreed on.
PLACES = {
    "pv": dict(box=(W, S, E, N), step=STEP, threshold=-1.0, pixel_m=80, total_scale_cm=60.0, default_total=False,
        place="Portuguese Bend, California", area_name="slide area", near_label="Roads",
        tracks={"asc": dict(track=34, frame=19), "desc": dict(track=13, frame=71)},
        cause_note="The location matches the known Portuguese Bend landslide, but this data alone does not show why the ground moved.",
        engineer=["Visual survey: new pavement cracks, offset curbs, tilted poles, broken utility lines.",
                  "GNSS or survey monuments: full 3D movement of fixed points. Projected onto the radar's viewing direction, this should match the map.",
                  "Inclinometers in boreholes: how deep the sliding surface is.",
                  "Crack gauges and extensometers: opening across the slide's edges.",
                  "Piezometers: groundwater pressure, which often drives slide speed.",
                  "Tiltmeters on structures: whether a building or wall is rotating."]),
    "mexico": dict(box=SITES["mexico"]["box"], step=SITES["mexico"]["step"], threshold=-0.5, pixel_m=330, total_scale_cm=20.0, default_total=True,
        place="Mexico City", area_name="sinking area", near_label="Major roads, rail lines and runways",
        tracks={"asc": dict(track=19, frame=11), "desc": dict(track=113, frame=79)},
        cause_note="Mexico City is widely reported to sink as groundwater is pumped from the soft lake clays beneath it, but this data alone does not show why the ground moved.",
        engineer=["Leveling surveys and GNSS stations: how fast each spot sinks in true vertical terms.",
                  "Extensometers and piezometers in wells: how much the clay layers are compacting and how groundwater pressure is changing.",
                  "Uneven sinking: where the speed changes over a short distance, check metro tracks, viaducts, pipes and runways for tilting and cracks.",
                  "Drain and sewer slopes: sinking can flatten or reverse the slope that makes water flow.",
                  "Tilt and foundation surveys on buildings near the edges of the sinking area."]),
}
VIEW = {"asc": {"pass": "Ascending", "from": "east", "toward": "west"}, "desc": {"pass": "Descending", "from": "west", "toward": "east"}}

def export(sid, cfg, series):
    w, s, e, n = cfg["box"]; step = cfg["step"]
    nc_, nr_ = round((e - w) / step), round((n - s) / step)
    # The moving area = largest connected area where the typical (median) ascending pair is beyond the threshold.
    typical = np.nanmedian(np.array([p["los"] for p in series["asc"]]), axis=0)
    lab, k = ndimage.label(typical <= cfg["threshold"])
    area = lab == (1 + int(np.argmax(np.bincount(lab.ravel())[1:])))
    day = lambda d: datetime.strptime(d, "%Y%m%d")
    def pack(p):
        los = p["los"]; inside_ = los[area]; outside_ = los[~area]
        has_in, has_out = np.isfinite(inside_).any(), np.isfinite(outside_).any()
        return {"d1": p["d1"], "d2": p["d2"], "days": (day(p["d2"]) - day(p["d1"])).days,
                "mm": [int(round(v * 10)) if np.isfinite(v) else None for v in los.ravel()],
                "med_cm": round(float(np.nanmedian(inside_)), 2) if has_in else 0,
                "mask_good_pct": int(round(100 * np.isfinite(inside_).mean())),
                "noise_cm": [round(float(v), 1) for v in np.nanpercentile(outside_, [5, 95])] if has_out else [0, 0]}
    def total(pairs):
        # Sum of all pairs per cell. A cell needs data in at least 80% of pairs; its missing pairs count as zero.
        stack = np.array([p["los"] for p in pairs]); have = np.isfinite(stack).sum(axis=0)
        tot = np.where(have >= 0.8 * len(pairs), np.nansum(stack, axis=0), np.nan)
        return [int(round(v * 10)) if np.isfinite(v) else None for v in tot.ravel()]
    out = {
        "bounds": [[n - nr_ * step, w], [n, w + nc_ * step]], "step_deg": step, "nrows": nr_, "ncols": nc_,
        "mask": [int(v) for v in np.flatnonzero(area)],
        "meta": {"product": "NISAR L2 GUNW (collection NISAR_L2_GUNW_PROVISIONAL_V1)", "status": "Provisional",
                 "wavelength_m": round(wavelength, 4), "coherence_min": COH_MIN, "color_scale_cm": 4.0,
                 "total_scale_cm": cfg["total_scale_cm"], "default_total": cfg["default_total"], "pixel_m": cfg["pixel_m"],
                 "place": cfg["place"], "area_name": cfg["area_name"], "near_label": cfg["near_label"],
                 "cause_note": cfg["cause_note"], "engineer": cfg["engineer"],
                 "sign": "positive = toward satellite, negative = away",
                 "reference": "each date pair is relative to the median of its good pixels in the map area",
                 "mask_rule": "the largest connected area where the typical (median) ascending date pair moved more than %s cm away from the satellite" % abs(cfg["threshold"]),
                 "generated": "2026-11-14, NASA Space Apps Challenge. Processing code written with AI assistance (Claude, Anthropic); values come from NISAR data."},
        "tracks": [dict(VIEW[t], **cfg["tracks"][t], id=t, look=series[t][0]["look"], pairs=[pack(p) for p in series[t]], total_mm=total(series[t])) for t in ["asc", "desc"]],
    }
    txt = "// NISAR line-of-sight displacement, %s. Processing code written with AI assistance (Claude, Anthropic).\nconst TS = " % cfg["place"] + json.dumps(out, allow_nan=False)
    open("ts-%s.js" % sid, "w").write(txt)
    rr, cc = np.where(area)
    print(sid, "| area cells", int(area.sum()), "| bbox S,W,N,E: %.4f,%.4f,%.4f,%.4f" % (n - (rr.max() + 1) * step, w + cc.min() * step, n - rr.min() * step, w + (cc.max() + 1) * step), "| bytes", len(txt))
    for t in out["tracks"]:
        ok = [p for p in t["pairs"] if p["mask_good_pct"] >= 50]
        tot = [v for i, v in enumerate(t["total_mm"]) if v is not None and area.ravel()[i]]
        print("   ", t["id"], "pairs used", len(ok), "of", len(t["pairs"]), "| per pair cm:", min(p["med_cm"] for p in ok), "to", max(p["med_cm"] for p in ok),
              "| running total cm:", round(sum(p["med_cm"] for p in ok), 1), "| pairs moving away:", sum(1 for p in ok if p["med_cm"] < 0),
              "| total-map cm in area (min, median):", round(min(tot) / 10, 1), round(float(np.median(tot)) / 10, 1))

for sid, cfg in PLACES.items():
    export(sid, cfg, SERIES[sid])

# Step 10. What this does: downloads the data files for the web map.
files.download("ts-pv.js")
files.download("ts-mexico.js")

# Step 11. What this does: saves every step above into one readable file, nisar_processing.py, and downloads it.
import re
steps = {}
for src in In:                                   # In holds the code of every cell that was run
    m = re.match(r"# (NISAR ground movement|Step (\d+))", src)
    if m: steps[int(m.group(2) or 1)] = src      # keep the latest version of each step
header = ("# These are the cells run in Google Colab on 2026-11-14 for NASA Space Apps 2026, in order.\n"
          "# AI acknowledgment: this code was written with Claude (Anthropic), an AI assistant.\n"
          "# Paste each step into its own Colab cell. get_ipython().system(...) lines are Colab shell commands.\n\n")
open("nisar_processing.py", "w").write(header + "\n\n".join(steps[k] for k in sorted(steps)))
print("steps saved:", sorted(steps))
files.download("nisar_processing.py")