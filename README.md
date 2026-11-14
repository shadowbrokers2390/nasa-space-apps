# Ground Movement Near Infrastructure

An interactive map of ground movement under roads, buildings, and rail lines, measured by the
NASA-ISRO NISAR radar satellite. Built for the 2026 NASA Space Apps Challenge, "Dancing with the
SARs", during the event on November 14-15, 2026.

**Live map:** https://shadowbrokers2390.github.io/nasa-space-apps/
(Mexico City: https://shadowbrokers2390.github.io/nasa-space-apps/?site=mexico)

## Two places, two kinds of surface change

| Place | Change | What the data shows (line of sight, "at least") |
|---|---|---|
| Portuguese Bend, California | Landslide | About 2 km² moved away from the satellite in all 25 ascending date pairs, about 57 cm in total from Nov 2025 to Oct 2026. 222 mapped buildings touch it. |
| Mexico City | Sinking ground | About 73 km² east of the city center moved away by about 16 cm from Oct 2025 to Oct 2026. The two viewing directions agree (about 16 cm and 15 cm), which is what vertical sinking should look like. |

## What the app does

- Shows line-of-sight ground displacement for each date pair, with a date slider and a Play button.
- Shows the total of all date pairs as one map.
- Charts the moving area's running total over time.
- Compares two satellite viewing directions.
- Outlines the moving area and the infrastructure that crosses it.
- States the limits of the measurement and what an engineer would check on the ground.

The NISAR data is **Provisional**. Radar measures motion only toward or away from the satellite,
and fast-moving ground can be underestimated, so every number is "at least". See the Limits panel
in the app.

## What I tried that did not work

I also processed California's Central Valley near Corcoran, where farmland is known to sink. In
this data no sinking patch was consistent from one date pair to the next, and the two viewing
directions did not agree, so I left that site out instead of showing a result I could not trust.

## Files

| File | What it is |
|---|---|
| `index.html` | The web map (Leaflet). Add `?site=mexico` to the address for Mexico City. |
| `ts-pv.js`, `ts-mexico.js` | Displacement data exported from Google Colab, one file per place. |
| `near-pv.js`, `near-mexico.js` | OpenStreetMap buildings, roads, rail lines, and runways near each moving area. |
| `nisar_processing.py` | Every step run in Google Colab, in order. |

## Data and credits

- NISAR L2 GUNW, collection NISAR_L2_GUNW_PROVISIONAL_V1, via the Alaska Satellite Facility DAAC
  and NASA Earthdata. Portuguese Bend: track 34 frame 19 (ascending), track 13 frame 71
  (descending). Mexico City: track 19 frame 11 (ascending), track 113 frame 79 (descending).
- Buildings, roads, rail lines, and runways: © OpenStreetMap contributors.
- Base map: Esri World Street Map. Map library: Leaflet.

## Use of AI

The code in this repository (the web page and the Python processing) and the text in the app were
written with Claude (Anthropic), an AI assistant, working under my direction. Each code and data
file says so in a comment at the top. The project idea and the choice of the first site are mine.
The displacement values come directly from NASA's NISAR data, not from an AI model. No AI-generated
images, video, or audio are used.

Before the event I explored NISAR data for Portuguese Bend to learn how it works. The project in
this repository was built during the event.
