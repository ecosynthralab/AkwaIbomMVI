# Akwa Ibom MVI: Earth Engine App version

app/akwa_ibom_mvi_app.js   paste into the Earth Engine Code Editor and publish as an App
site/                      index.html + site-config.js for your Netlify repo

1. Code Editor (code.earthengine.google.com): New script, paste the file, click Run. Check the panel appears.
2. Click Apps (top right of the script panel) -> New App. Choose your Cloud project, source = current script,
   permissions = Anyone can view. Publish. Copy the app URL.
3. Put the URL in site/site-config.js as window.MVI_APP, commit and push site/ files. Netlify redeploys.

Edit the weights in the script (DEFAULT_W) or via the sliders in the app. Add ROADS / FACILITIES assets in the script
configuration to include the access criteria.

## Suitability and risk dashboard (index.html)

`index.html` is the planner-facing site. It reads static files from `data/`:

| File | Built by | Needs Earth Engine |
|---|---|---|
| `akwa_ibom_burden.json` | `python scripts/build_burden.py` (2023 cases CSV + GRID3 facilities) | no |
| `akwa_ibom_suitability.json` | `python scripts/build_data.py` (5-layer equal-weight suitability, population, facility access) | yes (`EE_PROJECT`, optional `EE_SERVICE_ACCOUNT_JSON`) |

Settings (case year, equal weights, access threshold) are in `config.json`. The *Refresh Earth Engine data* workflow
runs both scripts. Until `akwa_ibom_suitability.json` exists the site shows burden and facility counts and a notice.
The older 13-criterion Earth Engine App above is separate from this dashboard.
