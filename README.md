# Akwa Ibom Malaria Vulnerability Dashboard

Static site + an Earth Engine build job. Visitors never need Earth Engine: the job pulls everything
through the Earth Engine API and writes two small static files that the page loads.

```
weights.json ─┐
config.json  ─┼─► scripts/build_data.py (Earth Engine API) ─► data/akwa_ibom_mvi.json
cases_2024.csv┘                                              ─► data/akwa_ibom_lga.geojson ─► index.html
```

## One-time setup
1. **Earth Engine project.** Use a Google Cloud project registered for Earth Engine (non-commercial registration
   is fine for a public research dashboard). Enable the *Earth Engine API* on it.
2. **Service account.** IAM & Admin → Service accounts → create one → grant roles
   *Earth Engine Resource Viewer* and *Service Usage Consumer*. Register the service account for Earth Engine
   if your project requires it. Create a JSON key.
3. **GitHub repo.** Push this folder. Add repository secrets:
   - `EE_PROJECT` = your Cloud project id
   - `EE_SERVICE_ACCOUNT_JSON` = the full contents of the key file
4. **Netlify.** New site from the repo, no build command, publish directory `.` (already in `netlify.toml`).
5. Run **Actions → Refresh Earth Engine data → Run workflow**. It commits the data files and Netlify redeploys.

Never commit the key file. The browser only ever sees the static JSON.

## Day-to-day
- **Change weights:** on the dashboard's *AHP weights* tab click *Download weights.json*, replace `weights.json`,
  commit. The workflow reruns automatically and class populations are recomputed.
- **Add DHIS2 cases:** fill `data/cases_2024.csv` (`lga,cases`, names matched ignoring case/punctuation), commit.
- **Roads / health facilities:** upload OSM roads and GRID3 facilities to Earth Engine assets and put the asset ids in
  `config.json` (`roads_asset`, `facilities_asset`). Until then those two criteria are dropped and their weight is
  redistributed, and the dashboard's access quadrant is hidden.
- **GRID3 LGA boundaries:** set `lga_asset` and `lga_name_field` in `config.json`. The default is FAO GAUL 2015,
  and the build prints a warning if it doesn't find 31 LGAs.

## Run locally
```
pip install earthengine-api
earthengine authenticate
EE_PROJECT=your-project python scripts/build_data.py
python -m http.server 8000      # open http://localhost:8000
```

## Why not live per-request?
Each run takes minutes and the inputs are annual, so every visitor would get the same answer. A live endpoint
would add latency, cost, quota exposure and a public credential for no benefit. The per-visitor interactivity
(changing weights, ranks, sensitivity) already runs in the browser, because a weighted linear combination of LGA
mean scores equals the LGA mean of the combination. Only the pixel classes and class populations depend on the
weights used at build time, and the dashboard shows a notice when the visitor's weights differ from the build's.
