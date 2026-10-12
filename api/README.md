---
title: Akwa Ibom MVI API
emoji: 🦟
sdk: docker
app_port: 7860
---
# Akwa Ibom MVI analysis API

Runs the malaria vulnerability model in Google Earth Engine on demand (area + period + weights) and returns
statistics, previews and GeoTIFF download links. The public site calls this service; visitors need no Earth Engine account.

Endpoints: `GET /health`, `GET /meta`, `POST /jobs`, `GET /jobs/{id}` (jobs are asynchronous; the page polls).

## Environment variables
| Name | Kind | Value |
|---|---|---|
| `EE_PROJECT` | secret | your Earth Engine Cloud project id |
| `EE_SERVICE_ACCOUNT_JSON` | secret | full contents of the service-account key JSON |
| `ALLOWED_ORIGINS` | variable | your site, e.g. `https://akwaibommvi.netlify.app` (comma-separate several) |
| `MAX_JOBS_PER_IP_PER_HOUR` | variable, optional | default 8 |
| `MAX_CONCURRENT_JOBS` | variable, optional | default 2 |

## Deploy on Hugging Face Spaces (free, no card)
1. huggingface.co → New Space → SDK **Docker** → Blank → **Public** (secrets stay private).
2. Upload every file in this folder to the Space (Files → Add file → Upload).
3. Settings → Variables and secrets → add the table above.
4. The Space builds (a few minutes). Its URL is `https://USER-SPACENAME.hf.space`. Open `/health`: it should say `"ready": true` after about a minute.
5. Put that URL in the site's `site-config.js` as `window.MVI_API` and commit.

Free Spaces sleep when idle and wake on the next visit (about a minute). For always-on, use Cloud Run:
`gcloud run deploy mvi-api --source . --region europe-west1 --allow-unauthenticated --set-env-vars EE_PROJECT=...,ALLOWED_ORIGINS=...`
and add `EE_SERVICE_ACCOUNT_JSON` as a secret.

## Notes
- Allowed period: 2016-01-01 to about 6 weeks ago; 30 days to about 3 years.
- Land cover uses the 12 months ending at the period end; population is WorldPop 2020.
- Normalisation and class breaks are relative to Akwa Ibom for the chosen period.
- Set `roads_asset` / `facilities_asset` in `config.json` to include the access criteria; weights are redistributed until then.
- Rate limiting is in memory per instance. For heavy public traffic add a CDN or auth in front.
