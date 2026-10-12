"""HTTP API for on-demand MVI analysis. Jobs are asynchronous: POST /jobs, then poll GET /jobs/{id}."""
import datetime as dt
import json
import os
import threading
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

import mvi_core as core

app = FastAPI(title='Akwa Ibom MVI API')
origins = [o.strip() for o in os.getenv('ALLOWED_ORIGINS', '*').split(',') if o.strip()]
app.add_middleware(CORSMiddleware, allow_origins=origins, allow_methods=['GET', 'POST'], allow_headers=['*'])

JOBS: Dict[str, dict] = {}
CACHE: Dict[str, tuple] = {}          # key -> (timestamp, result)
RATE: Dict[str, list] = {}
POOL = ThreadPoolExecutor(max_workers=int(os.getenv('MAX_CONCURRENT_JOBS', '2')))
CACHE_TTL = 50 * 60                   # download links expire, so don't serve them for long
MAX_PER_HOUR = int(os.getenv('MAX_JOBS_PER_IP_PER_HOUR', '8'))
STATE = {'ready': False, 'error': None}


def _warm():
    try:
        core.init()
        core.static()
        STATE['ready'] = True
    except Exception as e:  # surfaced via /meta
        traceback.print_exc()
        STATE['error'] = str(e)[:300]


@app.on_event('startup')
def startup():
    threading.Thread(target=_warm, daemon=True).start()


class JobReq(BaseModel):
    start: str
    end: str
    lga: str = 'ALL'
    weights: Optional[Dict[str, float]] = None


@app.get('/health')
def health():
    return {'ok': True, 'ready': STATE['ready'], 'error': STATE['error']}


@app.get('/meta')
def meta():
    if STATE['error']:
        raise HTTPException(500, 'Earth Engine initialisation failed: ' + STATE['error'])
    if not STATE['ready']:
        raise HTTPException(503, 'Service is warming up, retry in about a minute.')
    lo, hi = core.date_bounds()
    return {'lgas': core.static()['names'], 'min_date': lo.isoformat(), 'max_date': hi.isoformat(),
            'min_days': 30, 'max_days': 1100, 'criteria': core.KEYS}


def _validate(req: JobReq):
    try:
        s, e = dt.date.fromisoformat(req.start), dt.date.fromisoformat(req.end)
    except ValueError:
        raise HTTPException(400, 'Dates must be YYYY-MM-DD.')
    lo, hi = core.date_bounds()
    if s < lo or e > hi:
        raise HTTPException(400, 'Dates must be between %s and %s.' % (lo, hi))
    if not 30 <= (e - s).days + 1 <= 1100:
        raise HTTPException(400, 'Period must be between 30 days and about 3 years.')
    if req.lga != 'ALL' and req.lga not in core.static()['names']:
        raise HTTPException(400, 'Unknown LGA.')
    w = {}
    for k, v in (req.weights or {}).items():
        if k not in core.KEYS or not (0 <= v <= 1):
            raise HTTPException(400, 'Bad weight: ' + k)
        w[k] = round(float(v), 4)
    if w and sum(w.values()) <= 0:
        raise HTTPException(400, 'Weights must not all be zero.')
    return w


def _ip(request: Request):
    fwd = request.headers.get('x-forwarded-for', '')
    return (fwd.split(',')[0].strip() or (request.client.host if request.client else 'x'))


@app.post('/jobs')
def create_job(req: JobReq, request: Request):
    if not STATE['ready']:
        raise HTTPException(503, 'Service is warming up, retry in about a minute.')
    w = _validate(req)
    key = json.dumps([req.start, req.end, req.lga, w], sort_keys=True)
    now = time.time()
    for k in [k for k, (t, _) in CACHE.items() if now - t > CACHE_TTL]:
        CACHE.pop(k, None)
    for k in [k for k, j in JOBS.items() if now - j['created'] > 2 * 3600]:
        JOBS.pop(k, None)
    jid = uuid.uuid4().hex[:12]
    if key in CACHE:
        JOBS[jid] = dict(status='done', message='Done (cached)', result=CACHE[key][1], created=now)
        return {'id': jid}
    ip = _ip(request)
    hits = [t for t in RATE.get(ip, []) if now - t < 3600]
    if len(hits) >= MAX_PER_HOUR:
        raise HTTPException(429, 'Rate limit reached (%d analyses per hour). Please try again later.' % MAX_PER_HOUR)
    hits.append(now)
    RATE[ip] = hits
    JOBS[jid] = dict(status='queued', message='Queued…', created=now)
    POOL.submit(_run, jid, req, w, key)
    return {'id': jid}


def _run(jid, req: JobReq, w, key):
    job = JOBS[jid]

    def report(m):
        job['message'] = m
    try:
        job['status'] = 'running'
        job['message'] = 'Starting…'
        res = core.run_analysis(req.start, req.end, req.lga, w, report)
        CACHE[key] = (time.time(), res)
        job.update(status='done', message='Done', result=res)
    except Exception as e:
        traceback.print_exc()
        job.update(status='error', error=str(e)[:300] or 'Analysis failed.')


@app.get('/jobs/{jid}')
def get_job(jid: str):
    j = JOBS.get(jid)
    if not j:
        raise HTTPException(404, 'Unknown or expired job.')
    return {'status': j['status'], 'message': j.get('message'), 'error': j.get('error'),
            'result': j.get('result')}
