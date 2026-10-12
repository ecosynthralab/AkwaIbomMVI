"""
On-demand Earth Engine pipeline for the Akwa Ibom Malaria Vulnerability Index (MVI).
Method adapted from Abbah et al. (2026), Discover Public Health 23:551.

Scores are 0-1 where 1 = MOST vulnerable. Time-varying layers (rainfall, LST, NDVI) use the requested
period; land cover uses the 12 months ending at the period end; terrain and hydrology are static.
Normalisation limits and class breaks are the 2nd-98th percentile range across Akwa Ibom for that period,
so results describe RELATIVE vulnerability within the state for the chosen period.
"""
import datetime as dt
import json
import os
import pathlib
import re
import threading

import ee

HERE = pathlib.Path(__file__).resolve().parent
CFG = json.loads((HERE / 'config.json').read_text())
DEFAULT_W = json.loads((HERE / 'weights.json').read_text())['weights']
KEYS = ['stream', 'wetland', 'precip', 'twi', 'lst', 'cropland', 'elev', 'slope', 'ndvi',
        'health', 'road', 'builtup', 'bare']
SCALE, CRS, CAP = CFG['scale_m'], CFG['crs'], CFG['distance_cap_m']
MIN_DATE = dt.date(2016, 1, 1)      # Dynamic World starts mid-2015
LAG_DAYS = 45                       # leave time for CHIRPS / MODIS to be published
NODATA_MVI = -9999


def date_bounds():
    return MIN_DATE, dt.date.today() - dt.timedelta(days=LAG_DAYS)


def init():
    project = os.environ.get('EE_PROJECT')
    if not project:
        raise RuntimeError('EE_PROJECT is not set')
    key = os.environ.get('EE_SERVICE_ACCOUNT_JSON')
    if key:
        email = json.loads(key)['client_email']
        ee.Initialize(ee.ServiceAccountCredentials(email, key_data=key), project=project)
    else:
        ee.Initialize(project=project)


# --------------------------------------------------------------------------- static objects
_static = {}
_slock = threading.Lock()


def _limits(img, name, aoi, scale=250):
    p = img.reduceRegion(reducer=ee.Reducer.percentile([2, 98]), geometry=aoi, scale=scale,
                         maxPixels=int(1e10), bestEffort=True, tileScale=4).getInfo()
    lo, hi = p[name + '_p2'], p[name + '_p98']
    return lo, (hi if hi > lo else lo + 1e-9)


def static():
    """Time-invariant objects and normalisation limits, created once per process."""
    with _slock:
        if _static:
            return _static
        name_field = CFG.get('lga_name_field') or 'ADM2_NAME'
        if CFG.get('lga_asset'):
            lgas = ee.FeatureCollection(CFG['lga_asset'])
        else:
            lgas = (ee.FeatureCollection('FAO/GAUL/2015/level2')
                    .filter(ee.Filter.eq('ADM0_NAME', 'Nigeria'))
                    .filter(ee.Filter.eq('ADM1_NAME', 'Akwa Ibom')))
            name_field = 'ADM2_NAME'
        names = sorted(lgas.aggregate_array(name_field).getInfo())
        aoi = lgas.geometry()

        merit = ee.Image('MERIT/Hydro/v1_0_1')
        upa = merit.select('upa')
        srtm = ee.Image('USGS/SRTMGL1_003')
        elev = srtm.select('elevation').rename('elev')
        slope = ee.Terrain.slope(srtm).rename('slope')
        slope_m = ee.Terrain.slope(merit.select('elv')).multiply(3.141592653589793 / 180).max(0.0017)
        twi = upa.multiply(1e6).divide(90).divide(slope_m.tan()).log().rename('twi')

        lim = {'elev': _limits(elev, 'elev', aoi), 'slope': _limits(slope, 'slope', aoi),
               'twi': _limits(twi, 'twi', aoi)}
        roads = ee.FeatureCollection(CFG['roads_asset']) if CFG.get('roads_asset') else None
        facil = ee.FeatureCollection(CFG['facilities_asset']) if CFG.get('facilities_asset') else None
        pop = (ee.ImageCollection('WorldPop/GP/100m/pop')
               .filter(ee.Filter.eq('country', 'NGA')).filter(ee.Filter.eq('year', 2020))
               .mosaic().select('population'))
        _static.update(dict(lgas=lgas, names=names, name_field=name_field, aoi=aoi, upa=upa,
                            elev=elev, slope=slope, twi=twi, lim=lim, roads=roads, facil=facil, pop=pop))
        return _static


def effective_weights(w):
    st = static()
    W = {k: float((w or {}).get(k, DEFAULT_W[k])) for k in KEYS}
    if st['roads'] is None:
        W['road'] = 0.0
    if st['facil'] is None:
        W['health'] = 0.0
    tot = sum(W.values())
    return {k: v / tot for k, v in W.items()}


# --------------------------------------------------------------------------- period context
_ctx_cache = {}
_clock = threading.Lock()


def _nrm(img, lim):
    lo, hi = lim
    return img.subtract(lo).divide(hi - lo).clamp(0, 1)


def _layers(st, start, end):
    d_end = (dt.date.fromisoformat(end) + dt.timedelta(days=1)).isoformat()
    dw_start = (dt.date.fromisoformat(end) - dt.timedelta(days=364)).isoformat()
    aoi = st['aoi']
    proj = ee.Projection(CRS).atScale(SCALE)

    dw_col = (ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filterDate(dw_start, d_end)
              .filterBounds(aoi).select('label'))
    dw = dw_col.mode().setDefaultProjection(ee.Image(dw_col.first()).projection())
    water, flooded, crops, built, bare = dw.eq(0), dw.eq(3), dw.eq(4), dw.eq(6), dw.eq(7)

    def present(mask):
        return (mask.reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024)
                .reproject(proj).gt(0.02))

    def dist(pres):
        return pres.selfMask().distance(ee.Kernel.euclidean(CAP, 'meters')).unmask(CAP).min(CAP)

    precip = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate(start, d_end).sum().rename('precip')
    lst_c = (ee.ImageCollection('MODIS/061/MOD11A2').filterDate(start, d_end).select('LST_Day_1km')
             .mean().multiply(0.02).subtract(273.15))
    ndvi = (ee.ImageCollection('MODIS/061/MOD13Q1').filterDate(start, d_end).select('NDVI')
            .mean().multiply(0.0001).rename('ndvi'))
    return dict(
        land=water.Not(),
        d_stream=dist(present(st['upa'].gt(2)).Or(present(water))),
        d_wet=dist(present(flooded)), d_crop=dist(present(crops)),
        d_built=dist(present(built)), d_bare=dist(present(bare)),
        precip=precip, lst=lst_c, ndvi=ndvi)


def get_context(start, end, W, report=lambda m: None):
    key = start + '|' + end + '|' + json.dumps({k: round(v, 5) for k, v in W.items()}, sort_keys=True)
    with _clock:
        if key in _ctx_cache:
            return _ctx_cache[key]
        if len(_ctx_cache) > 12:
            _ctx_cache.pop(next(iter(_ctx_cache)))
    st = static()
    aoi = st['aoi']
    lay = _layers(st, start, end)
    lim = {}
    for nm in ('precip', 'ndvi'):
        report('Normalising ' + ('rainfall' if nm == 'precip' else 'NDVI') + ' for the period…')
        lim[nm] = _limits(lay[nm], nm, aoi, scale=500)

    s_health = (st['facil'].distance(searchRadius=20000, maxError=50).unmask(20000).min(20000).divide(20000)
                if st['facil'] is not None else ee.Image(0))
    s_road = (st['roads'].distance(searchRadius=10000, maxError=50).unmask(10000).min(10000).divide(10000)
              if st['roads'] is not None else ee.Image(0))
    S = ee.Image.cat([
        ee.Image(1).subtract(lay['d_stream'].divide(CAP)).rename('s_stream'),
        ee.Image(1).subtract(lay['d_wet'].divide(CAP)).rename('s_wetland'),
        _nrm(lay['precip'], lim['precip']).rename('s_precip'),
        _nrm(st['twi'], st['lim']['twi']).rename('s_twi'),
        # temperature suitability peaks near 25 C (Mordecai et al. 2013); zero near 16 and 34 C
        ee.Image(1).subtract(lay['lst'].subtract(25).abs().divide(9)).clamp(0, 1).rename('s_lst'),
        ee.Image(1).subtract(lay['d_crop'].divide(CAP)).rename('s_cropland'),
        ee.Image(1).subtract(_nrm(st['elev'], st['lim']['elev'])).rename('s_elev'),
        ee.Image(1).subtract(_nrm(st['slope'], st['lim']['slope'])).rename('s_slope'),
        _nrm(lay['ndvi'], lim['ndvi']).rename('s_ndvi'),
        s_health.rename('s_health'), s_road.rename('s_road'),
        lay['d_built'].divide(CAP).rename('s_builtup'),
        lay['d_bare'].divide(CAP).rename('s_bare'),
    ])
    mvi = ee.Image(0)
    for k in KEYS:
        mvi = mvi.add(S.select('s_' + k).multiply(W[k]))
    mvi = mvi.updateMask(lay['land']).rename('mvi')
    valid = mvi.mask()

    report('Computing class breaks across Akwa Ibom…')
    lo, hi = _limits(mvi, 'mvi', aoi, scale=500)
    step = (hi - lo) / 3
    breaks = [lo + step, hi - step]
    cls = (ee.Image(2).where(mvi.lt(breaks[0]), 1).where(mvi.gte(breaks[1]), 3)
           .updateMask(valid).rename('cls'))   # 1 low, 2 medium, 3 high
    ctx = dict(mvi=mvi, cls=cls, valid=valid, lo=lo, hi=hi, breaks=breaks, lim=lim)
    with _clock:
        _ctx_cache[key] = ctx
    return ctx


# --------------------------------------------------------------------------- analysis
def _safe(s):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)


def run_analysis(start, end, lga, weights, report=lambda m: None):
    st = static()
    W = effective_weights(weights)
    ctx = get_context(start, end, W, report)
    mvi, cls, valid = ctx['mvi'], ctx['cls'], ctx['valid']
    names = st['names'] if lga == 'ALL' else [lga]
    if lga != 'ALL' and lga not in st['names']:
        raise ValueError('Unknown LGA: ' + lga)
    sel = st['lgas'] if lga == 'ALL' else st['lgas'].filter(ee.Filter.eq(st['name_field'], lga))
    geom = sel.geometry()
    region = geom.bounds()

    area = ee.Image.pixelArea().divide(1e6)
    pop = st['pop']
    sum_img = mvi.rename('mvi').addBands([
        ee.Image(1).updateMask(valid).rename('n'),
        pop.updateMask(valid).rename('pop_total'),
        pop.updateMask(cls.eq(3)).rename('pop_high'),
        pop.updateMask(cls.eq(2)).rename('pop_med'),
        pop.updateMask(cls.eq(1)).rename('pop_low'),
        area.updateMask(cls.eq(3)).rename('area_high_km2'),
        area.updateMask(cls.eq(2)).rename('area_med_km2'),
        area.updateMask(cls.eq(1)).rename('area_low_km2'),
    ])
    rows = []
    CH = 6
    for off in range(0, len(names), CH):
        chunk = names[off:off + CH]
        report('Zonal statistics %d/%d…' % (min(off + CH, len(names)), len(names)))
        part = st['lgas'].filter(ee.Filter.inList(st['name_field'], chunk))
        res = sum_img.reduceRegions(collection=part, reducer=ee.Reducer.sum(), scale=SCALE,
                                    crs=CRS, tileScale=8)
        res = res.map(lambda f: ee.Feature(None, ee.Feature(f).toDictionary()))
        for f in res.getInfo()['features']:
            p = f['properties']
            n = p.get('n') or 0
            rows.append(dict(
                name=p.get(st['name_field']) or p.get('ADM2_NAME'),
                mean_mvi=round((p.get('mvi') or 0) / n, 5) if n else None,
                pop_total=round(p.get('pop_total') or 0), pop_high=round(p.get('pop_high') or 0),
                pop_med=round(p.get('pop_med') or 0), pop_low=round(p.get('pop_low') or 0),
                area_high_km2=round(p.get('area_high_km2') or 0, 2),
                area_med_km2=round(p.get('area_med_km2') or 0, 2),
                area_low_km2=round(p.get('area_low_km2') or 0, 2)))
    rows.sort(key=lambda r: r['name'] or '')

    out = dict(period=dict(start=start, end=end), lga=lga,
               weights={k: round(v, 4) for k, v in W.items()},
               value_range=[ctx['lo'], ctx['hi']], breaks=ctx['breaks'], rows=rows,
               downloads={}, previews={}, notes=[])
    if st['roads'] is None or st['facil'] is None:
        out['notes'].append('Road and/or health-facility layers are not configured, so those criteria are '
                            'omitted and their weight is redistributed.')

    base = _safe('AkwaIbom_%s_%s_%s' % (lga, start, end))
    report('Preparing GeoTIFF downloads…')
    try:
        out['downloads']['mvi'] = mvi.clip(geom).unmask(NODATA_MVI).toFloat().getDownloadURL(
            dict(name=base + '_MVI', scale=SCALE, crs=CRS, region=region, format='GEO_TIFF'))
        out['downloads']['class'] = cls.clip(geom).unmask(0).toByte().getDownloadURL(
            dict(name=base + '_class', scale=SCALE, crs=CRS, region=region, format='GEO_TIFF'))
    except Exception as e:  # keep stats and previews even if the raster is too large for a direct download
        out['notes'].append('GeoTIFF download could not be prepared for this extent (%s). '
                            'Try a single LGA or a shorter period.' % str(e)[:160])
    report('Rendering previews…')
    try:
        out['previews']['mvi'] = mvi.clip(geom).visualize(
            min=ctx['lo'], max=ctx['hi'], palette=['1a9850', 'fee08b', 'd73027']
        ).getThumbURL(dict(region=region, dimensions=700, format='png'))
        out['previews']['class'] = cls.clip(geom).visualize(
            min=1, max=3, palette=['1a9850', 'ffffbf', 'd73027']
        ).getThumbURL(dict(region=region, dimensions=700, format='png'))
    except Exception as e:
        out['notes'].append('Preview could not be rendered (%s).' % str(e)[:120])
    return out
