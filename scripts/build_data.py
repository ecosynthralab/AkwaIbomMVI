#!/usr/bin/env python3
"""
Build the static data files for the Akwa Ibom Malaria Vulnerability Dashboard.

Runs the Earth Engine pipeline through the Python API (no copy-paste into the Code Editor),
then writes:
  data/akwa_ibom_mvi.json      per-LGA mean criterion scores, class areas, class populations, cases
  data/akwa_ibom_lga.geojson   simplified LGA polygons for the map

The public site only ever reads these two static files, so visitors never touch Earth Engine.
Scores are 0-1 where 1 = MOST vulnerable. Method adapted from Abbah et al. (2026), Discover Public Health 23:551.

Auth (either):
  * CI:    EE_SERVICE_ACCOUNT_JSON (full key JSON as a string) + EE_PROJECT
  * local: `earthengine authenticate`, then set EE_PROJECT
"""
import csv, datetime, json, os, pathlib, re, sys
import ee

ROOT = pathlib.Path(__file__).resolve().parent.parent
KEYS = ['stream', 'wetland', 'precip', 'twi', 'lst', 'cropland', 'elev', 'slope', 'ndvi',
        'health', 'road', 'builtup', 'bare']


def init():
    project = os.environ.get('EE_PROJECT')
    if not project:
        sys.exit('Set EE_PROJECT to your Earth Engine-enabled Google Cloud project id.')
    key = os.environ.get('EE_SERVICE_ACCOUNT_JSON')
    if key:
        email = json.loads(key)['client_email']
        ee.Initialize(ee.ServiceAccountCredentials(email, key_data=key), project=project)
    else:
        ee.Initialize(project=project)


def norm(s):
    return re.sub(r'[^a-z0-9]', '', str(s or '').lower())


def main():
    cfg = json.loads((ROOT / 'config.json').read_text())
    weights = json.loads((ROOT / 'weights.json').read_text())['weights']
    init()

    YEAR, SCALE, CRS, CAP = cfg['year'], cfg['scale_m'], cfg['crs'], cfg['distance_cap_m']
    name_field = cfg.get('lga_name_field') or 'ADM2_NAME'
    d0, d1 = f'{YEAR}-01-01', f'{YEAR + 1}-01-01'

    # ---- study area
    if cfg.get('lga_asset'):
        lgas = ee.FeatureCollection(cfg['lga_asset'])
    else:
        lgas = (ee.FeatureCollection('FAO/GAUL/2015/level2')
                .filter(ee.Filter.eq('ADM0_NAME', 'Nigeria'))
                .filter(ee.Filter.eq('ADM1_NAME', 'Akwa Ibom')))
        name_field = 'ADM2_NAME'
    n_lga = lgas.size().getInfo()
    print('LGAs found:', n_lga)
    if n_lga != 31:
        print('WARNING: expected 31 LGAs. Check boundaries/name filters.', file=sys.stderr)
    aoi = lgas.geometry()

    # ---- optional vector assets (roads, facilities); weights are redistributed if absent
    roads = ee.FeatureCollection(cfg['roads_asset']) if cfg.get('roads_asset') else None
    facil = ee.FeatureCollection(cfg['facilities_asset']) if cfg.get('facilities_asset') else None
    W = dict(weights)
    if roads is None:
        W['road'] = 0.0
    if facil is None:
        W['health'] = 0.0
    tot = sum(W.values())
    W = {k: v / tot for k, v in W.items()}
    print('Effective weights:', {k: round(v, 4) for k, v in W.items()})

    proj = ee.Projection(CRS).atScale(SCALE)

    # ---- land cover (Dynamic World 10 m, annual mode)
    dw_col = (ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filterDate(d0, d1)
              .filterBounds(aoi).select('label'))
    dw = dw_col.mode().setDefaultProjection(ee.Image(dw_col.first()).projection())
    water, flooded, crops, built, bare = dw.eq(0), dw.eq(3), dw.eq(4), dw.eq(6), dw.eq(7)
    land = water.Not()

    def present(mask):
        return (mask.reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024)
                .reproject(proj).gt(0.02))

    def dist(pres):
        return (pres.selfMask().distance(ee.Kernel.euclidean(CAP, 'meters'))
                .unmask(CAP).min(CAP))

    # ---- hydrology and terrain
    merit = ee.Image('MERIT/Hydro/v1_0_1')
    upa = merit.select('upa')
    d_stream = dist(present(upa.gt(2)).Or(present(water)))
    d_wet = dist(present(flooded))
    d_crop = dist(present(crops))
    d_built = dist(present(built))
    d_bare = dist(present(bare))

    srtm = ee.Image('USGS/SRTMGL1_003')
    elev = srtm.select('elevation').rename('elev')
    slope = ee.Terrain.slope(srtm).rename('slope')
    slope_m = ee.Terrain.slope(merit.select('elv')).multiply(3.141592653589793 / 180).max(0.0017)
    twi = upa.multiply(1e6).divide(90).divide(slope_m.tan()).log().rename('twi')

    # ---- climate and vegetation
    precip = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate(d0, d1).sum().rename('precip')
    lst_c = (ee.ImageCollection('MODIS/061/MOD11A2').filterDate(d0, d1).select('LST_Day_1km')
             .mean().multiply(0.02).subtract(273.15))
    ndvi = (ee.ImageCollection('MODIS/061/MOD13Q1').filterDate(d0, d1).select('NDVI')
            .mean().multiply(0.0001).rename('ndvi'))

    def pnorm(img, name):
        """2nd-98th percentile min-max, with the limits fetched once as plain floats."""
        p = img.reduceRegion(reducer=ee.Reducer.percentile([2, 98]), geometry=aoi, scale=250,
                             maxPixels=int(1e10), bestEffort=True, tileScale=4).getInfo()
        lo, hi = p[name + '_p2'], p[name + '_p98']
        print(f'  {name}: p2={lo:.4g} p98={hi:.4g}')
        return img.subtract(lo).divide(hi - lo).clamp(0, 1)

    print('Normalising layers...')
    s_health = (facil.distance(searchRadius=20000, maxError=50).unmask(20000).min(20000).divide(20000)
                if facil is not None else ee.Image(0))
    s_road = (roads.distance(searchRadius=10000, maxError=50).unmask(10000).min(10000).divide(10000)
              if roads is not None else ee.Image(0))

    S = ee.Image.cat([
        ee.Image(1).subtract(d_stream.divide(CAP)).rename('s_stream'),
        ee.Image(1).subtract(d_wet.divide(CAP)).rename('s_wetland'),
        pnorm(precip, 'precip').rename('s_precip'),
        pnorm(twi, 'twi').rename('s_twi'),
        # temperature suitability peaks near 25 C (Mordecai et al. 2013): 0 at about 16 and 34 C
        ee.Image(1).subtract(lst_c.subtract(25).abs().divide(9)).clamp(0, 1).rename('s_lst'),
        ee.Image(1).subtract(d_crop.divide(CAP)).rename('s_cropland'),
        ee.Image(1).subtract(pnorm(elev, 'elev')).rename('s_elev'),
        ee.Image(1).subtract(pnorm(slope, 'slope')).rename('s_slope'),
        pnorm(ndvi, 'ndvi').rename('s_ndvi'),
        s_health.rename('s_health'),
        s_road.rename('s_road'),
        d_built.divide(CAP).rename('s_builtup'),
        d_bare.divide(CAP).rename('s_bare'),
    ])

    # ---- weighted linear combination and classes
    mvi = ee.Image(0)
    for k in KEYS:
        mvi = mvi.add(S.select('s_' + k).multiply(W[k]))
    mvi = mvi.updateMask(land).rename('mvi')
    valid = mvi.mask()
    q = mvi.reduceRegion(reducer=ee.Reducer.percentile([2, 98]), geometry=aoi, scale=250,
                         maxPixels=int(1e10), bestEffort=True, tileScale=4).getInfo()
    lo, hi = q['mvi_p2'], q['mvi_p98']
    step = (hi - lo) / 3
    breaks = [lo + step, hi - step]
    print('Class breaks:', breaks)
    cls = (ee.Image(2).where(mvi.lt(breaks[0]), 1).where(mvi.gte(breaks[1]), 3)
           .updateMask(valid).rename('cls'))  # 1 low, 2 medium, 3 high

    # ---- population (WorldPop 2020) and area per class
    pop = (ee.ImageCollection('WorldPop/GP/100m/pop')
           .filter(ee.Filter.eq('country', 'NGA')).filter(ee.Filter.eq('year', 2020))
           .mosaic().select('population'))
    area = ee.Image.pixelArea().divide(1e6)
    sum_img = S.updateMask(valid).addBands([
        ee.Image(1).updateMask(valid).rename('n'),
        pop.updateMask(valid).rename('pop_total'),
        pop.updateMask(cls.eq(3)).rename('pop_high'),
        pop.updateMask(cls.eq(2)).rename('pop_med'),
        pop.updateMask(cls.eq(1)).rename('pop_low'),
        area.updateMask(cls.eq(3)).rename('area_high_km2'),
        area.updateMask(cls.eq(2)).rename('area_med_km2'),
        area.updateMask(cls.eq(1)).rename('area_low_km2'),
    ])

    # ---- zonal sums in chunks (keeps each request small); means computed here
    print('Zonal statistics...')
    lst_feats = lgas.toList(lgas.size())
    CH = 6
    rows, geo_feats = [], []
    for off in range(0, n_lga, CH):
        part = ee.FeatureCollection(lst_feats.slice(off, off + CH))
        res = sum_img.reduceRegions(collection=part, reducer=ee.Reducer.sum(), scale=SCALE,
                                    crs=CRS, tileScale=8)
        res = res.map(lambda f: ee.Feature(ee.Feature(f).geometry().simplify(100), ee.Feature(f).toDictionary()))
        feats = res.getInfo()['features']
        for f in feats:
            p = f['properties']
            n = p.get('n') or 0
            name = p.get(name_field) or p.get('ADM2_NAME') or p.get('name')
            rows.append({
                'name': name,
                's': [round((p.get('s_' + k) or 0) / n, 5) if n else 0 for k in KEYS],
                'pop_total': round(p.get('pop_total') or 0),
                'ph': round(p.get('pop_high') or 0), 'pm': round(p.get('pop_med') or 0),
                'pl': round(p.get('pop_low') or 0),
                'ah': round(p.get('area_high_km2') or 0, 2), 'am': round(p.get('area_med_km2') or 0, 2),
                'al': round(p.get('area_low_km2') or 0, 2),
                'cases': None,
            })
            geo_feats.append({'type': 'Feature', 'properties': {'lga': name},
                              'geometry': round_geom(f['geometry'])})
        print(f'  {min(off + CH, n_lga)}/{n_lga}')

    # ---- optional DHIS2 cases (data/cases_2024.csv: lga,cases)
    cp = ROOT / cfg.get('cases_csv', 'data/cases_2024.csv')
    if cp.exists():
        cases = {}
        for r in csv.DictReader(cp.open(encoding='utf-8-sig')):
            try:
                cases[norm(r['lga'])] = float(r['cases'])
            except (KeyError, ValueError):
                pass
        miss = []
        for r in rows:
            r['cases'] = cases.get(norm(r['name']))
        if cases:
            matched = sum(r['cases'] is not None for r in rows)
            print(f'Cases matched for {matched}/{len(rows)} LGAs')

    out = {
        'generated': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'year': YEAR, 'scale_m': SCALE, 'crs': CRS,
        'criteria': KEYS, 'weights': {k: round(v, 5) for k, v in W.items()},
        'class_breaks': breaks, 'rows': rows,
    }
    (ROOT / 'data').mkdir(exist_ok=True)
    (ROOT / 'data' / 'akwa_ibom_mvi.json').write_text(json.dumps(out, separators=(',', ':')))
    (ROOT / 'data' / 'akwa_ibom_lga.geojson').write_text(
        json.dumps({'type': 'FeatureCollection', 'features': geo_feats}, separators=(',', ':')))
    print('Wrote data/akwa_ibom_mvi.json and data/akwa_ibom_lga.geojson')


def round_geom(g, nd=4):
    """Polygon/MultiPolygon with rounded coordinates; flattens GeometryCollections."""
    def rc(c):
        return [round(c[0], nd), round(c[1], nd)] if isinstance(c[0], (int, float)) else [rc(x) for x in c]

    polys = []

    def walk(geom):
        if not geom:
            return
        t = geom.get('type')
        if t == 'Polygon':
            polys.append(rc(geom['coordinates']))
        elif t == 'MultiPolygon':
            polys.extend(rc(geom['coordinates']))
        elif t == 'GeometryCollection':
            for x in geom.get('geometries', []):
                walk(x)

    walk(g)
    if len(polys) == 1:
        return {'type': 'Polygon', 'coordinates': polys[0]}
    return {'type': 'MultiPolygon', 'coordinates': polys}


if __name__ == '__main__':
    main()
