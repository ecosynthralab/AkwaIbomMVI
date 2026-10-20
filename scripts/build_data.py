#!/usr/bin/env python3
"""
Build the malaria habitat suitability and access data for the Akwa Ibom dashboard (Earth Engine).

Writes data/akwa_ibom_suitability.json: per-LGA mean layer scores, suitability class areas/populations,
raw layer means, and health-facility access (distance and population beyond a threshold).

Habitat suitability combines FIVE layers with equal weights (0.2 each; provisional, see config.json):
  rainfall      CHIRPS annual total            p2-p98 min-max, higher = more suitable
  temperature   MODIS LST day mean (deg C)     peaks at 25 C, zero at about 16 and 34 C (Mordecai et al. 2013)
  ndvi          MODIS MOD13Q1 annual mean      p2-p98 min-max, higher = more suitable
  surface_water JRC Global Surface Water       1 - distance to water (occurrence >= 10%) / cap
  twi           Topographic Wetness Index      p2-p98 min-max, higher = more suitable
Scores are 0-1 where 1 = most suitable for Anopheles habitat. Suitability is NOT disease risk: population and
health-facility access are reported separately and never enter the suitability score.

Zones come from data/akwa_ibom_lga_malaria_2023.geojson (adm2_pcode is the join key).
Facilities come from data/grid3_health_facilities.geojson.

Auth (either):
  * CI:    EE_SERVICE_ACCOUNT_JSON (full key JSON as a string) + EE_PROJECT
  * local: `earthengine authenticate`, then set EE_PROJECT
"""
import datetime, json, os, pathlib, sys
import ee

ROOT = pathlib.Path(__file__).resolve().parent.parent
KEYS = ['rainfall', 'temperature', 'ndvi', 'surface_water', 'twi']


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


def main():
    cfg = json.loads((ROOT / 'config.json').read_text())
    W = cfg['weights']
    if abs(sum(W.values()) - 1) > 1e-6 or set(W) != set(KEYS):
        sys.exit('config.json weights must cover exactly ' + ', '.join(KEYS) + ' and sum to 1.')
    init()

    YEAR, SCALE, CRS, CAP = cfg['year'], cfg['scale_m'], cfg['crs'], cfg['distance_cap_m']
    ACCESS_KM = cfg['access_threshold_km']
    d0, d1 = f'{YEAR}-01-01', f'{YEAR + 1}-01-01'

    # ---- zones and facilities from the repo (no assets to maintain)
    lga_geo = json.loads((ROOT / cfg['lga_geojson']).read_text(encoding='utf-8'))
    lgas = ee.FeatureCollection([
        ee.Feature(ee.Geometry(f['geometry'], None, False),
                   {'pcode': f['properties']['adm2_pcode'], 'name': f['properties']['lga_name']})
        for f in lga_geo['features']])
    n_lga = len(lga_geo['features'])
    if n_lga != 31:
        print('WARNING: expected 31 LGAs.', file=sys.stderr)
    aoi = lgas.geometry()

    fac_geo = json.loads((ROOT / cfg['facilities_geojson']).read_text(encoding='utf-8'))
    facil = ee.FeatureCollection([
        ee.Feature(ee.Geometry.Point(f['geometry']['coordinates'][:2]))
        for f in fac_geo['features']])
    print('Facilities:', len(fac_geo['features']))

    proj = ee.Projection(CRS).atScale(SCALE)

    # ---- surface water: JRC occurrence, long-term; permanent water bodies are masked out of the land
    gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0)
    land = gsw.lt(50)
    water_present = (gsw.gte(10).selfMask().reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024)
                     .reproject(proj).gt(0.02))
    d_water = (water_present.selfMask().distance(ee.Kernel.euclidean(CAP, 'meters'))
               .unmask(CAP).min(CAP))

    # ---- terrain / TWI
    merit = ee.Image('MERIT/Hydro/v1_0_1')
    upa = merit.select('upa')
    slope_m = ee.Terrain.slope(merit.select('elv')).multiply(3.141592653589793 / 180).max(0.0017)
    twi = upa.multiply(1e6).divide(90).divide(slope_m.tan()).log().rename('twi')

    # ---- climate and vegetation for the case year
    precip = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate(d0, d1).sum().rename('rainfall')
    lst_c = (ee.ImageCollection('MODIS/061/MOD11A2').filterDate(d0, d1).select('LST_Day_1km')
             .mean().multiply(0.02).subtract(273.15).rename('temperature'))
    ndvi = (ee.ImageCollection('MODIS/061/MOD13Q1').filterDate(d0, d1).select('NDVI')
            .mean().multiply(0.0001).rename('ndvi'))

    def pnorm(img, name):
        """2nd-98th percentile min-max, limits fetched once as plain floats."""
        p = img.updateMask(land).reduceRegion(reducer=ee.Reducer.percentile([2, 98]), geometry=aoi,
                                              scale=250, maxPixels=int(1e10), bestEffort=True,
                                              tileScale=4).getInfo()
        lo, hi = p[name + '_p2'], p[name + '_p98']
        print(f'  {name}: p2={lo:.4g} p98={hi:.4g}')
        return img.subtract(lo).divide(hi - lo).clamp(0, 1)

    print('Normalising layers...')
    S = ee.Image.cat([
        pnorm(precip, 'rainfall').rename('s_rainfall'),
        ee.Image(1).subtract(lst_c.subtract(25).abs().divide(9)).clamp(0, 1).rename('s_temperature'),
        pnorm(ndvi, 'ndvi').rename('s_ndvi'),
        ee.Image(1).subtract(d_water.divide(CAP)).rename('s_surface_water'),
        pnorm(twi, 'twi').rename('s_twi'),
    ]).updateMask(land)

    suit = ee.Image(0)
    for k in KEYS:
        suit = suit.add(S.select('s_' + k).multiply(W[k]))
    suit = suit.rename('suit').updateMask(land)
    valid = suit.mask()
    q = suit.reduceRegion(reducer=ee.Reducer.percentile([2, 98]), geometry=aoi, scale=250,
                          maxPixels=int(1e10), bestEffort=True, tileScale=4).getInfo()
    lo, hi = q['suit_p2'], q['suit_p98']
    step = (hi - lo) / 3
    breaks = [lo + step, hi - step]
    print('Class breaks:', breaks)
    cls = (ee.Image(2).where(suit.lt(breaks[0]), 1).where(suit.gte(breaks[1]), 3)
           .updateMask(valid).rename('cls'))  # 1 low, 2 moderate, 3 high suitability

    # ---- population (WorldPop 2020, the latest available) and facility access
    pop = (ee.ImageCollection('WorldPop/GP/100m/pop')
           .filter(ee.Filter.eq('country', 'NGA')).filter(ee.Filter.eq('year', 2020))
           .mosaic().select('population'))
    d_fac = facil.distance(searchRadius=60000, maxError=50).unmask(60000).min(60000)
    area = ee.Image.pixelArea().divide(1e6)
    sum_img = S.addBands([
        suit,
        ee.Image(1).updateMask(valid).rename('n'),
        precip.updateMask(valid), lst_c.updateMask(valid), ndvi.updateMask(valid),
        gsw.updateMask(valid).rename('water_occ'),
        twi.updateMask(valid),
        pop.updateMask(valid).rename('pop_total'),
        pop.updateMask(cls.eq(3)).rename('pop_high'),
        pop.updateMask(cls.eq(2)).rename('pop_med'),
        pop.updateMask(cls.eq(1)).rename('pop_low'),
        area.updateMask(cls.eq(3)).rename('area_high_km2'),
        area.updateMask(cls.eq(2)).rename('area_med_km2'),
        area.updateMask(cls.eq(1)).rename('area_low_km2'),
        d_fac.updateMask(valid).rename('d_fac'),
        pop.multiply(d_fac).updateMask(valid).rename('pop_x_dist'),
        pop.updateMask(d_fac.gt(ACCESS_KM * 1000)).updateMask(valid).rename('pop_far'),
        pop.updateMask(cls.eq(3)).updateMask(d_fac.gt(ACCESS_KM * 1000)).rename('pop_high_far'),
    ])

    print('Zonal statistics...')
    lst_feats = lgas.toList(lgas.size())
    CH = 6
    rows = []
    for off in range(0, n_lga, CH):
        part = ee.FeatureCollection(lst_feats.slice(off, off + CH))
        res = sum_img.reduceRegions(collection=part, reducer=ee.Reducer.sum(), scale=SCALE,
                                    crs=CRS, tileScale=8)
        for f in res.getInfo()['features']:
            p = f['properties']
            n = p.get('n') or 0
            pop_t = p.get('pop_total') or 0

            def mean(key, nd=4):
                return round((p.get(key) or 0) / n, nd) if n else None

            rows.append({
                'pcode': p['pcode'], 'name': p['name'],
                's': [mean('s_' + k, 5) for k in KEYS],
                'suit': mean('suit', 5),
                'raw': {'rainfall_mm': mean('rainfall', 1), 'temperature_c': mean('temperature', 2),
                        'ndvi': mean('ndvi', 3), 'water_occurrence_pct': mean('water_occ', 1),
                        'twi': mean('twi', 2)},
                'pop_total': round(pop_t),
                'ph': round(p.get('pop_high') or 0), 'pm': round(p.get('pop_med') or 0),
                'pl': round(p.get('pop_low') or 0),
                'ah': round(p.get('area_high_km2') or 0, 2), 'am': round(p.get('area_med_km2') or 0, 2),
                'al': round(p.get('area_low_km2') or 0, 2),
                'access': {'mean_dist_km': round(mean('d_fac', 1) / 1000, 2) if n else None,
                           'pop_weighted_dist_km': round((p.get('pop_x_dist') or 0) / pop_t / 1000, 2)
                           if pop_t else None,
                           'pop_beyond': round(p.get('pop_far') or 0),
                           'pop_high_beyond': round(p.get('pop_high_far') or 0)},
            })
        print(f'  {min(off + CH, n_lga)}/{n_lga}')

    out = {
        'generated': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'year': YEAR, 'scale_m': SCALE, 'crs': CRS, 'population_year': 2020,
        'layers': KEYS, 'weights': W, 'access_threshold_km': ACCESS_KM,
        'class_breaks': breaks, 'rows': rows,
    }
    (ROOT / 'data').mkdir(exist_ok=True)
    (ROOT / 'data' / 'akwa_ibom_suitability.json').write_text(json.dumps(out, separators=(',', ':')))
    print('Wrote data/akwa_ibom_suitability.json')


if __name__ == '__main__':
    main()
