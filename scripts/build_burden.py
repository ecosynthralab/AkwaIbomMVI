#!/usr/bin/env python3
"""
Offline build (no Earth Engine) of observed malaria burden and health-facility counts per LGA.

Inputs  (data/):
  akwa_ibom_lga_malaria_2023.csv      2023 Ministry of Health facility-reported cases, keyed by adm2_pcode
  akwa_ibom_lga_malaria_2023.geojson  LGA polygons with the same pcodes
  grid3_health_facilities.geojson     GRID3 v2.0 health facilities (NHFR 2024)
Output:
  data/akwa_ibom_burden.json

Burden metric: RDT test positivity rate (TPR) = RDT positives / fever patients tested by RDT.
Age groups are combined by summing counts, never by averaging rates.
"""
import csv, datetime, json, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
LEVELS = {'Primary': 0, 'Secondary': 1, 'Tertiary': 2}


def num(v):
    v = (v or '').replace(',', '').strip()
    return int(float(v)) if v else None


def ring_contains(ring, x, y):
    inside = False
    j = len(ring) - 1
    for i in range(len(ring)):
        xi, yi = ring[i][0], ring[i][1]
        xj, yj = ring[j][0], ring[j][1]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / (yj - yi) + xi:
            inside = not inside
        j = i
    return inside


def polygon_contains(poly, x, y):
    return ring_contains(poly[0], x, y) and not any(ring_contains(h, x, y) for h in poly[1:])


def geom_contains(geom, x, y):
    polys = [geom['coordinates']] if geom['type'] == 'Polygon' else geom['coordinates']
    return any(polygon_contains(p, x, y) for p in polys)


def bbox(geom):
    polys = [geom['coordinates']] if geom['type'] == 'Polygon' else geom['coordinates']
    pts = [c for p in polys for c in p[0]]
    xs, ys = [c[0] for c in pts], [c[1] for c in pts]
    return min(xs), min(ys), max(xs), max(ys)


def rate(pos, tested):
    if pos is None or not tested:
        return None
    return round(pos / tested, 4)


def add(*vals):
    got = [v for v in vals if v is not None]
    return sum(got) if got else None


def main():
    lga_geo = json.loads((DATA / 'akwa_ibom_lga_malaria_2023.geojson').read_text(encoding='utf-8'))
    zones = [(f['properties']['adm2_pcode'], f['geometry'], bbox(f['geometry'])) for f in lga_geo['features']]

    fac_geo = json.loads((DATA / 'grid3_health_facilities.geojson').read_text(encoding='utf-8'))
    fac = {code: {'total': 0, 'primary': 0, 'secondary': 0, 'tertiary': 0, 'public': 0, 'private': 0}
           for code, _, _ in zones}
    points, outside = [], 0
    for f in fac_geo['features']:
        x, y = f['geometry']['coordinates'][:2]
        p = f['properties']
        code = next((c for c, g, b in zones
                     if b[0] <= x <= b[2] and b[1] <= y <= b[3] and geom_contains(g, x, y)), None)
        if code is None:
            outside += 1
            continue
        c = fac[code]
        c['total'] += 1
        level = (p.get('facility_level') or '').lower()
        if level in c:
            c[level] += 1
        if (p.get('ownership') or '').lower() == 'public':
            c['public'] += 1
        elif (p.get('ownership') or '').lower() == 'private':
            c['private'] += 1
        points.append([round(x, 4), round(y, 4), LEVELS.get(p.get('facility_level'), 3)])

    rows = []
    with (DATA / 'akwa_ibom_lga_malaria_2023.csv').open(encoding='utf-8-sig') as fh:
        for r in csv.DictReader(fh):
            pc = r['adm2_pcode']
            t5, to, p5, po = (num(r['fever_tested_rdt_u5']), num(r['fever_tested_rdt_o5']),
                              num(r['positive_rdt_u5']), num(r['positive_rdt_o5']))
            u5, uo = num(r['uncomplicated_malaria_u5']), num(r['uncomplicated_malaria_o5'])
            tested, pos, unc = add(t5, to), add(p5, po), add(u5, uo)
            area = float(r['area_sqkm'])
            flags = []
            if (p5 or 0) > (t5 or 0) or (po or 0) > (to or 0):
                flags.append('positives_exceed_tested')
            if pos and unc and abs(pos - unc) / max(pos, unc) > 0.2:
                flags.append('rdt_vs_uncomplicated_gap_20pct')
            rows.append({
                'pcode': pc, 'name': r['lga_name'], 'area_km2': round(area, 1),
                'tested_u5': t5, 'tested_o5': to, 'pos_u5': p5, 'pos_o5': po,
                'tested': tested, 'pos': pos,
                'tpr': rate(pos, tested), 'tpr_u5': rate(p5, t5), 'tpr_o5': rate(po, to),
                'unc_u5': u5, 'unc_o5': uo, 'unc': unc,
                'flags': flags,
                'fac': fac.get(pc, {}),
                'fac_per_100km2': round(fac.get(pc, {}).get('total', 0) / area * 100, 2),
            })

    missing = [r['name'] for r in rows if r['pcode'] not in fac]
    if missing:
        raise SystemExit(f'pcodes in cases CSV missing from the GeoJSON: {missing}')

    out = {
        'generated': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'year': 2023,
        'sources': {'cases': 'Ministry of Health, Akwa Ibom State, 2023 (facility-reported)',
                    'facilities': 'GRID3 Nigeria health facilities v2.0 (NHFR 2024)'},
        'facilities_outside_lgas': outside,
        'rows': rows,
        'facilities': points,
    }
    (DATA / 'akwa_ibom_burden.json').write_text(json.dumps(out, separators=(',', ':')))
    n_fac = sum(r['fac']['total'] for r in rows)
    print(f'{len(rows)} LGAs, {n_fac} facilities inside LGAs, {outside} outside')
    for r in rows:
        if r['flags']:
            print(' flag', r['name'], r['flags'])


if __name__ == '__main__':
    main()
