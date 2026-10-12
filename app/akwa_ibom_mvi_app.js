// =====================================================================================
// Akwa Ibom Malaria Vulnerability Index (MVI): Earth Engine App
// AHP-WLC environmental suitability + population exposure, adapted from
// Abbah et al. (2026) Discover Public Health 23:551.
// Visitors pick an LGA and a period, run the model, and download GeoTIFFs. No GEE account needed.
// Scores are 0-1, 1 = MOST vulnerable. Classes are relative to Akwa Ibom for the chosen period.
// =====================================================================================

// ---------------------------------------------------------------- configuration
var CRS = 'EPSG:32632';      // UTM 32N (Akwa Ibom is about 7-8.2 E)
var SCALE = 100;             // m
var CAP = 5000;              // distance cap (m) for habitat-proximity layers
var ROADS = null;            // e.g. ee.FeatureCollection('projects/YOUR_PROJECT/assets/akwaibom_roads_osm')
var FACILITIES = null;       // e.g. ee.FeatureCollection('projects/YOUR_PROJECT/assets/akwaibom_facilities_grid3')
var MIN_DATE = '2016-01-01'; // Dynamic World starts mid-2015
var LAG_DAYS = 45;           // CHIRPS / MODIS publication lag
var NODATA = -9999;

var KEYS = ['stream','wetland','precip','twi','lst','cropland','elev','slope','ndvi','health','road','builtup','bare'];
var LABELS = {stream:'Distance to stream', wetland:'Distance to wetland', precip:'Precipitation',
  twi:'Topographic wetness index', lst:'Land surface temperature', cropland:'Distance from cropland',
  elev:'Elevation', slope:'Slope', ndvi:'NDVI', health:'Distance to health facilities',
  road:'Distance to roads', builtup:'Distance from built-up', bare:'Distance from bare ground'};
// Default weights: AHP priority vector from the Bayelsa pairwise matrix. Replace with Akwa Ibom expert weights.
var DEFAULT_W = {stream:0.2218, wetland:0.2022, precip:0.1295, twi:0.0980, lst:0.0804, cropland:0.0711,
  elev:0.0480, slope:0.0311, ndvi:0.0311, health:0.0253, road:0.0253, builtup:0.0253, bare:0.0109};

// ---------------------------------------------------------------- study area and static layers
var lgas = ee.FeatureCollection('FAO/GAUL/2015/level2')
  .filter(ee.Filter.eq('ADM0_NAME', 'Nigeria')).filter(ee.Filter.eq('ADM1_NAME', 'Akwa Ibom'));
var aoi = lgas.geometry();
var proj = ee.Projection(CRS).atScale(SCALE);

var merit = ee.Image('MERIT/Hydro/v1_0_1');
var upa = merit.select('upa');
var srtm = ee.Image('USGS/SRTMGL1_003');
var elev = srtm.select('elevation').rename('elev');
var slope = ee.Terrain.slope(srtm).rename('slope');
var slopeM = ee.Terrain.slope(merit.select('elv')).multiply(Math.PI / 180).max(0.0017);
var twi = upa.multiply(1e6).divide(90).divide(slopeM.tan()).log().rename('twi');
var pop = ee.ImageCollection('WorldPop/GP/100m/pop')
  .filter(ee.Filter.eq('country', 'NGA')).filter(ee.Filter.eq('year', 2020))
  .mosaic().select('population');
var area = ee.Image.pixelArea().divide(1e6);
var STATIC = null;           // 2nd/98th percentile limits for twi, elev, slope (filled at start-up)

// ---------------------------------------------------------------- model
function periodLayers(start, end) {
  var s = ee.Date(start), e = ee.Date(end).advance(1, 'day');
  var dwCol = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
    .filterDate(ee.Date(end).advance(-364, 'day'), e).filterBounds(aoi).select('label');
  var dw = dwCol.mode().setDefaultProjection(ee.Image(dwCol.first()).projection());
  var water = dw.eq(0);
  function present(m) {
    return m.reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024}).reproject(proj).gt(0.02);
  }
  function dist(p) {
    return p.selfMask().distance(ee.Kernel.euclidean(CAP, 'meters')).unmask(CAP).min(CAP);
  }
  return {
    land: water.not(),
    dStream: dist(present(upa.gt(2)).or(present(water))),
    dWet: dist(present(dw.eq(3))), dCrop: dist(present(dw.eq(4))),
    dBuilt: dist(present(dw.eq(6))), dBare: dist(present(dw.eq(7))),
    precip: ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate(s, e).sum().rename('precip'),
    lst: ee.ImageCollection('MODIS/061/MOD11A2').filterDate(s, e).select('LST_Day_1km')
      .mean().multiply(0.02).subtract(273.15),
    ndvi: ee.ImageCollection('MODIS/061/MOD13Q1').filterDate(s, e).select('NDVI')
      .mean().multiply(0.0001).rename('ndvi')
  };
}
function nrm(img, lim) { return img.subtract(lim[0]).divide(lim[1] - lim[0]).clamp(0, 1); }

function buildMvi(L, lim, W) {
  var sHealth = FACILITIES ? FACILITIES.distance({searchRadius: 20000, maxError: 50}).unmask(20000).min(20000).divide(20000) : ee.Image(0);
  var sRoad = ROADS ? ROADS.distance({searchRadius: 10000, maxError: 50}).unmask(10000).min(10000).divide(10000) : ee.Image(0);
  var S = ee.Image.cat([
    ee.Image(1).subtract(L.dStream.divide(CAP)).rename('s_stream'),
    ee.Image(1).subtract(L.dWet.divide(CAP)).rename('s_wetland'),
    nrm(L.precip, lim.precip).rename('s_precip'),
    nrm(twi, STATIC.twi).rename('s_twi'),
    // temperature suitability peaks near 25 C (Mordecai et al. 2013); zero near 16 and 34 C
    ee.Image(1).subtract(L.lst.subtract(25).abs().divide(9)).clamp(0, 1).rename('s_lst'),
    ee.Image(1).subtract(L.dCrop.divide(CAP)).rename('s_cropland'),
    ee.Image(1).subtract(nrm(elev, STATIC.elev)).rename('s_elev'),
    ee.Image(1).subtract(nrm(slope, STATIC.slope)).rename('s_slope'),
    nrm(L.ndvi, lim.ndvi).rename('s_ndvi'),
    sHealth.rename('s_health'), sRoad.rename('s_road'),
    L.dBuilt.divide(CAP).rename('s_builtup'), L.dBare.divide(CAP).rename('s_bare')
  ]);
  var mvi = ee.Image(0);
  KEYS.forEach(function (k) { mvi = mvi.add(S.select('s_' + k).multiply(W[k])); });
  return mvi.updateMask(L.land).rename('mvi');
}

// ---------------------------------------------------------------- UI
var panel = ui.Panel({style: {width: '380px', padding: '10px'}});
panel.add(ui.Label('Akwa Ibom Malaria Vulnerability', {fontSize: '20px', fontWeight: 'bold'}));
panel.add(ui.Label('Environmental suitability for malaria transmission and the population exposed, by LGA and period. ' +
  'Method: AHP weights + weighted linear combination (after Abbah et al. 2026). ' +
  'Reflects ecology, not reported disease burden.', {fontSize: '12px', color: '#555'}));

var lgaSelect = ui.Select({items: ['All of Akwa Ibom'], value: 'All of Akwa Ibom', style: {stretch: 'horizontal'}});
panel.add(ui.Label('1. Area', {fontWeight: 'bold'}));
panel.add(lgaSelect);

var today = new Date();
var maxD = new Date(today.getTime() - LAG_DAYS * 86400000);
function iso(d) { return d.toISOString().slice(0, 10); }
var years = [];
for (var y = maxD.getUTCFullYear(); y >= 2016; y--) { years.push(String(y)); }
var yearSel = ui.Select({items: years, value: years[Math.min(1, years.length - 1)]});
var seasonSel = ui.Select({items: ['Full year', 'Wet season (Apr-Oct)', 'Dry season (Nov-Mar)', 'Custom dates'], value: 'Full year'});
var startBox = ui.Textbox({value: '2024-01-01', style: {width: '110px'}});
var endBox = ui.Textbox({value: '2024-12-31', style: {width: '110px'}});
function applyPreset() {
  var yr = Number(yearSel.getValue()), s = seasonSel.getValue(), a, b;
  if (s === 'Custom dates') { return; }
  if (s === 'Full year') { a = new Date(Date.UTC(yr, 0, 1)); b = new Date(Date.UTC(yr, 11, 31)); }
  else if (s.indexOf('Wet') === 0) { a = new Date(Date.UTC(yr, 3, 1)); b = new Date(Date.UTC(yr, 9, 31)); }
  else { a = new Date(Date.UTC(yr - 1, 10, 1)); b = new Date(Date.UTC(yr, 2, 31)); }
  if (b > maxD) { b = maxD; }
  if (a < new Date(MIN_DATE)) { a = new Date(MIN_DATE); }
  startBox.setValue(iso(a)); endBox.setValue(iso(b));
}
yearSel.onChange(applyPreset); seasonSel.onChange(applyPreset);
panel.add(ui.Label('2. Period (30 days to about 3 years, from 2016 to ' + iso(maxD) + ')', {fontWeight: 'bold'}));
panel.add(ui.Panel([yearSel, seasonSel], ui.Panel.Layout.flow('horizontal')));
panel.add(ui.Panel([ui.Label('Start'), startBox, ui.Label('End'), endBox], ui.Panel.Layout.flow('horizontal')));
applyPreset();

// weights
var sliders = {};
var wPanel = ui.Panel({style: {shown: false}});
KEYS.forEach(function (k) {
  sliders[k] = ui.Slider({min: 0, max: 0.3, step: 0.005, value: DEFAULT_W[k], style: {stretch: 'horizontal'}});
  wPanel.add(ui.Label(LABELS[k], {fontSize: '11px', margin: '4px 0 0 4px'}));
  wPanel.add(sliders[k]);
});
wPanel.add(ui.Label('Weights are renormalised to sum to 1. Roads and health facilities are ignored until their assets are set in the script.', {fontSize: '11px', color: '#777'}));
var wBtn = ui.Button({label: '3. AHP weights (show / hide)', onClick: function () { wPanel.style().set('shown', !wPanel.style().get('shown')); }});
panel.add(wBtn); panel.add(wPanel);

var runBtn = ui.Button({label: 'Run analysis', disabled: true, style: {stretch: 'horizontal'}});
var status = ui.Label('Starting up…', {color: '#555'});
var results = ui.Panel();
panel.add(runBtn); panel.add(status); panel.add(results);
panel.add(ui.Label('Classes are equal thirds of the 2nd-98th percentile MVI range across Akwa Ibom for the chosen period, so they describe relative vulnerability. ' +
  'Population: WorldPop 2020. Land cover: Dynamic World, 12 months to the period end.', {fontSize: '11px', color: '#777'}));
ui.root.insert(0, panel);
Map.setOptions('SATELLITE');
Map.centerObject(lgas, 9);
Map.addLayer(lgas.style({color: 'ffffff', fillColor: '00000000', width: 1}), {}, 'LGAs');

// populate LGA list and static limits
lgas.aggregate_array('ADM2_NAME').sort().evaluate(function (names) {
  if (names) { lgaSelect.items().reset(['All of Akwa Ibom'].concat(names)); lgaSelect.setValue('All of Akwa Ibom'); }
});
function p2p98(img, geom, scale) {
  return img.reduceRegion({reducer: ee.Reducer.percentile([2, 98]), geometry: geom, scale: scale,
    maxPixels: 1e10, bestEffort: true, tileScale: 4});
}
p2p98(ee.Image.cat([elev, slope, twi]), aoi, 250).evaluate(function (r, err) {
  if (err || !r) { status.setValue('Start-up failed: ' + err); return; }
  STATIC = {elev: [r.elev_p2, r.elev_p98], slope: [r.slope_p2, r.slope_p98], twi: [r.twi_p2, r.twi_p98]};
  runBtn.setDisabled(false); status.setValue('Ready.');
});

// ---------------------------------------------------------------- run
var runId = 0;
function fail(msg) { status.setValue('Error: ' + msg); runBtn.setDisabled(false); }
function fmt(x) { return (Math.round(x)).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ','); }

runBtn.onClick(function () {
  var start = startBox.getValue(), end = endBox.getValue(), lga = lgaSelect.getValue();
  var re = /^\d{4}-\d{2}-\d{2}$/;
  if (!re.test(start) || !re.test(end)) { status.setValue('Dates must be YYYY-MM-DD.'); return; }
  var ds = new Date(start), de = new Date(end), span = (de - ds) / 86400000 + 1;
  if (isNaN(ds) || isNaN(de) || ds < new Date(MIN_DATE) || de > maxD) { status.setValue('Dates must be between ' + MIN_DATE + ' and ' + iso(maxD) + '.'); return; }
  if (span < 30 || span > 1100) { status.setValue('Period must be 30 days to about 3 years.'); return; }

  var W = {}, tot = 0;
  KEYS.forEach(function (k) { W[k] = sliders[k].getValue(); });
  if (!ROADS) { W.road = 0; }
  if (!FACILITIES) { W.health = 0; }
  KEYS.forEach(function (k) { tot += W[k]; });
  if (tot <= 0) { status.setValue('Weights must not all be zero.'); return; }
  KEYS.forEach(function (k) { W[k] = W[k] / tot; });

  var myRun = ++runId;
  runBtn.setDisabled(true); results.clear();
  status.setValue('Normalising rainfall and NDVI for the period…');
  var L = periodLayers(start, end);
  p2p98(ee.Image.cat([L.precip, L.ndvi]), aoi, 500).evaluate(function (r, err) {
    if (myRun !== runId) { return; }
    if (err || !r) { fail(err); return; }
    var lim = {precip: [r.precip_p2, r.precip_p98], ndvi: [r.ndvi_p2, r.ndvi_p98]};
    var mvi = buildMvi(L, lim, W);
    status.setValue('Computing class breaks across Akwa Ibom…');
    p2p98(mvi, aoi, 500).evaluate(function (q, err2) {
      if (myRun !== runId) { return; }
      if (err2 || !q) { fail(err2); return; }
      var lo = q.mvi_p2, hi = q.mvi_p98, step = (hi - lo) / 3, b1 = lo + step, b2 = hi - step;
      var cls = ee.Image(2).where(mvi.lt(b1), 1).where(mvi.gte(b2), 3).updateMask(mvi.mask()).rename('cls');
      var geom = (lga === 'All of Akwa Ibom') ? aoi : lgas.filter(ee.Filter.eq('ADM2_NAME', lga)).geometry();
      var name = (lga === 'All of Akwa Ibom' ? 'AkwaIbom' : lga).replace(/[^A-Za-z0-9]+/g, '_') + '_' + start + '_' + end;

      Map.layers().reset();
      Map.addLayer(mvi.clip(geom), {min: lo, max: hi, palette: ['1a9850', 'fee08b', 'd73027']}, 'MVI (high = vulnerable)');
      Map.addLayer(cls.clip(geom), {min: 1, max: 3, palette: ['1a9850', 'ffffbf', 'd73027']}, 'MVI class', false);
      Map.addLayer(lgas.style({color: 'ffffff', fillColor: '00000000', width: 1}), {}, 'LGAs');
      Map.centerObject(geom);

      status.setValue('Computing statistics…');
      var valid = mvi.mask();
      var sumImg = mvi.addBands([
        ee.Image(1).updateMask(valid).rename('n'),
        pop.updateMask(valid).rename('pop_total'),
        pop.updateMask(cls.eq(3)).rename('pop_high'), pop.updateMask(cls.eq(2)).rename('pop_med'), pop.updateMask(cls.eq(1)).rename('pop_low'),
        area.updateMask(cls.eq(3)).rename('a_high'), area.updateMask(cls.eq(2)).rename('a_med'), area.updateMask(cls.eq(1)).rename('a_low')
      ]);
      sumImg.reduceRegion({reducer: ee.Reducer.sum(), geometry: geom, scale: SCALE, crs: CRS, maxPixels: 1e10, tileScale: 8})
        .evaluate(function (s, err3) {
          if (myRun !== runId) { return; }
          if (err3 || !s) { fail(err3); return; }
          var T = s.pop_total || 1, A = (s.a_high + s.a_med + s.a_low) || 1;
          function line(t) { results.add(ui.Label(t, {margin: '2px 8px', fontSize: '13px'})); }
          results.add(ui.Label(lga + ', ' + start + ' to ' + end, {fontWeight: 'bold'}));
          line('Mean MVI: ' + (s.mvi / s.n).toFixed(3) + '  (class breaks ' + b1.toFixed(3) + ' / ' + b2.toFixed(3) + ')');
          line('High: ' + fmt(s.a_high) + ' km² (' + (100 * s.a_high / A).toFixed(1) + '%), ' + fmt(s.pop_high) + ' people (' + (100 * s.pop_high / T).toFixed(1) + '%)');
          line('Medium: ' + fmt(s.a_med) + ' km², ' + fmt(s.pop_med) + ' people (' + (100 * s.pop_med / T).toFixed(1) + '%)');
          line('Low: ' + fmt(s.a_low) + ' km², ' + fmt(s.pop_low) + ' people (' + (100 * s.pop_low / T).toFixed(1) + '%)');
          line('Total population: ' + fmt(T));
          status.setValue('Preparing GeoTIFF links…');
          var region = geom.bounds();
          mvi.clip(geom).unmask(NODATA).toFloat().getDownloadURL({name: name + '_MVI', scale: SCALE, crs: CRS, region: region, format: 'GEO_TIFF'},
            function (url, e) {
              if (myRun !== runId) { return; }
              if (e || !url) { results.add(ui.Label('MVI GeoTIFF unavailable for this extent (try a single LGA or a shorter period).', {color: '#a33', fontSize: '12px'})); }
              else { results.add(ui.Label({value: 'Download MVI GeoTIFF (NoData = -9999)', targetUrl: url, style: {margin: '6px 8px'}})); }
            });
          cls.clip(geom).unmask(0).toByte().getDownloadURL({name: name + '_class', scale: SCALE, crs: CRS, region: region, format: 'GEO_TIFF'},
            function (url, e) {
              if (myRun !== runId) { return; }
              if (!e && url) { results.add(ui.Label({value: 'Download class GeoTIFF (1 low, 2 medium, 3 high; NoData = 0)', targetUrl: url, style: {margin: '6px 8px'}})); }
              status.setValue('Done. Links expire after about an hour.');
              runBtn.setDisabled(false);
            });
        });
    });
  });
});
