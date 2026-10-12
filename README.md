# Akwa Ibom MVI: Earth Engine App version

app/akwa_ibom_mvi_app.js   paste into the Earth Engine Code Editor and publish as an App
site/                      index.html + site-config.js for your Netlify repo

1. Code Editor (code.earthengine.google.com): New script, paste the file, click Run. Check the panel appears.
2. Click Apps (top right of the script panel) -> New App. Choose your Cloud project, source = current script,
   permissions = Anyone can view. Publish. Copy the app URL.
3. Put the URL in site/site-config.js as window.MVI_APP, commit and push site/ files. Netlify redeploys.

Edit the weights in the script (DEFAULT_W) or via the sliders in the app. Add ROADS / FACILITIES assets in the script
configuration to include the access criteria.
