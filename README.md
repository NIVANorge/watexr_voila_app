# WateXr: Seasonal forecasting for Lake Vansjø

Creating seasonal forecasts of water quality for Lake Vansjø (Vanemfjorden) near Moss. Part of the [WateXr](https://watexr.eu/) project.

Forecasts are produced for the following variables:

 * Mean total phosphorus concentration (mg/l)
 * Mean chlorophyll-a concentration (mg/l)
 * Mean water colour concentration (mgPt/l
 * Maximum cyanobacteria concentration (mg/l)
 
Forecasts are based on a Bayesian network developed as part of the WateXr project, while the app itself uses [Voilà](https://github.com/voila-dashboards/voila) and Jupyter notebooks.
 
**Note:** The original intention was to include seasonal climate variables based on [SEAS5](https://gmd.copernicus.org/articles/12/1087/2019/) data from the [European Centre for Medium-Range Weather Forecasts (ECMWF)](https://www.ecmwf.int/). However, evaluation of these predictions indicates the weather models have little/no significant skill in this region. The precipitation and temperature forecasts have therefore been omitted for the present.

## Running the test app in Docker

 1. Clone this repository and `cd` into the folder containing the `Dockerfile`
 
 2. Build the image
 
            docker build -t voila_test .
 
 
 3. Run the image as a container
 
            docker run -ti --rm -p 8866:8866 voila_test
 
 
 4. Open a new browser tab and navigate to
 
          http://127.0.0.1:8866/
