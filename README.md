# WateXr: Seasonal forecasting for Lake Vansjø

Creating seasonal forecasts of weather and water quality for Lake Vansjø (Vanemfjorden) near Moss. Part of the [WateXr](https://watexr.eu/) project.

Forecasts are produced for the following variables:

 * Precipitation
 * Temperature
 * Total phosphorus
 * Chlorophyll-a
 * Colour
 * Cyanobacteria
 
Predictions for the weather/climate variables are based on [SEAS5](https://gmd.copernicus.org/articles/12/1087/2019/) data from the [European Centre for Medium-Range Weather Forecasts (ECMWF)](https://www.ecmwf.int/); water quality forecasts are based on a Bayesian network developed as part of the WateXr project.
 
The app itself uses [Voilà](https://github.com/voila-dashboards/voila) and Jupyter notebooks.
 
**Note:** Evaluation of the latest SEAS5 predictions for Lake Vansjø indicates the weather models have little/no significant skill in this region. The precipitation and temperature forecasts should therefore **not** be used for operational decision making.

## Running the test app in Docker

 1. Clone this repository and `cd` into the folder containing the `Dockerfile`
 
 2. Build the image
 
            docker build -t voila_test .
 
 
 3. Run the image as a container
 
            docker run -ti --rm -p 8866:8866 voila_test
 
 
 4. Open a new browser tab and navigate to
 
          http://127.0.0.1:8866/
