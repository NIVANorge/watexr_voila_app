import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from scipy.stats import boxcox


def get_obs_chem_prev_summer(year):
    """ Returns observed water chemistry (TP, chl-a, colour and cyanobacteria) for the previous
        summer (i.e. year - 1).
        
    Args:
        year: Int. Forecast year of interest
        
    Returns:
        Data series of summer water chemistry for (year - 1).
    """
    df = pd.read_excel(
        "./chemistry_data/seasonal_obs_vanemfjorden.xlsx",
        sheet_name="data",
        index_col=0,
    )

    df = df.loc[year - 1]

    return df


def bayes_net_predict_operational(
    rfile_fpath, year, chla_prev_summer, colour_prev_summer, tp_prev_summer,
):
    """ Make predictions for TP, colour and cyano given the evidence provided (based on the 
        pre-fitted Bayesian network). This function is just a thin "wrapper" around the R 
        function named 'bayes_net_predict_operational' in 'app_utils.R'.
               
    Args:
        rfile_fpath:       Str. Filepath to fitted BNLearn network object (.rds file)
        year:              Int. Year for prediction
        chla_prevSummer:   Float. Chl-a measured from the previous summer  (mg/l)
        colour_prevSummer: Float. Colour measured from the previous summer (mg Pt/l)
        TP_prevSummer:     Float. Total P measured from the previous summer (mg/l)
    
    Returns:
        Dataframe with columns 'year', 'node', 'threshold', 'prob_below_threshold',
       'prob_above_threshold', 'expected_value'
    """
    # Load R script
    ro.r.source("app_utils.R")

    # Call R function with user-specified evidence
    res = ro.r["bayes_net_predict_operational"](
        rfile_fpath, year, chla_prev_summer, colour_prev_summer, tp_prev_summer,
    )

    # Convert back to Pandas df
    with localconverter(ro.default_converter + pandas2ri.converter):
        df = ro.conversion.rpy2py(res)

    # Add 'year' to results as unique identifier
    df["year"] = int(year)
    df.reset_index(drop=True, inplace=True)

    # Add predicted WFD class
    df["wfd_class"] = (df["threshold"] < df["expected_value"]).astype(int)

    return df


def forecast_summer_chem(year):
    """ Produces a forecast of "summer" water chemistry (May to October) for the specified
        year. Uses the Bayesian netowrk for TP, cyanobacteria and colour and a naiive 
        forecast for chla-a.
        
    Args:
        year: Int. Forecast year of interest. Note that observed data from the previous year
              must be available in ./chemistry_data/seasonal_obs_vanemfjorden.xlsx
              
    Returns:
        Dataframe.
    """
    prev_summer_chem = get_obs_chem_prev_summer(year)

    # Use BN for TP, cyano and colour
    fcst_df = bayes_net_predict_operational(
        "vansjo_fitted_bn_1981-2019.rds",
        year,
        float(prev_summer_chem["mean_chla_mgpl"]),
        float(prev_summer_chem["mean_colour_mgPtpl"]),
        float(prev_summer_chem["mean_tp_mgpl"]),
    )

    fcst_df = fcst_df[
        [
            "year",
            "node",
            "threshold",
            "prob_below_threshold",
            "prob_above_threshold",
            "expected_value",
            "wfd_class",
        ]
    ]

    # Use "naiive" forecast for chl-a
    chla_fcst = prev_summer_chem["mean_chla_mgpl"]
    chla_thresh = 20
    chla_class = chla_fcst > chla_thresh

    fcst_df.loc[len(fcst_df) + 1] = [
        year,
        "chla",
        chla_thresh,
        np.NaN,
        np.NaN,
        chla_fcst,
        chla_class,
    ]

    fcst_df = fcst_df.set_index("node")

    return fcst_df
