import matplotlib.pyplot as plt
import matplotlib.table
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from matplotlib.collections import LineCollection
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from scipy.stats import boxcox

# Dicts defining options
# Colours for class probabilities
colour_dict = {
    "V": "lightsalmon",  # Very low
    "L": "yellow",  # Low
    "M": "yellowgreen",  # Medium
    "H": "lightsteelblue",  # High
}

# WFD class labels (var, class) => label
# Note that lower conc => better water quality, so labels are a bit confusing
wfd_class_dict = {
    ("tp", "lower"): "Upper Moderate or better\n(< 29.5 µg/l)",
    ("tp", "upper"): "Lower Moderate or worse\n(≥ 29.5 µg/l)",
    ("cyano", "lower"): "Good or better\n(< 1 mg/l)",
    ("cyano", "upper"): "Moderate or worse\n(≥ 1 mg/l)",
    ("colour", "lower"): "Low\n(< 48 mg Pt/l)",
    ("colour", "upper"): "High\n(≥ 48 mg Pt/l)",
    ("chla", "lower"): "Moderate or better\n(< 20 mg/l)",
    ("chla", "upper"): "Poor or worse\n(≥ 20 mg/l)",
}

# Full var names
names_dict = {
    "tp": "Total phosphorus",
    "cyano": "Cyanobacteria",
    "chla": "Chlorophyll-a",
    "colour": "Colour",
}

# Var units
units_dict = {"tp": "µg/l", "cyano": "mg/l", "colour": "mg Pt/l", "chla": "mg/l"}

# Precision for rounding (for display in tables)
prec_dict = {
    "expected_value": 1,
    "rmse": 1,
    "class_error": 1,
    "mcc": 2,
}


class ChemistryDataError(Exception):
    pass


class MyCell(matplotlib.table.CustomCell):
    """ Extending matplotlib tables.
        
        Adapted from https://stackoverflow.com/a/53573651/505698
    """

    def __init__(self, *args, visible_edges, **kwargs):
        super().__init__(*args, visible_edges=visible_edges, **kwargs)
        seg = np.array(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]
        ).reshape(-1, 1, 2)
        segments = np.concatenate([seg[:-1], seg[1:]], axis=1)
        self.edgelines = LineCollection(segments, edgecolor=kwargs.get("edgecolor"))
        self._text.set_zorder(2)
        self.set_zorder(1)

    def set_transform(self, trans):
        self.edgelines.set_transform(trans)
        super().set_transform(trans)

    def draw(self, renderer):
        c = self.get_edgecolor()
        self.set_edgecolor((1, 1, 1, 0))
        super().draw(renderer)
        self.update_segments(c)
        self.edgelines.draw(renderer)
        self.set_edgecolor(c)

    def update_segments(self, color):
        x, y = self.get_xy()
        w, h = self.get_width(), self.get_height()
        seg = np.array(
            [[x, y], [x + w, y], [x + w, y + h], [x, y + h], [x, y]]
        ).reshape(-1, 1, 2)
        segments = np.concatenate([seg[:-1], seg[1:]], axis=1)
        self.edgelines.set_segments(segments)
        self.edgelines.set_linewidth(self.get_linewidth())
        colors = [
            color if edge in self._visible_edges else (1, 1, 1, 0)
            for edge in self._edges
        ]
        self.edgelines.set_edgecolor(colors)

    def get_path(self):
        codes = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
        return Path(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]],
            codes,
            readonly=True,
        )


# def get_months(season):
#     """ Get the start and end months for the specified season.

#     Args:
#         season: Str. One of ['spring', early_summer', 'late_summer', 'winter']

#     Returns:
#         List. [start_month, end_month]
#     """
#     if season == "winter":
#         return ["November", "January"]
#     elif season == "spring":
#         return ["February", "April"]
#     elif season == "early_summer":
#         return ["May", "July"]
#     elif season == "late_summer":
#         return ["August", "October"]


def get_obs_chem_prev_summer(year):
    """ Returns observed water chemistry (TP, chl-a and colour) for the previous summer 
        (i.e. year - 1).
        
    Args:
        year: Int. Forecast year of interest
        
    Returns:
        Data series of aggregated summer water chemistry for (year - 1).
    """
    df = pd.read_excel(
        "./chemistry_data/seasonal_obs_vanemfjorden.xlsx",
        sheet_name="data",
        index_col=0,
    )

    if df.index.max() < (year - 1):
        raise ChemistryDataError(
            f"Historic water chemistry data is not yet available for {year - 1}. "
            "Please update 'seasonal_obs_vanemfjorden.xlsx'."
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
        TP_prevSummer:     Float. Total P measured from the previous summer (ug/l)
    
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
        year. Uses the Bayesian network for TP, cyanobacteria and colour and a naiive 
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
        float(prev_summer_chem["mean_tp_ugpl"]),
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


def get_bn_hindcast_skill():
    """ Reads hindcast skills (1981 to 2019) from cross-validation of the BN.

    Args:
        None

    Returns:
        Dataframe.
    """
    df = pd.read_csv(
        r"./hindcast_skill_scores/all_bn_skill_scores_1981-2019.csv", index_col=0
    )

    # Convert 'class_error' to percent
    df["class_error"] = df["class_error"] * 100

    return df


# def assign_likelihood_label_climate(tercile_prob):
#     """ Returns an appropriate label based on tercile probabilities.

#     Args:
#         tercile_prob: Float. Between 0 and 1. Proportion of members in terceile

#     Returns:
#         Str. Label for tercile
#     """
#     assert 0 <= tercile_prob <= 1, "'tercile_prob' must be between 0 and 1."

#     if tercile_prob < 0.35:
#         return f"Very low ({tercile_prob*100:.0f}%)"
#     elif 0.35 <= tercile_prob < 0.5:
#         return f"Low ({tercile_prob*100:.0f}%)"
#     elif 0.5 <= tercile_prob < 0.65:
#         return f"Medium ({tercile_prob*100:.0f}%)"
#     else:
#         return f"High ({tercile_prob*100:.0f}%)"


def assign_likelihood_label_water_quality(prob):
    """ Returns an appropriate label based on BN probabilities.
    
    Args:
        prob: Float. Between 0 and 1. BN probability of class
        
    Returns:
        Str. Label for tercile
    """
    assert 0 <= prob <= 1, "'prob' must be between 0 and 1."

    if prob < 0.25:
        return f"Very low\n({prob*100:.0f}%)"
    elif 0.25 <= prob < 0.50:
        return f"Low\n({prob*100:.0f}%)"
    elif 0.50 <= prob < 0.75:
        return f"Medium\n({prob*100:.0f}%)"
    else:
        return f"High\n({prob*100:.0f}%)"


def assign_mcc_class(mcc):
    """ Returns labels and colours for the Matthews Correlation Coefficient, as defined in the
        guidance document.

    Args:
        mcc: Float. Matthews Correlation Coefficient

    Returns:
        Tuple of str (label, colour_name)
    """
    assert mcc <= 1, "'mcc' must be less than or equal to 1."

    if mcc < 0.2:
        return ("None", "lightgrey")
    elif 0.2 <= mcc < 0.4:
        return ("Low", "yellow")
    elif 0.4 <= mcc < 0.6:
        return ("Medium", "yellowgreen")
    else:
        return ("High", "lightsteelblue")


def get_overall_confidence(lik_label, mcc_label):
    """ Assigns overall confidence label for TP, cyano and colour based on class likelihood
        and historic skill. See water quality forecast guidance document for details.

    Args:
        lik_label: Str. Label for most likely class, as returned by 
                   assign_likelihood_label_water_quality(most_likely_class_prob)
        mcc_label: Str. Label for MCC historic skill, as returned by 
                   assign_mcc_class(mcc)[0]

    Returns:
        Str. Overall confidence for forecast summary.
    """
    if lik_label[0] == "H":
        return mcc_label
    elif (lik_label[0] == "M") and (mcc_label == "High"):
        return "Medium"
    else:
        return "Low"


# def make_climate_forecast_table(variable, data_dict, png_path):
#     """ Make a summary table for the climate forecast variables (temperature and
#         precipitation).

#     Args:
#         variable:  Str. One of ['precipitation', 'temperature']
#         data_dict: Dict. One key per tercile ['lower', 'normal', 'upper']. Values should
#                    be two-item lists of [tercile_likelihood, hindcast_rocss]
#         png_path:  Raw str. Path for PNG to be created

#     Returns:
#         None. Table is saved as a PNG to the specified path
#     """
#     matplotlib.table.CustomCell = MyCell

#     # Get correct adjectives for terciles
#     assert variable in [
#         "precipitation",
#         "temperature",
#     ], "'variable' must be 'precipitation' or 'temperature'."

#     adj_dict = {
#         "precipitation": ["Drier", "Wetter"],
#         "temperature": ["Cooler", "Warmer"],
#     }

#     # Data for table grid. Note: historic skill is "None" for all climate vars so hard-coded here
#     data = [
#         [
#             f"{adj_dict[variable][0]} than\nnormal",
#             assign_likelihood_label_climate(data_dict["lower"][0]),
#             f"None (ROCSS = {data_dict['lower'][1]:.2f})",
#             "",
#         ],
#         [
#             "Normal",
#             assign_likelihood_label_climate(data_dict["normal"][0]),
#             f"None (ROCSS = {data_dict['normal'][1]:.2f})",
#             "",
#         ],
#         [
#             f"{adj_dict[variable][1]} than\nnormal",
#             assign_likelihood_label_climate(data_dict["upper"][0]),
#             f"None (ROCSS = {data_dict['upper'][1]:.2f})",
#             "",
#         ],
#     ]

#     # Build table
#     col_labels = ["Tercile", "Likelihood", "Historic skill*", "Forecast summary"]
#     table = plt.table(
#         cellText=data, colLabels=col_labels, loc="center", cellLoc="center",
#     )

#     # Layout and basic formatting
#     table.auto_set_font_size(False)
#     table.auto_set_column_width(col=range(len(col_labels)))
#     table.scale(1, 5.5)
#     for (row, col), cell in table.get_celld().items():
#         if row == 0:
#             cell.set_text_props(fontproperties=FontProperties(weight="bold"))

#     # Set colours for columns
#     for row, tercile in enumerate(["lower", "normal", "upper"], start=1):
#         table[(row, 1)].set_facecolor(
#             colour_dict[assign_likelihood_label_climate(data_dict[tercile][0])[0]]
#         )
#         table[(row, 2)].set_facecolor("lightgrey")  # Skill is always None

#     # Fake "merged cells"
#     h = table.get_celld()[(0, 0)].get_height()
#     w = table.get_celld()[(0, 0)].get_width()
#     header = [
#         table.add_cell(
#             pos, 3, w, h, loc="center", facecolor="lightgrey"
#         )  # Skill is always None
#         for pos in [1, 2, 3]
#     ]
#     header[0].visible_edges = "TLR"
#     header[1].visible_edges = "LR"
#     header[2].visible_edges = "BLR"
#     header[1].get_text().set_text(
#         "Forecast confidence is $\\bf{too\\ low}$\nto make %s predictions \nor the coming season"
#         % variable
#     )
#     table.set_fontsize(20)
#     plt.axis("off")
#     plt.savefig(png_path, dpi=300, bbox_inches="tight")


# def make_forecast_period_table(png_path):
#     """ Create table describing when forecasts are produced. This table is static.

#     Args:
#         png_path: Raw str. Path for PNG to be created

#     Returns:
#         None. Table is saved as a PNG to the specified path
#     """
#     matplotlib.table.CustomCell = MyCell

#     # Data for table grid (historic skill is "None" for all climate vars)
#     data = [
#         [
#             "April",
#             "Summer water quality\nEarly summer climate",
#             "May - October\nMay - July",
#         ],
#         ["July", "Late summer climate", "August - October"],
#         ["October", "Early winter climate", "November - January"],
#         ["January", "Late winter climate", "February - April"],
#     ]

#     # Build table
#     col_labels = [
#         "Issued",
#         "Forecast season\nand variable",
#         "Months in forecast",
#     ]
#     table = plt.table(
#         cellText=data, colLabels=col_labels, loc="center", cellLoc="center",
#     )

#     # Layout and basic formatting
#     table.auto_set_font_size(False)
#     table.auto_set_column_width(col=range(len(col_labels)))
#     table.scale(1, 3)
#     for (row, col), cell in table.get_celld().items():
#         if row == 0:
#             cell.set_text_props(fontproperties=FontProperties(weight="bold"))
#             cell.set_height(0.2)
#         if row == 1:
#             cell.set_height(0.2)

#     table.set_fontsize(16)
#     plt.axis("off")
#     plt.savefig(png_path, dpi=300, bbox_inches="tight")


# def make_climate_confidence_score_guide(png_path):
#     """ Create table explaining confidence scores for climate forecasts. This table is static.

#     Args:
#         png_path: Raw str. Path for PNG to be created

#     Returns:
#         None. Table is saved as a PNG to the specified path
#     """
#     matplotlib.table.CustomCell = MyCell

#     # Data for table grid (historic skill is "None" for all climate vars)
#     data = [
#         ["Very low (less than 35%)", "Some skill", "Very low", "High",],
#         ["Low (35% - 50%)", "Some skill", "Low", "Medium",],
#         ["Medium (50% - 65%)", "Some skill", "Medium", "Low",],
#         ["High (65% or more)", "Some skill", "High", "Very low",],
#         ["Any value", "None", "None", "None",],
#     ]

#     cell_colours = [
#         ["lightsalmon", "yellowgreen", "lightsalmon", "lightsteelblue",],
#         ["yellow", "yellowgreen", "yellow", "yellowgreen",],
#         ["yellowgreen", "yellowgreen", "yellowgreen", "yellow",],
#         ["lightsteelblue", "yellowgreen", "lightsteelblue", "lightsalmon",],
#         ["lightgrey", "lightgrey", "lightgrey", "lightgrey",],
#     ]

#     # Build table
#     col_labels = [
#         "Likelihood (% of members\nwhich predict the tercile)",
#         "Historic skill\n(ROCSS)",
#         "Confidence\nthat tercile will\nhappen",
#         "Confidence\nthat tercile won't\nhappen",
#     ]
#     table = plt.table(
#         cellText=data,
#         colLabels=col_labels,
#         cellColours=cell_colours,
#         loc="center",
#         cellLoc="center",
#     )

#     # Layout and basic formatting
#     table.auto_set_font_size(False)
#     table.auto_set_column_width(col=range(len(col_labels)))
#     table.scale(1, 2.5)
#     for (row, col), cell in table.get_celld().items():
#         if row == 0:
#             cell.set_text_props(fontproperties=FontProperties(weight="bold"))
#             cell.set_height(0.4)

#     table.set_fontsize(20)
#     plt.axis("off")
#     plt.savefig(png_path, dpi=300, bbox_inches="tight")


def make_water_quality_confidence_score_guide(png_path):
    """ Create table explaining confidence scores for water quality forecasts. This table is 
        static.

    Args:
        png_path: Raw str. Path for PNG to be created

    Returns:
        None. Table is saved as a PNG to the specified path
    """
    matplotlib.table.CustomCell = MyCell

    # Data for table grid
    data = [
        ["None", "-", "< 0.2"],
        ["Very low", "< 25", "-"],
        ["Low", "25 - 50", "0.2 - 0.4"],
        ["Medium", "50 - 75", "0.4 - 0.6"],
        ["High", "> 75", "> 0.6"],
    ]

    cell_colours = [
        ["lightgrey", "white", "lightgrey"],
        ["lightsalmon", "lightsalmon", "white"],
        ["yellow", "yellow", "yellow"],
        ["yellowgreen", "yellowgreen", "yellowgreen"],
        ["lightsteelblue", "lightsteelblue", "lightsteelblue"],
    ]

    # Build table
    col_labels = [
        "Label",
        "Likelihood\n(%)",
        "Historic skill\n(MCC³)",
    ]
    table = plt.table(
        cellText=data,
        colLabels=col_labels,
        cellColours=cell_colours,
        loc="center",
        cellLoc="center",
    )

    # Layout and basic formatting
    table.auto_set_font_size(False)
    table.auto_set_column_width(col=range(len(col_labels)))
    table.scale(1, 2.5)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(fontproperties=FontProperties(weight="bold"))
            cell.set_height(0.25)

    table.set_fontsize(20)
    plt.axis("off")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")


def make_tp_cyano_colour_forecast_table(
    variable, lik_low, lik_hi, fcst, rmse, cls_err, mcc, png_path
):
    """ Make a summary table for TP, cyano or colour.

    Args:
        variable: Str. Variable of interest. One of ['tp', 'cyano', 'colour']
        lik_low:  Float. Likelihood of lower class (probability between 0 and 1)
        lik_hi:   Float. Likelihood of upper class (probability between 0 and 1)
        fcst:     Float. Predicted value
        rmse:     Float. Hindcast root mean squared error
        cls_err:  Float. Hindcast classification error (%)
        mcc:      Float. Hindcast Matthews Correlation Coefficient
        png_path: Raw str. Path for PNG to be created

    Returns:
        None. Table is saved as a PNG to the specified path
    """
    assert 0 <= lik_low <= 1, "'lik_low' must be a probability between 0 and 1."
    assert 0 <= lik_hi <= 1, "'lik_hi' must be a probability between 0 and 1."
    assert 0 <= cls_err <= 100, "'cls_err' must be between 0 and 100."
    assert mcc <= 1, "'mcc' cannot be greater than 1."

    matplotlib.table.CustomCell = MyCell

    # Get most likely class
    if lik_low > lik_hi:
        pred_class = "lower"
        max_lik = lik_low
    else:
        pred_class = "upper"
        max_lik = lik_hi

    # Data for table grid (historic skill is "None" for all climate vars)
    data = [
        [
            wfd_class_dict[(variable, "lower")],
            assign_likelihood_label_water_quality(lik_low),
            "",
            "",
            "",
            "",
            "",
        ],
        [
            wfd_class_dict[(variable, "upper")],
            assign_likelihood_label_water_quality(lik_hi),
            "",
            "",
            "",
            "",
            "",
        ],
    ]

    # Build table
    col_labels = [
        "WFD class",
        "Likelihood\nof class",
        f"Forecasted\nvalue\n({units_dict[variable]})",
        f"RMSE\n({units_dict[variable]})¹",
        "Classification\nerror (%)²",
        "MCC³",
        "Forecast summary",
    ]
    table = plt.table(
        cellText=data, colLabels=col_labels, loc="center", cellLoc="center",
    )

    # Layout and basic formatting
    table.auto_set_font_size(False)
    table.auto_set_column_width(col=range(len(col_labels)))
    table.scale(1, 6)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(fontproperties=FontProperties(weight="bold"))

    table[(1, 1)].set_facecolor(
        colour_dict[assign_likelihood_label_water_quality(lik_low)[0]]
    )
    table[(2, 1)].set_facecolor(
        colour_dict[assign_likelihood_label_water_quality(lik_hi)[0]]
    )

    # Fake merged cells
    # Historic skill
    h = table.get_celld()[(0, 0)].get_height()
    w = table.get_celld()[(0, 0)].get_width()
    header = [table.add_cell(-1, pos, w, h, loc="center") for pos in [3, 4, 5]]
    [
        cell.set_text_props(fontproperties=FontProperties(weight="bold"))
        for cell in header
    ]
    [cell.set_height(0.15) for cell in header]
    header[0].visible_edges = "TBL"
    header[1].visible_edges = "TB"
    header[2].visible_edges = "TBR"
    header[1].get_text().set_text("Historic skill*")

    # Stat cols
    for col, stat in enumerate([fcst, rmse, cls_err, mcc], start=2):
        patch = [table.add_cell(pos, col, w, h, loc="center") for pos in [1, 2]]
        patch[0].visible_edges = "TLR"
        patch[1].visible_edges = "BLR"
        patch[0].get_text().set_text(f"\n\n{stat}")

    # Set colour just for MCC
    [cell.set_facecolor(assign_mcc_class(mcc)[1]) for cell in patch]

    # Forecast summary
    class_label = wfd_class_dict[(variable, pred_class)].split("(")[0][:-1]
    summary = (
        f"{names_dict[variable]} is expected\nto be "
        + "$\\bf{"
        + class_label.replace(" ", "\\ ")
        + "}$"
    )
    conf = get_overall_confidence(
        assign_likelihood_label_water_quality(max_lik), assign_mcc_class(mcc)[0]
    )
    conf = "Confidence level: $\\bf{" + conf + "}$"
    patch = [table.add_cell(pos, 6, w, h, loc="center") for pos in [1, 2]]
    patch[0].visible_edges = "TLR"
    patch[1].visible_edges = "BLR"
    patch[0].get_text().set_text(summary)
    patch[1].get_text().set_text(conf)

    table.set_fontsize(20)
    plt.axis("off")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")


def make_chla_forecast_table(fcst, rmse, cls_err, mcc, png_path):
    """ Make a summary table for chl-a.

    Args:
        fcst:     Float. Predicted value
        rmse:     Float. Hindcast root mean squared error
        cls_err:  Float. Hindcast classification error (%)
        mcc:      Float. Hindcast Matthews Correlation Coefficient
        png_path: Raw str. Path for PNG to be created

    Returns:
        None. Table is saved as a PNG to the specified path
    """
    assert 0 <= cls_err <= 100, "'cls_err' must be between 0 and 100."
    assert mcc <= 1, "'mcc' cannot be greater than 1."

    matplotlib.table.CustomCell = MyCell

    # Get most likely class
    if fcst < 20:
        pred_class = "lower"
    else:
        pred_class = "upper"

    # Data for table grid (historic skill is "None" for all climate vars)
    data = [
        [wfd_class_dict[("chla", "lower")], "", "", "", "", "", "",],
        [wfd_class_dict[("chla", "upper")], "", "", "", "", "", "",],
    ]

    # Build table
    col_labels = [
        "WFD class",
        "Likelihood\nof class",
        f"Forecasted\nvalue\n({units_dict['chla']})",
        "RMSE\n(mg/l)¹",
        "Classification\nerror (%)²",
        "MCC³",
        "Forecast summary",
    ]
    table = plt.table(
        cellText=data, colLabels=col_labels, loc="center", cellLoc="center",
    )

    # Layout and basic formatting
    table.auto_set_font_size(False)
    table.auto_set_column_width(col=range(len(col_labels)))
    table.scale(1, 6)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(fontproperties=FontProperties(weight="bold"))

    # Fake merged cells
    # Historic skill
    h = table.get_celld()[(0, 0)].get_height()
    w = table.get_celld()[(0, 0)].get_width()
    header = [table.add_cell(-1, pos, w, h, loc="center") for pos in [3, 4, 5]]
    [
        cell.set_text_props(fontproperties=FontProperties(weight="bold"))
        for cell in header
    ]
    [cell.set_height(0.15) for cell in header]
    header[0].visible_edges = "TBL"
    header[1].visible_edges = "TB"
    header[2].visible_edges = "TBR"
    header[1].get_text().set_text("Historic skill*")

    # Stat cols
    lik = "Not available"
    for col, stat in enumerate([lik, fcst, rmse, cls_err, mcc], start=1):
        patch = [table.add_cell(pos, col, w, h, loc="center") for pos in [1, 2]]
        patch[0].visible_edges = "TLR"
        patch[1].visible_edges = "BLR"
        patch[0].get_text().set_text(f"\n\n{stat}")

    # Set colour just for MCC
    [cell.set_facecolor(assign_mcc_class(mcc)[1]) for cell in patch]

    # Forecast summary
    class_label = wfd_class_dict[("chla", pred_class)].split("(")[0][:-1]
    summary = (
        f"{names_dict['chla']} is expected\nto be "
        + "$\\bf{"
        + class_label.replace(" ", "\\ ")
        + "}$"
    )
    conf = "Confidence level: $\\bf{Medium}$"  # Always 'Medium'; See guidance doc
    patch = [table.add_cell(pos, 6, w, h, loc="center") for pos in [1, 2]]
    patch[0].visible_edges = "TLR"
    patch[1].visible_edges = "BLR"
    patch[0].get_text().set_text(summary)
    patch[1].get_text().set_text(conf)

    table.set_fontsize(20)
    plt.axis("off")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")


def make_forecast(year):
    """ Get a water quality forecast for the specified year from the Bayesian network and make
        all the tables needed to display the summary output (TP, chl-a, colour and cyano).

    Args:
        year: Int. Year of interest for forecast

    Returns:
        None. Forecast tables as saved as PNGs for display in the app.
    """
    # Get forecast and historic skill
    print("Running forecast...")
    fcst_df = forecast_summer_chem(year).round(prec_dict)
    skill_df = get_bn_hindcast_skill().round(prec_dict)

    # Make tables
    print("Building summary tables...")
    make_tp_cyano_colour_forecast_table(
        "tp",
        fcst_df.loc["TP", "prob_below_threshold"],
        fcst_df.loc["TP", "prob_above_threshold"],
        fcst_df.loc["TP", "expected_value"],
        skill_df.loc["TP", "rmse"],
        skill_df.loc["TP", "class_error"],
        skill_df.loc["TP", "mcc"],
        "./images/tp_forecast_summary.png",
    )
    plt.close()

    make_tp_cyano_colour_forecast_table(
        "cyano",
        fcst_df.loc["cyano", "prob_below_threshold"],
        fcst_df.loc["cyano", "prob_above_threshold"],
        fcst_df.loc["cyano", "expected_value"],
        skill_df.loc["cyano", "rmse"],
        skill_df.loc["cyano", "class_error"],
        skill_df.loc["cyano", "mcc"],
        "./images/cyano_forecast_summary.png",
    )
    plt.close()

    make_tp_cyano_colour_forecast_table(
        "colour",
        fcst_df.loc["colour", "prob_below_threshold"],
        fcst_df.loc["colour", "prob_above_threshold"],
        fcst_df.loc["colour", "expected_value"],
        skill_df.loc["colour", "rmse"],
        skill_df.loc["colour", "class_error"],
        skill_df.loc["colour", "mcc"],
        "./images/colour_forecast_summary.png",
    )
    plt.close()

    make_chla_forecast_table(
        fcst_df.loc["chla", "expected_value"],
        skill_df.loc["chla", "rmse"],
        skill_df.loc["chla", "class_error"],
        skill_df.loc["chla", "mcc"],
        "./images/chla_forecast_summary.png",
    )
    plt.close()

    print("Done!")

