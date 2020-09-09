bayes_net_predict_operational <- function(rfile_fpath,
                                          year,
                                          chla_prevSummer,
                                          colour_prevSummer,
                                          TP_prevSummer) {
    #' Make predictions given the evidence provided based on the pre-fitted Bayesian network (saved as
    #' 'Vansjo_fitted_seasonal_GaussianBN_1981-2019.rds'. Drop met nodes for operational forecast,
    #' and only forecast for TP, colour and cyano (not chla)). Remove standard deviation info.
    #'
    #' Args:
    #'     rfile_fpath:       Str. Filepath to fitted BNLearn network object (.rds file)
    #'     year:              Int. Year for prediction
    #'     chla_prevSummer:   Float. Chl-a measured from the previous summer (mg/l)
    #'     colour_prevSummer: Float. Colour measured from the previous summer (mg Pt/l)
    #'     TP_prevSummer:     Float. Total P measured from the previous summer (mg/l)
    #'
    #' Returns:
    #'     R dataframe

    suppressMessages(library(tidyverse))
    suppressMessages(library(bnlearn))

    # Build dataframe from input data
    year <- c(year)
    chla_prevSummer <- c(chla_prevSummer)
    colour_prevSummer <- c(colour_prevSummer)
    TP_prevSummer <- c(TP_prevSummer)

    driving_data <- data.frame(chla_prevSummer, colour_prevSummer, TP_prevSummer)
    row.names(driving_data) <- year

    # Convert any integer cols to numeric
    driving_data[1:ncol(driving_data)] <- lapply(driving_data[1:ncol(driving_data)], as.numeric)

    # Load fitted Bayesian network and pre-calculated std. devs.
    fitted_BN <- readRDS(rfile_fpath)

    # Nodes to make predictions for. Must match nodes present in the fitted BN.
    # Add check that list is sorted alphabetically, as concatenation of final df assumes this
    nodes_to_predict <- sort(c("colour", "cyano", "TP"))

    # Nodes to use when making predictions
    nodes_to_use <- list(
        "TP" = c("TP_prevSummer", "chla_prevSummer"),
        "colour" = c("colour_prevSummer"),
        "cyano" = c("chla_prevSummer", "TP_prevSummer", "colour_prevSummer")
    )

    set.seed(1)

    # Loop over nodes and derive expected values
    expectedValue_li <- vector(mode = "list", length = 0)

    for (node in nodes_to_predict) {
        pred <- predict(fitted_BN,
            data = driving_data,
            node = node,
            method = "bayes-lw",
            from = unlist(nodes_to_use[node], use.names = FALSE),
            n = 10000
        )

        # If node is cyano, then remove the boxcox transformation before adding expected value to list
        if (node == "cyano") {
            pred <- (pred * 0.1 + 1)**(1 / 0.1) # 0.1 is lambda value chosen in transformation
        }
        expectedValue_li[[node]] <- pred # Update list with value for this node
    }

    # Sort alphabetically
    expectedValue_li <- expectedValue_li[order(names(expectedValue_li))]

    # Thresholds to use in classification
    boundaries_list <- list(
        "TP" = 29.5, # Middle of 'Moderate' class
        "colour" = 48.0, # 66th percentile (i.e. upper tercile). No management implications
        "cyano" = 1.0 # M-P boundary is 2.0, but there were only 2 values in this class. Plenty above 2 tho
    )

    # Sort alphabetically
    boundaries_list <- boundaries_list[order(names(boundaries_list))]

    # Data for evidence, converted to named list
    evidence_li <- as.list(driving_data)

    # Empty list to be populated with probability of being below boundary
    prob_li <- vector(mode = "list", length = 0)

    # Loop over nodes and predict probabilities of being within various classes
    for (node in nodes_to_predict) {
        boundary <- unlist(boundaries_list[node], use.names = FALSE)

        # If cyanomax, apply boxcox transformation with lambda=0.1 to boundary first
        if (node == "cyano") {
            boundary <- (boundary^0.1 - 1) / 0.1
        }

        # Horrible hack to make 'query' available in the global namespace. See
        #     https://stackoverflow.com/q/19260580/505698
        query <<- paste("(", node, " < ", boundary, ")", sep = "")
        prob <- cpquery(fitted_BN,
            event = eval(parse(text = query)),
            evidence = evidence_li,
            method = "lw"
        )

        # Round to 2 d.p. Below this, cpquery returns variable results over diff calls
        # Even with rounding, still get some variability in results
        prob <- round(prob, digits = 2)

        prob_li[[node]] <- prob
    }

    # Sort alphabetically
    prob_li <- prob_li[order(names(prob_li))] # Sort alphabetically

    # Build dataframe
    prob_df <- data.frame(
        node = nodes_to_predict,
        threshold = unlist(boundaries_list, use.names = FALSE),
        prob_below_threshold = unlist(prob_li, use.names = FALSE),
        prob_above_threshold = 1 - unlist(prob_li, use.names = FALSE),
        expected_value = signif(unlist(expectedValue_li, use.names = FALSE), 3) # Round to 3 s.f
    )

    return(prob_df)
}