# Libraries ---------------------------------------------------------------
packages <- c('here', 'tidyverse', 'ggplot2', 'depmixS4', 'tidyr', 'dplyr', 'zoo', 
              'tseries','urca', 'rugarch', 'VineCopula', 'stats', 'GGally', 
              'reshape2', 'pbapply', 'e1071', 'PerformanceAnalytics', 'knitr', 
              'forecast', 'fGarch', 'quadprog')

instld <- sapply(packages, require, character.only = T)
if (sum(instld)<length(instld)) install.packages(packages[!instld])
invisible(lapply(packages, library, character.only = TRUE))
rm(packages, instld)

# Sourcing Data -----------------------------------------------------------
data_setup <- function(prices){
  # Setup the Data and the working directory
  # Args:
  #   prices: data frame of prices/spreads
  
  if (!is.null(colname_replace_pattern)) {
    # clean up bank names if needed
    colnames(prices) <- gsub(colname_replace_pattern, colname_replacement, colnames(prices))
  }
  
  if (exists("bank_abbr")) {
    # rename columns using tickers
    colnames(prices)[-1] <- sapply(colnames(prices)[-1], bank_abbr)
  }
  
  # setup directories for plots and tables
  if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)
  if (!dir.exists(table_folder)) dir.create(table_folder, recursive = TRUE)
}

calculate_returns <- function(prices, price_columns){
  # Calculate Log Returns from the `prices` data frame, only for `price_columns` columns
  # Args:
  #   prices: data frame of prices/spreads
  #   price_columns: vector of numeric columns
  # Returns:
  #   Log Returns data frame
  returns <- as.data.frame(lapply(prices[, price_columns], function(x) diff(log(x))))
  returns <- cbind(prices$Date[-1], returns)
  
  colnames(returns) <- colnames(prices)
  returns$Date <- as.Date(returns$Date)
  return(returns)
}

# Missing Values ----------------------------------------------------------
# Plot of missing values
missing_values_plot <- function(df_long){
  ggplot(df_long, aes(x = Date, y = Ticker, fill = is.na(Return))) +
    geom_tile() +
    scale_fill_manual(values = c("darkgreen", "red")) + 
    scale_x_date(date_labels = "%Y", date_breaks = "3 year") + 
    labs(title = "Missing Values in Return Columns", x = "Date", y = "Ticker")
}

# Pivot long returns and filter dates
pivot_returns_long <- function(returns_long, first_date) {
  message("Pivoting returns_long to wide format and filtering by first_date...")
  returns_long %>%
    pivot_wider(names_from = "Ticker", values_from = "Return") %>%
    filter(Date >= floor_date(first_date, 'month'))
}

# Compute share of missing values
compute_missing_share <- function(ret_wide, returns) {
  message("Computing share of missing values...")
  colSums(is.na(ret_wide[-1])) / nrow(returns)
}

# Exclude banks with too many missing values
exclude_high_missing <- function(miss_values, threshold = 0.1) {
  message("Identifying banks to exclude with missing share above threshold = ", threshold)
  names(miss_values[miss_values >= threshold])
}

# Drop excluded banks and update objects
update_returns <- function(returns, ret_with_missing, banks_to_exclude, price_columns, bank_names) {
  # Before excluding
  miss_values_before <- compute_missing_share(ret_with_missing, returns)
  miss_values_before <- data.frame(
    Bank = names(miss_values_before),
    MissingShare = as.numeric(miss_values_before),
    Status = ifelse(names(miss_values_before) %in% banks_to_exclude, "Excluded", "Kept")
  )
  
  miss_values_before_table <- kable(
    miss_values_before, format = "latex", booktabs = TRUE, digits = 3,
    caption = "Share of Missing Values Before Exclusion"
  )
  writeLines(miss_values_before_table, paste0(table_folder, data_type, "_miss_values_before_table.tex"))
  message("Saved LaTeX table of missing values *before exclusion* to: ",
          paste0(table_folder, data_type, "_miss_values_before_table.tex"))
  
  
  message("Excluding banks with too many missing values and updating objects...")
  returns <- returns %>% dplyr::select(!any_of(banks_to_exclude))
  ret_with_missing <- ret_with_missing %>% dplyr::select(!any_of(banks_to_exclude))
  price_columns <- price_columns[!price_columns %in% banks_to_exclude]
  bank_names <- bank_names[!bank_names %in% banks_to_exclude]
  
  miss_values <- colSums(is.na(ret_with_missing[-1])) / nrow(returns)
  miss_values_table <- kable(miss_values, format = "latex", booktabs = TRUE,
                             caption = "Share of Missing Values", digits = 3)
  writeLines(miss_values_table, paste0(table_folder, data_type, "_miss_values_table.tex"))
  message("Saved LaTeX table of missing values to: ", paste0(table_folder, data_type, "_miss_values_table.tex"))
  
  list(
    returns = returns,
    ret_with_missing = ret_with_missing,
    price_columns = price_columns,
    bank_names = bank_names
  )
}

# Spline Interpolation
spline_interp <- function(df){
  x <- 1:nrow(df)
  for (col in names(df)[-1]){
    na_indices <- which(is.na(df[[col]]))
    if (length(na_indices) > 0) {
      df[na_indices, col] <- spline(x, df[[col]], xout = na_indices)$y
      message("Interpolated missing values in column: ", col)
    }
  }
  df
}

# Interpolate missing values using cubic spline
interpolate_missing <- function(returns, data_type, first_date, special_bank) {
  message("Interpolating missing values using cubic splines...")
  
  # Handle the special bank column separately
  message("Handling special bank: ", special_bank)
  lastcol <- returns %>% dplyr::select(Date, all_of(special_bank)) %>% filter(Date >= first_date)
  lastcol <- spline_interp(lastcol)
  
  main <- returns %>% dplyr::select(-all_of(special_bank))
  interpolated <- spline_interp(main) %>%
    filter(Date >= first_date) %>%
    left_join(lastcol, by = "Date")
  
  message("Interpolation complete.")
  interpolated
}

# Stationarity ------------------------------------------------------------

# Check the Stationarity of Columns
check_stationarity <- function(data) {
  message("Running stationarity checks (ADF and KPSS)...")
  results <- list()
  cols <- colnames(data)[-1]  # Assuming the first column is 'Date' or similar
  
  for (col in cols) {
    x <- data[[col]]
    
    if (all(is.na(x))) {
      message("Column ", col, " is all NA — skipping.")
      next
    }
    
    # ADF Test
    adf <- tryCatch(
      adf.test(x, alternative = "stationary"),
      error = function(e) {
        message("ADF test failed for ", col, ": ", e$message)
        return(NULL)
      }
    )
    
    # KPSS Test
    kpss <- tryCatch(
      ur.kpss(x, type = "tau"),
      error = function(e) {
        message("KPSS test failed for ", col, ": ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(adf) && !is.null(kpss)) {
      results[[col]] <- data.frame(
        Column = col,
        ADF_Statistic = adf$statistic,
        ADF_p.value = adf$p.value,
        ADF_Result = ifelse(adf$p.value < 0.05, "Stationary", "Non-Stationary"),
        KPSS_Statistic = kpss@teststat,
        KPSS_Result = ifelse(kpss@teststat > kpss@cval[1], "Non-Stationary", "Stationary")
      )
      message("Tested column: ", col, " → ADF: ", results[[col]]$ADF_Result,
              ", KPSS: ", results[[col]]$KPSS_Result)
    }
  }
  
  res <- do.call(rbind, results)
  rownames(res) <- NULL
  message("Stationarity check complete.")
  return(res)
}

# Function to save LaTeX table
save_stationarity_table <- function(results, caption_prefix, table_folder, data_type) {
  message("Saving stationarity results to LaTeX...")
  caption <- paste(caption_prefix, "-", "Stationarity")
  stationarity_table <- kable(results, format = "latex", booktabs = TRUE,
                              caption = caption, digits = 3)
  filepath <- paste0(table_folder, data_type, "_stationarity_table.tex")
  writeLines(stationarity_table, filepath)
  message("LaTeX table saved to: ", filepath)
  return(results)
}


# VIX ---------------------------------------------------------------------

# Align VIX with the returns data
prepare_vix <- function(vix_raw, returns, spline_interp_func = spline_interp) {
  message("Preparing VIX data...")
  
  # Align dates
  first_date <- max(index(vix_raw)[1], returns$Date[1])
  message("First date for alignment: ", as.character(first_date))
  
  returns <- returns %>% filter(Date >= first_date)
  vix <- vix_raw[index(vix_raw) >= first_date]
  
  # Convert to data frame and join
  vix_df <- data.frame(vix = as.numeric(vix), Date = index(vix))
  rownames(vix_df) <- NULL
  
  vix_df <- left_join(returns %>% dplyr::select(Date), vix_df, by = "Date")
  
  # Interpolate missing VIX values
  vix_df <- spline_interp_func(vix_df)
  missing_count <- sum(is.na(vix_df))
  if (missing_count > 0) {
    warning("There are still ", missing_count, " missing VIX values after interpolation.")
  } else {
    message("All VIX values successfully interpolated.")
  }
  
  return(vix_df)
}

# Dataset Split -----------------------------------------------------------
# Split the `returns` and `vix` datasets into train and backtest sets
split_backtest_set <- function(returns, vix, test_days = 500) {
  message("Splitting returns and VIX into training and backtest sets...")
  
  # Ensure returns are sorted by Date
  returns <- arrange(returns, Date)
  
  # Create backtest and training sets
  returns_backtest <- tail(returns, test_days)
  returns_train <- head(returns, nrow(returns) - test_days)
  
  # Match VIX by Date
  vix_backtest <- filter(vix, Date %in% returns_backtest$Date)
  vix_train <- filter(vix, Date %in% returns_train$Date)
  
  # Sanity checks
  if (nrow(returns_train) != nrow(vix_train)) {
    warning("Training sets misaligned: returns (", nrow(returns_train), 
            ") ≠ VIX (", nrow(vix_train), ")")
  } else {
    message("Training sets aligned: ", nrow(returns_train), " observations.")
  }
  
  if (nrow(returns_backtest) != nrow(vix_backtest)) {
    warning("Backtest sets misaligned: returns (", nrow(returns_backtest), 
            ") ≠ VIX (", nrow(vix_backtest), ")")
  } else {
    message("Backtest sets aligned: ", nrow(returns_backtest), " observations.")
  }
  
  return(list(
    returns = returns_train,
    returns_backtest = returns_backtest,
    vix = vix_train,
    vix_backtest = vix_backtest
  ))
}

# Descriptive Statistics and Plots ---------------------------------------------------
compute_descriptive_stats <- function(returns, returns_backtest, price_columns, long_caption, table_folder, data_type) {
  message("Computing descriptive statistics...")
  
  combined <- rbind(dplyr::select(returns, Date, all_of(price_columns)), returns_backtest)
  
  stats <- summarise(combined, across(all_of(price_columns), list(
    Mean = ~ mean(.),
    SD = ~ sd(.),
    Min = ~ min(.),
    Max = ~ max(.),
    Skewness = ~ moments::skewness(.),
    Kurtosis = ~ moments::kurtosis(.),
    N = ~ sum(!is.na(.))
  ), .names = "{.col}_{.fn}"))
  
  stats_long <- stats %>%
    pivot_longer(cols = everything(),
                 names_to = c("Variable", "Statistic"),
                 names_sep = "_") %>%
    pivot_wider(names_from = Statistic, values_from = value)
  
  latex_table <- kable(stats_long, format = "latex", booktabs = TRUE,
                       caption = paste(long_caption, "-", "Descriptive Statistics"), digits = 3)
  
  path <- paste0(table_folder, data_type, "_descriptive_stats_table.tex")
  writeLines(latex_table, path)
  
  message("Descriptive statistics table saved to: ", path)
  return(stats_long)
}

volatility_clustering <- function(returns, returns_backtest, price_columns, banks_to_plot) {
  message("Plotting volatility clustering for ", data_type, " data...")
  
  # Combine returns and convert to xts
  combined_df <- rbind(dplyr::select(returns, Date, all_of(price_columns)), returns_backtest)
  combined_xts <- xts(combined_df[, banks_to_plot], order.by = combined_df$Date)
 
  return(combined_xts)
}

# VIX ---------------------------------------------------------------------

# 2-state HMM for State Detection
run_vix_regime_detection <- function(vix, seed, returns, data_type, save_path = "Saved Results/Regimes.Rds") {
  message("Fitting 2-state HMM to VIX...")
  
  # Fit 2-state Gaussian HMM
  mod <- depmix(response = vix ~ 1, family = gaussian(), nstates = 2,
                          type = 'viterbi', data = data.frame(coredata(vix)))
  set.seed(seed)
  fit <- fit(mod)
  regimes <- posterior(fit, type = 'viterbi')$state
  
  # Attach regime info to VIX
  vix$regime <- as.factor(regimes)
  vix$segment <- cumsum(c(TRUE, diff(as.numeric(vix$regime)) != 0))
  
  # Compute average VIX by regime
  regime_means <- vix %>%
    group_by(regime) %>%
    summarise(mean_vix = mean(vix), .groups = "drop")
  
  # Assign labels automatically: higher mean VIX = Crisis
  crisis_state <- regime_means$regime[which.max(regime_means$mean_vix)]
  tranquility_state <- regime_means$regime[which.min(regime_means$mean_vix)]
  
  vix$regime_label <- factor(vix$regime,
                             levels = c(crisis_state, tranquility_state),
                             labels = c("Crisis", "Tranquility"))
  
  # Assign regimes to returns
  returns$regime <- as.factor(regimes)
  
  return(list(vix = vix, returns = returns, regimes = regimes))
}

# Generate and Save Plot with VIX States
plot_vix_with_regimes <- function(vix, data_type, plot_folder,
                                  show_title = TRUE,
                                  legend_text_size = 12) {
  message("Plotting VIX with regime coloring...")
  
  vix_plot <- ggplot(vix, aes(x = Date, y = vix, color = regime_label, group = segment)) +
    geom_line() +
    scale_color_manual(values = c("Tranquility" = "forestgreen", "Crisis" = "firebrick")) +
    labs(
      title = if (show_title) "VIX Colored by Regime" else NULL,
      y = "VIX", x = "", color = "Regime"
    ) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = legend_text_size + 1),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(2, "lines"),
      axis.text = element_text(size = 18)
    )
  
  plot_path <- paste0(plot_folder, data_type, "_vix_regime_plot.png")
  ggsave(plot_path, plot = vix_plot, width = 8, height = 5, dpi = 200, bg = "transparent")
  message("Saved regime ggplot to: ", plot_path)
  
  return(vix_plot)
}

# GARCH  ------------------------------------------------------------------

# Fitting GARCH Models to all columns
fit_garch_models <- function(returns, price_columns, 
                             armaorder, dist_grch, garchorder = c(1, 1),
                             special_banks = NULL,
                             special_armaorder = NULL, 
                             special_garchorder = c(1, 1), 
                             log = FALSE, log_file = "log.txt") {
  message("Fitting GARCH models with ARMA(", armaorder[1], ",", armaorder[2], ") and ", dist_grch, " distribution.")
  
  garch_spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = garchorder),
    mean.model = list(armaOrder = armaorder, include.mean = TRUE),
    distribution.model = dist_grch
  )
  
  if (!is.null(special_banks)) {
    garch_spec_indiv <- rugarch::ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = special_garchorder),
      mean.model = list(armaOrder = special_armaorder, include.mean = TRUE),
      distribution.model = dist_grch
    )
  }
  
  garch_fits <- list()
  
  for (col in price_columns) {
    if (log) cat(sprintf("Fitting: %s\n", col), file = log_file, append = TRUE)
    
    spec <- if (!is.null(special_banks) && col %in% special_banks) garch_spec_indiv else garch_spec
    
    fit <- tryCatch({
      rugarch::ugarchfit(spec = spec, data = returns[[col]])
    }, warning = function(w) {
      warning_msg <- sprintf("Warning while fitting %s: %s", col, w$message)
      message(warning_msg)
      if (log) cat(warning_msg, file = log_file, append = TRUE, sep = "\n")
      return(NULL)
    }, error = function(e) {
      error_msg <- sprintf("Error while fitting %s: %s", col, e$message)
      message(error_msg)
      if (log) cat(error_msg, file = log_file, append = TRUE, sep = "\n")
      return(NULL)
    })
    
    # Ensure every column has an entry
    garch_fits[[col]] <- fit
  }
  
  # Extract residuals safely
  std_residuals <- as.data.frame(sapply(price_columns, function(col) {
    fit <- garch_fits[[col]]
    if (!is.null(fit)) {
      residuals(fit, standardize = TRUE)
    } else {
      rep(NA, nrow(returns))
    }
  }))
  colnames(std_residuals) <- price_columns
  
  return(list(fits = garch_fits, residuals = std_residuals))
}

# Diagnistics for residuals, autocorr.
check_residual_diagnostics <- function(std_residuals) {
  message("Running Ljung–Box tests on standardized residuals...")
  
  resid_check <- function(resid) {
    pval_mean <- Box.test(resid, lag = 10, type = "Ljung-Box")$p.value
    pval_vol  <- Box.test(resid^2, lag = 10, type = "Ljung-Box")$p.value
    return(c(pval_mean, pval_vol))
  }
  
  check_df <- as.data.frame(t(sapply(std_residuals, resid_check)))
  colnames(check_df) <- c("pval_mean", "pval_vol")
  check_df$Bank <- colnames(std_residuals)
  message("Done.")
  return(dplyr::select(check_df, Bank, everything()))
}

# GARCH Fit ---------------------------------------------------------------

# QQplot for std. residuals vs. theoretical quantiles of sstd dist.
qqplot_garch_sstd <- function(bank, garch_fits, std_residuals){
  n <- nrow(std_residuals)
  xi <- coef(garch_fits[[bank]])["skew"]
  nu <- coef(garch_fits[[bank]])["shape"]
  qqplot(qsstd(ppoints(n), nu = nu, xi = xi), std_residuals[[bank]],
         main = paste0("QQ Plot: ", bank, " Residuals vs Skewed t-distribution"), 
         ylab = paste("St. Residuals of", bank),
         xlab = sprintf("Theoretical quantiles (xi = %.3f, nu = %.3f)", xi, nu))
  abline(0, 1)
}

# QQplot plot and save for a particular bank
plot_bank_qqplot <- function(garch_fits, std_residuals, data_type, plot_folder, bank, qqplot_function){
  b_name <- paste0("_qqplot_garch_", tolower(bank), ".png")
  png(paste0(plot_folder, data_type, b_name), width = 1500, height = 1200,
      res = 200, bg = "transparent")
  qqplot_function(bank, garch_fits = garch_fits, std_residuals = std_residuals)
  dev.off()
  qqplot_function(bank, garch_fits = garch_fits, std_residuals = std_residuals)
}

# KS test
kolmogorov_smirnov_test <- function(bank, garch_fits, std_residuals, B = 2000, seed = 7) {
  # Extract skewness and shape parameters
  xi <- coef(garch_fits[[bank]])["skew"]
  nu <- coef(garch_fits[[bank]])["shape"]
  N <- nrow(std_residuals)
  
  # Empirical test statistic
  test_st <- ks.test(std_residuals[[bank]], "psstd", nu = nu, xi = xi)$statistic
  
  # Bootstrap under null
  set.seed(seed)
  bootstrap_stats <- replicate(B, {
    sample <- rsstd(N, nu = nu, xi = xi)
    ks.test(sample, "psstd", nu = nu, xi = xi)$statistic
  })
  p_value <- mean(bootstrap_stats >= test_st)
  
  # Histogram of bootstrap distribution
  plot_data <- data.frame(ks_stat = bootstrap_stats)
  plot <- ggplot(plot_data, aes(x = ks_stat)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "white") +
    geom_vline(xintercept = test_st, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = "Bootstrap KS Statistic Distribution",
         subtitle = bank,
         x = "KS Statistic", y = "Frequency") +
    theme_minimal()
  
  return(list(
    ks_statistic = as.numeric(test_st),
    p_value = p_value,
    bootstrap_distribution = bootstrap_stats,
    plot = plot
  ))
}

save_ks_plot <- function(ks_test_results, bank, data_type, plot_folder) {
  plot <- ks_test_results[[bank]]$plot
  filename <- paste0(plot_folder, data_type, "_kstest_", gsub(" ", "_", tolower(bank)), ".png")
  
  ggsave(filename, plot = plot,
         width = 10, height = 8, dpi = 300, bg = "transparent")
  
  message("Saved KS plot for ", bank, " to: ", filename)
  plot
}

save_ks_test_table <- function(ks_test_results, data_type, table_folder) {
  # Extract stats and p-values into a data frame
  ks_test_df <- data.frame(
    Bank = names(ks_test_results),
    KS_stat = sapply(ks_test_results, function(x) x$ks_statistic),
    P_value = sapply(ks_test_results, function(x) x$p_value)
  )
  
  # Create LaTeX table
  ks_test_table <- kable(ks_test_df, format = "latex", booktabs = TRUE, 
                         caption = "Kolmogorov–Smirnov Test for Standardized Residuals",
                         digits = 3, row.names = FALSE)
  
  # Save
  file_path <- paste0(table_folder, data_type, "_ks_test_table.tex")
  writeLines(ks_test_table, file_path)
  
  message("Saved KS test table to: ", file_path)
  ks_test_df
}


# Pseudo-observations -----------------------------------------------------

# Empirical CDF on columns
pit <- function(df){
  # calculate ecdf(value = stand. residuals)
  cols <- colnames(df)
  helper <- function(col){
    sapply(df[[col]], emp_cdfs[[col]])
  }
  as.data.frame(sapply(cols, helper))
}

save_pairplots <- function(u, data_type, plot_folder, banks_per_plot = 5) {
  message("Creating pair plots with ", banks_per_plot, " banks per panel...")
  
  n_cols <- ncol(u)
  indices <- split(1:n_cols, ceiling(seq_along(1:n_cols) / banks_per_plot))
  
  # Combine last group if it has only 1 column and more than 1 group exists
  if (length(indices) > 1 && length(indices[[length(indices)]]) == 1) {
    # Combine last two groups
    indices[[length(indices) - 1]] <- c(indices[[length(indices) - 1]], indices[[length(indices)]])
    indices <- indices[-length(indices)]
  }
  
  for (i in seq_along(indices)) {
    idx <- indices[[i]]
    col_subset <- u[, idx]
    
    pair_plot <- GGally::ggpairs(col_subset,
                                 upper = list(continuous = wrap("cor", method = "kendall")),
                                 lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5))) +
      theme_minimal()
    
    col_names <- colnames(u)[idx]
    filename <- paste0(data_type, "_pairplot_", col_names[1], "_", tail(col_names, 1), ".png")
    file_path <- file.path(plot_folder, filename)
    
    ggsave(file_path, plot = pair_plot, width = 10, height = 10, dpi = 300, bg = "transparent")
    message("Saved pairplot: ", file_path)
  }
}


# RVINE Fitting -----------------------------------------------------------

# Custom plotting of the tree structure
custom_rvine_plot <- function(rvine, tree, col, type = 0){
  labelcex <- ifelse(type == 0, 0.8, 0.6)
  plot(rvine, tree = tree, vertex.cex = 1.2, type = type,
       pad = 0.8, # Smaller nodes
       label.cex = labelcex,       
       label.col = 'black', # Smaller label font
       boxed.labels = FALSE,    # No boxes around labels
       edge.col = "gray50",     # Lighter edges for clarity
       vertex.col = col,  # Node color
       edge.lwd = 1,            # Edge thickness
       usearrows = FALSE,       # For undirected trees
       mode = "fruchtermanreingold")
}

# Saving the tree plot to a specified folder with a specified name
save_rvine_tree_plot <- function(vine_object, tree = 1, type, col,
                                 plot_folder, data_type,
                                 filename_prefix) {
  # Construct file name
  filename <- paste0(plot_folder, data_type, "_", filename_prefix, "_tree", tree, ".png")
  
  # Save the plot
  png(filename, width = 1600, height = 1200, res = 200, bg = "transparent")
  custom_rvine_plot(vine_object, tree = tree, type = type, col = col)
  dev.off()
  
  message("Saved tree plot: ", filename)
}

# Structurised Summary 
vine_trees_summary <- function(summary_df) {
  df <- as.data.frame(summary_df)
  df <- df[, c("tree", "edge", "cop", "family", "par", "par2", "tau", "utd", "ltd")]
  df$edge <- gsub(";", " | ", df$edge)
  df$copula <- BiCopName(df$family, short = FALSE)
  
  df <- df[, c("tree", "edge", "copula", "par", "par2", "tau", "utd", "ltd")]
  colnames(df) <- c("Tree", "Edge", "Copula", "Param. 1", "Param. 2",
                    "Tau", "Upp. Tail", "Low. Tail")
  rownames(df) <- NULL
  
  # Split by tree
  tree_list <- split(df, df$Tree)
  return(tree_list)
}

# Save summaries for different trees
save_vine_tree_tables <- function(vine_model, data_type, model_name, table_folder, 
                                  trees_to_export = 1, digits = 3, tree_to_print = NULL) {
  # Compute the summary and clean it
  s <- summary(vine_model)
  tree_summaries <- vine_trees_summary(s)
  
  # For each tree requested
  for (tree in trees_to_export) {
    latex_table <- kable(tree_summaries[[tree]],
                         format = "latex", booktabs = TRUE,
                         caption = paste(model_name, "Vine Copula, Tree", tree),
                         digits = digits, row.names = FALSE)
    
    file_path <- paste0(table_folder, data_type, "_", model_name, "_summary_tree", tree, ".tex")
    writeLines(latex_table, file_path)
    message("\nSaved: ", file_path)
  }
  
  if (!is.null(tree_to_print)){
    tree_summaries[[tree_to_print]]
  }
}

# Compare Tails and Tau across different trees
compare_tree_dependence <- function(vine1, vine2, tree = 1) {
  t1 <- vine_trees_summary(summary(vine1))[[tree]]
  t2 <- vine_trees_summary(summary(vine2))[[tree]]
  
  stats <- c("Tau", "Upp. Tail", "Low. Tail")
  means1 <- sapply(stats, function(stat) mean(t1[[stat]], na.rm = TRUE))
  means2 <- sapply(stats, function(stat) mean(t2[[stat]], na.rm = TRUE))
  
  result <- data.frame(
    Statistic = stats,
    Model1 = round(means1, 3),
    Model2 = round(means2, 3),
    Difference = round(means1 - means2, 3)
  )
  return(result)
}

# Information Criteria -----------------------------------------------------

# IC
compute_ic_values <- function(u, u_regime1, u_regime2,
                              model_single, model_regime1, model_regime2,
                              model_cvine, regimes, returns) {
  
  # Log-likelihoods
  ll_single <- RVineLogLik(u, model_single, separate = TRUE)$loglik
  ll_r1 <- RVineLogLik(u_regime1, model_regime1, separate = TRUE)$loglik
  ll_r2 <- RVineLogLik(u_regime2, model_regime2, separate = TRUE)$loglik
  
  # Degrees of freedom
  df1 <- (RVineAIC(u_regime1, model_regime1)$AIC + 2 * model_regime1$logLik) / 2
  df2 <- (RVineAIC(u_regime2, model_regime2)$AIC + 2 * model_regime2$logLik) / 2
  df_single <- (RVineAIC(u, model_single)$AIC + 2 * model_single$logLik) / 2
  
  aic_single <- RVineAIC(u, model_single)$AIC
  bic_single <- RVineBIC(u, model_single)$BIC
  
  aic_regime <- -2 * (model_regime1$logLik + model_regime2$logLik) + 2 * (df1 + df2)
  bic_regime <- -2 * (model_regime1$logLik + model_regime2$logLik) +
    (log(nrow(u_regime1)) * df1 + log(nrow(u_regime2)) * df2)
  
  aic_cvine <- RVineAIC(u, model_cvine)$AIC
  bic_cvine <- RVineBIC(u, model_cvine)$BIC
  df_cvine <- (aic_cvine + 2 * model_cvine$logLik) / 2
  
  # Loglik for display
  loglik_vec <- c(model_single$logLik, model_regime1$logLik + model_regime2$logLik, model_cvine$logLik)
  
  ic_df <- data.frame(
    Model = c("Single R-Vine", "Regime R-Vine", "C-Vine"),
    LogLik = round(loglik_vec, 3),
    AIC = round(c(aic_single, aic_regime, aic_cvine), 3),
    BIC = round(c(bic_single, bic_regime, bic_cvine), 3)
  )
  
  return(list(ic_df = ic_df, df = list(df1 = df1, df2 = df2, df_single = df_single, 
                                       df_cvine = df_cvine), 
              ll = list(ll_single = ll_single, ll_r1 = ll_r1, ll_r2 = ll_r2)))
}

save_ic_table <- function(ic_df, path, caption = "Information Criteria") {
  ic_table <- kable(ic_df, format = "latex", booktabs = TRUE,
                    caption = caption, digits = 3)
  writeLines(ic_table, path)
}


# Vuong Test --------------------------------------------------------------

vuong_test <- function(loglik_r, loglik_s, df_r, df_s, alpha = 0.05, verbose = TRUE) {
  ll_diff <- loglik_r - loglik_s
  mean_diff <- mean(ll_diff)
  sd_diff <- sd(ll_diff)
  n <- length(ll_diff)
  
  vuong_stat <- sqrt(n) * mean_diff / sd_diff
  p_value <- 2 * (1 - pnorm(abs(vuong_stat)))
  
  penalty <- (df_r - df_s) * log(n) / (2 * n)
  vuong_stat_adj <- sqrt(n) * (mean_diff - penalty) / sd_diff
  p_value_adj <- 2 * (1 - pnorm(abs(vuong_stat_adj)))
  
  result <- ifelse(p_value < alpha,
                   ifelse(vuong_stat > 0, "Regimes model preferred (NON-adjusted)", 
                          "Single model preferred (NON-adjusted)"),
                   "No significant difference (NON-adjusted)")
  
  result_adj <- ifelse(p_value_adj < alpha,
                       ifelse(vuong_stat_adj > 0, "Regimes model preferred (BIC-adjusted)", 
                              "Single model preferred (BIC-adjusted)"),
                       "No significant difference (BIC-adjusted)")
  
  res <- data.frame(
    `p-value` = round(c(p_value, p_value_adj), 4),
    `Vuong Statistic` = round(c(vuong_stat, vuong_stat_adj), 3),
    `Model Preference` = c(result, result_adj)
  )
  rownames(res) <- c("Non-adjusted", "Adjusted")
  
  if (verbose) {
    cat("Mean obs-level difference in LogLik (Regimes - Single):", round(mean_diff, 3), "\n")
    cat("Standard deviation of differences:", round(sd_diff, 3), "\n")
  }
  
  return(res)
}

save_vuong_table <- function(vuong_result, table_path, caption = "Vuong Test for Non-nested Models") {
  vuong_table <- kable(vuong_result, format = "latex", booktabs = TRUE,
                       caption = caption, digits = 3)
  writeLines(vuong_table, table_path)
}

# Rolling Vines -----------------------------------------------------------

# Estimate rolling r vines
rolling_rvines <- function(returns, rolling_window, step_size, u=NULL, 
                           compute_u = FALSE, armaorder = NULL, 
                           garchorder = c(1, 1), special_banks = NULL, 
                           special_armaorder = NULL, 
                           special_garchorder = c(1, 1), dist_grch = "sstd", 
                           price_columns = NULL){
  # NO PRE_COMPUTED U
  if (compute_u){
    # General GARCH specification
    garch_spec <- rugarch::ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = garchorder),
      mean.model = list(armaOrder = armaorder, include.mean = TRUE),
      distribution.model = dist_grch
    )
    
    # Optional: special spec for selected banks
    if (!is.null(special_banks)) {
      garch_spec_indiv <- rugarch::ugarchspec(
        variance.model = list(model = "sGARCH", garchOrder = special_garchorder),
        mean.model = list(armaOrder = special_armaorder, include.mean = TRUE),
        distribution.model = dist_grch
      )
    }
    
    # Fit models
    garch_fits <- lapply(price_columns, function(col) {
      spec <- if (!is.null(special_banks) && col %in% special_banks) garch_spec_indiv else garch_spec
      rugarch::ugarchfit(spec, data = returns[[col]])
    })
    names(garch_fits) <- price_columns
    
    std_residuals <- as.data.frame(sapply(garch_fits, residuals, standardize = TRUE))
    colnames(std_residuals) <- price_columns
    
    emp_cdfs <- lapply(std_residuals, ecdf)
    emp_cdfs <- setNames(emp_cdfs, colnames(std_residuals))
    
    pit <- function(df){
      cols <- colnames(df)
      helper <- function(col){
        sapply(df[[col]], emp_cdfs[[col]])
      }
      as.data.frame(sapply(cols, helper))
    }
    
    u <- pit(std_residuals)
  }
  
  # IF U PRE-COMPUTED
  # Start indices
  n_obs <- nrow(u)
  start_indices <- seq(1, n_obs - rolling_window + 1, by = step_size)
  # Estimating RVines
  helper <- function(start_idx){
    u_roll <- u[start_idx:(start_idx + rolling_window - 1), ]
    rvine <- RVineStructureSelect(u_roll, familyset = 0:6)
    list(start_date = returns$Date[start_idx], 
         end_date = returns$Date[start_idx + rolling_window - 1], 
         rvine = rvine)
  }
  pblapply(start_indices, helper)
}

plot_roll_rvine <- function(roll_rvines, indices, plot_folder, data_type, tree = 1,
                            type = "contour", col = "black") {
  for (i in indices) {
    rvine_obj <- roll_rvines[[i]]
    start_date <- as.character(rvine_obj[[1]])
    end_date <- as.character(rvine_obj[[2]])
    cat("From", start_date, "to", end_date, "\n")
    
    file_name <- paste0(plot_folder, data_type, "_roll_rvine_", i, ".png")
    png(file_name, width = 2000, height = 1600, res = 300, bg = "transparent")
    custom_rvine_plot(rvine_obj[[3]], tree = tree, type = type, col = col)
    dev.off()
  }
}


# Simulations -------------------------------------------------------------

# If some estimations faileed - carry forward the previous results
carry_forward_results <- function(results_list) {
  null_indices <- which(sapply(results_list, is.null))
  null_count <- length(null_indices)
  
  if (null_count > 0) {
    message(sprintf("There are %d NULL entries in the passed list. The first NULL is at index %d.", 
                    null_count, null_indices[1]))
    message("Carrying forward last non-null indices...")
    last_valid <- NULL
    cleaned_list <- vector("list", length(results_list))
    replaced_indices <- integer(0)
    
    for (i in seq_along(results_list)) {
      if (is.null(results_list[[i]])) {
        cleaned_list[[i]] <- last_valid
        replaced_indices <- c(replaced_indices, i)
      } else {
        cleaned_list[[i]] <- results_list[[i]]
        last_valid <- results_list[[i]]
      }
    }
    message("Done.")
    return(list(
      cleaned = cleaned_list,
      replaced_indices = replaced_indices
    ))
  } else {
    message("There are no NULL entries in the passed list.")
    return(list(cleaned = results_list, 
                replaced_indices = NULL))
  }
}

# Simulate log returns from an rr-vine
simulated_logret <- function(N, rvine, returns, price_columns, dist_grch,
                             armaorder, garchorder = c(1, 1), special_armaorder, 
                             special_garchorder = c(1, 1), special_banks = NULL, 
                             start_idx = NULL, rolling = TRUE,
                             rolling_window = NULL, log_file = NULL, seed = 7){
  
  if (rolling == TRUE){
    # Returns used to estimate GARCH
    returns_4garch <- returns[start_idx:(start_idx + rolling_window - 1), price_columns]
  } else {
    returns_4garch <- returns[,  price_columns]
  }
  
  # Define general spec
  garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = garchorder),
    mean.model = list(armaOrder = armaorder, include.mean = TRUE),
    distribution.model = dist_grch
  )
  
  # Define special spec
  if (!is.null(special_banks) && !is.null(special_armaorder)) {
    garch_spec_special <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = special_garchorder),
      mean.model = list(armaOrder = special_armaorder, include.mean = TRUE),
      distribution.model = dist_grch
    )
  }
  
  safe_ugarchfit <- function(spec, data, colname) {
    solvers <- c("solnp", "hybrid")
    
    for (solver in solvers) {
      fit <- tryCatch({
        ugarchfit(spec, data = data, solver = solver)
      }, warning = function(w) {
        return(NULL)
      }, error = function(e) {
        return(NULL)
      })
      
      # Check convergence
      if (!is.null(fit) && fit@fit$convergence == 0 && is.finite(fit@fit$LLH)) {
        return(fit)
      }
    }
    
    if (!is.null(log_file)){
      cat(sprintf("[%s] All solvers failed. Returning NULL.\n", colname), file = log_file, append = TRUE)
    }
    return(NULL)
  }
  
  # Fit GARCH models in parallel safely
  garch_fits <- mapply(function(x, colname) {
    spec <- if (!is.null(special_banks) && colname %in% special_banks) garch_spec_special else garch_spec
    safe_ugarchfit(spec, x, colname)
  }, x = returns_4garch, colname = price_columns, SIMPLIFY = FALSE)
  
  # Residuals
  std_residuals <- lapply(garch_fits, function(fit) {
    if (is.null(fit)) {
      rep(0, nrow(returns_4garch))  # fallback vector of zeros
    } else {
      residuals(fit, standardize = TRUE)
    }
  })
  std_residuals <- as.data.frame(std_residuals)
  colnames(std_residuals) <- price_columns
  
  # Simulate pseudo-observations
  set.seed(seed = seed)
  sim_u <- RVineSim(N = N, RVM = rvine)
  
  sim_u <- as.data.frame(matrix(sim_u, nrow = N))  # force 2D format
  colnames(sim_u) <- colnames(std_residuals)
  
  # Create inverse CDFs using interpolation
  inv_cdfs <- lapply(std_residuals, function(x) {
    ecdf_x <- ecdf(x)
    x_vals <- sort(x)
    u_vals <- ecdf_x(x_vals)
    approxfun(u_vals, x_vals, rule = 2, ties = "mean")  # extrapolate at edges
  })
  inv_cdfs <- setNames(inv_cdfs, colnames(std_residuals))
  
  # Get the standardised residuals from the pseudo-observation: the inverse cdf to preudo-obs
  sim_z <- as.data.frame(mapply(function(colname, u_col) {
    inv_cdfs[[colname]](u_col)
  }, colname = colnames(sim_u), u_col = as.data.frame(sim_u), SIMPLIFY = FALSE))
  
  if (rolling == FALSE){
    garch_estimates <- lapply(garch_fits, function(fit) {
      if (is.null(fit)) {
        list(sigma = 1, mean = 0)  # fallback values
      } else {
        list(sigma = tail(sigma(fit), 1), mean = tail(fitted(fit), 1))
      }
    })
    sigmas <- sapply(garch_estimates, `[[`, "sigma")
    means <- sapply(garch_estimates, `[[`, "mean")
  } else {
    # Forecast conditional sigma and mean from GARCH (rolling)
    forecasts <- lapply(garch_fits, function(fit) {
      if (is.null(fit)) {
        list(sigma = 1, mean = 0)  # fallback values
      } else {
        fc <- ugarchforecast(fit, n.ahead = 1)
        list(sigma = sigma(fc), mean = fitted(fc))
      }
    })
    sigmas <- sapply(forecasts, `[[`, "sigma")
    means <- sapply(forecasts, `[[`, "mean")
  }
  # Rescale standardized residuals and add the mean
  # THUS: we get one-day-ahead log-ret observations
  sim_returns <- sweep(sim_z, 2, sigmas, `*`)
  sim_returns <- sweep(sim_returns, 2, means, `+`)
  
  sim_returns
}

# Simulate log returns on a rolling basis
rolling_simulated_logret <- function(rvine, returns, returns_test, price_columns, test_window, rolling_window){
  # Cutting the observations earlier than the first obs. for the rolling window
  returns <- returns %>% dplyr::select(-"regime")
  returns <- tail(returns, rolling_window)
  
  # Adding the test observations
  returns <- rbind(returns, returns_test)
  # Test_window-number starting indices, for test_window ahead estimation
  start_indices <- seq(1, test_window)
  simul_returns_list <- pblapply(start_indices, function(x) simulated_logret(start_idx = x, rolling = TRUE, N = 1, 
                                                                             rvine = rvine, price_columns = price_columns, returns = returns, 
                                                                             rolling_window = rolling_window))
  do.call(rbind, simul_returns_list)
}

# HMM ---------------------------------------------------------------------

# Predict one-day-ahead state
HMM_predict <- function(vix, start_idx = NULL, rolling_window = NULL, rolling = TRUE, seed = 7){
  if (rolling == TRUE){
    # Returns used to estimate GARCH
    vix <- vix[start_idx:(start_idx + rolling_window - 1), ]
  } 
  
  mod <- depmix(response = vix ~ 1, family = gaussian(), nstates = 2, 
                type = 'viterbi', data = data.frame(coredata(vix)))
  set.seed(seed = seed)
  invisible(capture.output({ fit <- fit(mod, verbose = FALSE)}))
  regimes <- posterior(fit, type = 'viterbi')$state
  
  # extract the state-transition matrix
  transition_mat <- rbind(getpars(getmodel(fit,"transition",1)),
                          getpars(getmodel(fit,"transition",2)))
  
  # extract the probability of the states at the final time point in the data (t=T)
  # this will act as a "prior" to compute the forecasted state distributions
  prior_vec <- as.numeric(posterior(fit, type = 'viterbi')[length(regimes),-1])
  
  # for T + 1, the forecasted state distribution is (prior_vec %*% transition_mat)
  state_oneahead_prob <- prior_vec %*% transition_mat
  which.max(state_oneahead_prob)
}

HMM_rolling <- function(vix, vix_test, rolling_window, 
                        test_window, seed = 7){
  # Cutting the observations earlier than the first obs. for the rolling window
  vix <- vix %>% dplyr::select(c("Date", "vix"))
  vix <- tail(vix, rolling_window)
  
  # Adding the test observations
  vix <- rbind(vix, vix_test)
  # Test_window-number starting indices, for test_window ahead estimation
  start_indices <- seq(1, test_window)
  simul_states <- pbsapply(start_indices, function(x) HMM_predict(start_idx = x, vix = vix, seed = seed,
                                                                  rolling_window = rolling_window))
  simul_states
}

future_vix_plot <- function(vix_backtest_to_plot, plot_folder, data_type, vix_future_states,
                            show_title = TRUE,
                            legend_text_size = 12,
                            legend_text_face = "plain") {
  plt <- ggplot(vix_backtest_to_plot %>% 
                  mutate(segment = cumsum(c(TRUE, diff(as.numeric(regime)) != 0))),
                aes(x = Date, y = vix, color = regime_label, group = segment)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Tranquility" = "forestgreen", "Crisis" = "firebrick")) +
    labs(
      title = if (show_title) "VIX Colored by Regime" else NULL,
      y = "VIX", x = "", color = "Regime"
    ) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = legend_text_size + 1, face = legend_text_face),
      legend.text = element_text(size = legend_text_size, face = legend_text_face), 
      legend.key.size = unit(2, "lines")
    )
  
  plot_path <- paste0(plot_folder, data_type, "_vix_future_states.png")
  ggsave(plot_path, plot = plt,
         width = 8, height = 5, dpi = 200, bg = "transparent")
  message("Saved VIX future states plot to: ", plot_path)
  plt
}

# Risk Metrics: VaR and ES --------------------------------------------------------------

VaR <- function(sim_returns, weights, alpha){
  port_returns <- as.vector(as.matrix(sim_returns) %*% weights)
  quantile(port_returns, probs = alpha)
}

ES <- function(sim_returns, weights, alpha){
  port_returns <- as.vector(as.matrix(sim_returns) %*% weights)
  VaR <- VaR(sim_returns, weights, alpha)
  mean(port_returns[port_returns <= VaR])
}


# Portfolio Objectives ----------------------------------------------------

# Min var objective
min_var_weights <- function(sim_returns, alpha) {
  d <- ncol(sim_returns)
  init_weights <- rep(1/d, d)
  opt <- optim(
    init_weights,
    fn = function(w) -VaR(weights = w / sum(w), sim_returns = sim_returns, alpha = alpha),
    method = "L-BFGS-B",
    lower = rep(0, d),
    upper = rep(1, d)
  )
  opt_weights <- opt$par / sum(opt$par)
  return(list(weights_mvar = opt_weights, 
              VaR_mvar = VaR(weights = opt_weights, sim_returns = sim_returns, alpha = alpha), 
              ES_mvar = ES(weights = opt_weights, sim_returns = sim_returns, alpha = alpha)))
}

# MaxIR objective
max_ir_weights <- function(sim_returns, alpha) {
  d <- ncol(sim_returns)
  init_weights <- rep(1/d, d)
  opt <- optim(
    init_weights,
    fn = function(w) {
      w <- w / sum(w)
      port_ret <- as.vector(as.matrix(sim_returns) %*% w)
      -mean(port_ret) / abs(VaR(weights = w, sim_returns = sim_returns, alpha = alpha))
    },
    method = "L-BFGS-B",
    lower = rep(0, d),
    upper = rep(1, d)
  )
  opt_weights <- opt$par / sum(opt$par)
  
  return(list(weights_ir = opt_weights, 
              VaR_ir = VaR(weights = opt_weights, sim_returns = sim_returns, alpha = alpha), 
              ES_ir = ES(weights = opt_weights, sim_returns = sim_returns, alpha = alpha)))
}

# Backtesting -------------------------------------------------------------

# Optimise portfolio weights based on the simulated one-day-ahead returns
rolling_risk <- function(rvine, returns, returns_test, price_columns, test_window, 
                         rolling_window, N, solver_spec = NULL, seed = 7,
                         dist_grch, armaorder, special_armaorder = NULL, special_banks = NULL,
                         garchorder = c(1, 1), special_garchorder = c(1, 1), 
                         var_alpha = 0.01){
  # Cutting the observations earlier than the first obs. for the rolling window
  returns <- returns %>% dplyr::select(-"regime")
  returns <- tail(returns, rolling_window)
  
  # Adding the test observations
  returns <- rbind(returns, returns_test)
  # Test_window-number starting indices, for test_window ahead estimation
  start_indices <- seq(1, test_window)
  
  message("Initializing cluster with ", detectCores() - 1, " cores...")
  cl <- makeCluster(detectCores() - 1)
  on.exit({
    if (!is.null(cl)) stopCluster(cl)
  }, add = TRUE)
  
  message("Exporting variables to cluster...")
  clusterExport(cl, varlist = c("simulated_logret",
                                "VaR", "ES", "min_var_weights", "max_ir_weights"))
  clusterExport(cl, varlist = c("rvine", "returns", "returns_test", "price_columns", 
                                "N", "rolling_window", "solver_spec", "seed",
                                "dist_grch", "armaorder", "special_armaorder", "special_banks",
                                "garchorder", "special_garchorder", "var_alpha"), envir = environment())
  clusterEvalQ(cl, {
    library(rugarch)
    library(VineCopula)
    library(dplyr)
  })
  
  helper <- function(start_idx){
    log_file <- "log.txt"
    cat(sprintf("Running start_idx = %s\n", start_idx), file = log_file, append = TRUE)
    tryCatch({
      sim_returns <- simulated_logret(start_idx = start_idx, rolling = TRUE, N = N, 
                                      rvine = rvine, price_columns = price_columns, returns = returns, 
                                      rolling_window = rolling_window, log_file = log_file, 
                                      dist_grch = dist_grch, armaorder = armaorder, 
                                      special_armaorder = special_armaorder, special_banks = special_banks,
                                      garchorder = garchorder, special_garchorder = special_garchorder, seed = seed)
      
      weights_eq <- rep(1 / ncol(sim_returns), ncol(sim_returns))
      VaR_eq <- VaR(sim_returns, weights_eq, var_alpha)
      ES_eq <- ES(sim_returns, weights_eq, var_alpha)
      return(list(list(weights_eq = weights_eq, VaR_eq = VaR_eq, ES_eq = ES_eq), 
                  min_var_weights(sim_returns, var_alpha), max_ir_weights(sim_returns, var_alpha)))
    }, error = function(e) {
      error_msg <- sprintf("Error at start_idx = %s: %s", start_idx, e$message)
      message(error_msg)
      cat(error_msg, file = log_file, append = TRUE, sep = "\n")
      return(NULL)
    })
  }
  message("Starting parallel simulation with pblapply...")
  results <- pblapply(start_indices, helper, cl = cl)
  message("Rolling risk analysis complete.")
  results
}

ifdebug <- function(){
  roll_risk <- rolling_risk(rvine = rvine_single, returns = returns, returns_test = returns_backtest, 
                            price_columns, test_window = 500, rolling_window = 1000, N = 1000)
  roll_risk[[425]] #error at 425
  
  roll_debugging <- function(err_index){
    returns_debug <- returns %>% dplyr::select(-"regime")
    returns_debug <- tail(returns_debug, 1000)
    returns_debug <- rbind(returns_debug, returns_backtest)
    sim_returns_debug <- simulated_logret(start_idx = err_index, rolling = TRUE, N = 5000, 
                                          rvine = rvine_single, price_columns = price_columns, returns = returns_debug, 
                                          rolling_window = 1000)
    returns_4garch_debug <- returns_debug[err_index:(err_index + 1000 - 1), price_columns]
    
    garch_spec_debug <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
      distribution.model = "sstd"
    )
    
    safe_ugarchfit_debug <- function(spec, data, bank) {
      solvers <- c("solnp", "hybrid")
      
      log_file <- "log_garch.txt"
      write(sprintf("Fitting Bank %s", bank), file = log_file, append = TRUE)
      for (solver in solvers) {
        msg <- sprintf("Trying solver = %s", solver)
        write(msg, file = log_file, append = TRUE)
        
        fit <- tryCatch({
          ugarchfit(spec, data = data, solver = solver)
        }, warning = function(w) {
          write(sprintf("Warning using %s: %s", solver, w$message), file = log_file, append = TRUE)
          return(NULL)
        }, error = function(e) {
          write(sprintf("Error using %s: %s", solver, e$message), file = log_file, append = TRUE)
          return(NULL)
        })
        
        # Check convergence
        if (!is.null(fit) && fit@fit$convergence == 0 && is.finite(fit@fit$LLH)) {
          write(sprintf("Success with solver = %s", solver), file = log_file, append = TRUE)
          return(fit)
        } else {
          write(sprintf("Solver %s failed to converge or returned bad fit", solver), file = log_file, append = TRUE)
        }
      }
      
      write("All solvers failed. Returning NULL.\n", file = log_file, append = TRUE)
      return(NULL)
    }
    
    # Fit GARCH models in parallel safely
    garch_fits_debug <- mapply(function(x, name) {
      safe_ugarchfit_debug(garch_spec_debug, x, bank = name)
    }, returns_4garch_debug[, price_columns], price_columns, SIMPLIFY = FALSE)
  }
}

# Estimate weights ... using two models to choose from based on the predicted state
rolling_risk_regimes <- function(rvine1, rvine2, returns, returns_test, vix, vix_test, 
                                 price_columns, test_window, rolling_window, N, seed = 7, 
                                 fut_states = NULL, dist_grch, armaorder, 
                                 special_armaorder = NULL, special_banks = NULL,
                                 garchorder = c(1, 1), special_garchorder = c(1, 1), 
                                 var_alpha = 0.01
                                 ){
  # Cutting the observations earlier than the first obs. for the rolling window
  returns <- returns %>% dplyr::select(-"regime")
  returns <- tail(returns, rolling_window)
  
  # Adding the test observations
  returns <- rbind(returns, returns_test)
  # Test_window-number starting indices, for test_window ahead estimation
  start_indices <- seq(1, test_window)
  
  if (is.null(fut_states)){
    message("Running HMM regime rolling prediction...")
    fut_states <- HMM_rolling(vix = vix, vix_test = vix_test, seed = seed, 
                              test_window = test_window, rolling_window = rolling_window)
  } 
  
  message("HMM future regimes passed as an argument.")
  
  message("Initializing cluster with ", detectCores() - 1, " cores...")
  cl <- makeCluster(detectCores() - 1)
  on.exit({
    if (!is.null(cl)) stopCluster(cl)
  }, add = TRUE)
  
  message("Exporting variables to cluster...")
  clusterExport(cl, varlist = c("simulated_logret", "VaR", "ES", "min_var_weights", "max_ir_weights"))
  clusterExport(cl, varlist = c("rvine1", "rvine2", "returns", "returns_test", "price_columns", 
                                "rolling_window", "dist_grch", "armaorder", "special_armaorder", "special_banks", 
                                "fut_states", "garchorder", "special_garchorder", "seed", "var_alpha"), envir = environment())
  clusterEvalQ(cl, {
    library(rugarch)
    library(VineCopula)
    library(dplyr)
  })
  
  helper <- function(start_idx){
    log_file <- "log.txt"
    cat(sprintf("Running start_idx = %s\n", start_idx), file = log_file, append = TRUE)
    tryCatch({
      future_state <- fut_states[start_idx]
      if (future_state==1){
        # Generate from the first state RVine
        sim_returns <- simulated_logret(start_idx = start_idx, rolling = TRUE, N = N, 
                                        rvine = rvine1, price_columns = price_columns, returns = returns, 
                                        rolling_window = rolling_window, log_file = log_file, 
                                        dist_grch = dist_grch, armaorder = armaorder, special_banks = special_banks,
                                        special_armaorder = special_armaorder, seed = seed, special_garchorder = special_garchorder, 
                                        garchorder = garchorder)
      } else if (future_state==2){
        # Generate from the second state RVine
        sim_returns <- simulated_logret(start_idx = start_idx, rolling = TRUE, N = N, 
                                        rvine = rvine2, price_columns = price_columns, returns = returns, 
                                        rolling_window = rolling_window, log_file = log_file, 
                                        dist_grch = dist_grch, armaorder = armaorder, special_banks = special_banks,
                                        special_armaorder = special_armaorder, seed = seed, special_garchorder = special_garchorder, 
                                        garchorder = garchorder)
      }
      
      weights_eq <- rep(1 / ncol(sim_returns), ncol(sim_returns))
      VaR_eq <- VaR(sim_returns, weights_eq, var_alpha)
      ES_eq <- ES(sim_returns, weights_eq, var_alpha)
      list(list(weights_eq = weights_eq, VaR_eq = VaR_eq, ES_eq = ES_eq), 
           min_var_weights(sim_returns, var_alpha), max_ir_weights(sim_returns, var_alpha))
    }, error = function(e) {
      error_msg <- sprintf("Error at start_idx = %s: %s", start_idx, e$message)
      message(error_msg)
      cat(error_msg, file = log_file, append = TRUE, sep = "\n")
      return(NULL)
    })
  }
  message("Starting parallel simulation with pblapply...")
  results <- pblapply(start_indices, helper, cl = cl)
  message("Rolling regime-based risk analysis complete.")
  results
}

# Estimate weights ... using rolling vines
rolling_risk_rollrvines <- function(rvine_list, returns, returns_test, 
                                    price_columns, test_window, rolling_window, N, seed = 7, 
                                    dist_grch, armaorder, special_armaorder = NULL, special_banks = NULL,
                                    garchorder = c(1, 1), special_garchorder = c(1, 1), 
                                    var_alpha = 0.01){
  # Cutting the observations earlier than the first obs. for the rolling window
  returns <- returns %>% dplyr::select(-"regime")
  returns <- tail(returns, rolling_window)
  
  # Adding the test observations
  returns <- rbind(returns, returns_test)
  # Test_window-number starting indices, for test_window ahead estimation
  start_indices <- seq(1, test_window)
  
  # Which rvine from rvine_list to use?
  dif_rvines <- length(rvine_list)
  which_rvine_to_use <- rep(1:dif_rvines, each = test_window/dif_rvines)
  
  message("Initializing cluster with ", detectCores() - 1, " cores...")
  cl <- makeCluster(detectCores() - 1)
  on.exit({
    if (!is.null(cl)) stopCluster(cl)
  }, add = TRUE)
  
  message("Exporting variables to cluster...")
  clusterExport(cl, varlist = c("simulated_logret", "VaR", "ES", "min_var_weights", "max_ir_weights"))
  clusterExport(cl, varlist = c("rvine_list", "returns", "returns_test", "price_columns", 
                                "rolling_window", "dist_grch", "armaorder", "special_armaorder", "special_banks",
                                "garchorder", "special_garchorder", "seed", "var_alpha"), envir = environment())
  clusterEvalQ(cl, {
    library(rugarch)
    library(VineCopula)
    library(dplyr)
  })
  
  helper <- function(start_idx){
    log_file <- "log.txt"
    cat(sprintf("Running start_idx = %s\n", start_idx), file = log_file, append = TRUE)
    tryCatch({
      n_rvine_from_list <- which_rvine_to_use[start_idx]
      rvine_to_use <- rvine_list[[n_rvine_from_list]]
      sim_returns <- simulated_logret(start_idx = start_idx, rolling = TRUE, N = N, seed = seed,
                                      rvine = rvine_to_use, price_columns = price_columns, returns = returns, 
                                      rolling_window = rolling_window, log_file = log_file,
                                      dist_grch = dist_grch, armaorder = armaorder, 
                                      special_armaorder = special_armaorder, special_banks = special_banks,
                                      garchorder = garchorder, special_garchorder = special_garchorder)
      
      weights_eq <- rep(1 / ncol(sim_returns), ncol(sim_returns))
      VaR_eq <- VaR(sim_returns, weights_eq, var_alpha)
      ES_eq <- ES(sim_returns, weights_eq, var_alpha)
      list(list(weights_eq = weights_eq, VaR_eq = VaR_eq, ES_eq = ES_eq), 
           min_var_weights(sim_returns, var_alpha), max_ir_weights(sim_returns, var_alpha))
    }, error = function(e) {
      error_msg <- sprintf("Error at start_idx = %s: %s", start_idx, e$message)
      message(error_msg)
      cat(error_msg, file = log_file, append = TRUE, sep = "\n")
      return(NULL)
    })
  }
  message("Starting parallel simulation with pblapply...")
  results <- pblapply(start_indices, helper, cl = cl)
  message("Rolling R-vine risk simulation complete.")
  results
}


# Portfolio Weights Summary -----------------------------------------------

compute_weights_and_risk <- function(sim_returns, alpha = 0.01) {
  minvar <- min_var_weights(sim_returns, alpha)
  maxir  <- max_ir_weights(sim_returns, alpha)
  
  list(
    weights_minvar = minvar$weights_mvar,
    weights_maxir  = maxir$weights_ir,
    VaR_minvar     = minvar$VaR_mvar,
    VaR_maxir      = maxir$VaR_ir,
    ES_minvar      = minvar$ES_mvar,
    ES_maxir       = maxir$ES_ir
  )
}

# normalise weights
norm_near_zero <- function(col_vector, threshold = 1e-10){
  col_vector <- sapply(col_vector, function(x) ifelse(abs(x) < threshold, 0, x))
  extra_weight <- 1 - sum(col_vector)
  if (extra_weight > 0) {
    col_vector <- col_vector/sum(col_vector)
  } 
  col_vector
}


# Weights: Visualisations -------------------------------------------------

# Plot weights depending on the strategy
plot_weights5000 <- function(port_weights_5000, sim_returns, plot_folder, data_type, 
                             show_title = TRUE, 
                             legend_text_size = 12) {
  message("Creating weights5000 plot...")
  
  weights_plot <- ggplot(port_weights_5000 %>%
                           pivot_longer(cols = starts_with("Weights_"),
                                        names_to = "Strategy",
                                        values_to = "Weight") %>% 
                           filter(Strategy != "Weights_eq") %>% 
                           mutate(Type = case_when(
                             str_detect(Strategy, "minvar") ~ "MinVaR",
                             str_detect(Strategy, "maxir") ~ "MaxIR",
                             TRUE ~ "Other"
                           ),
                           RVine = case_when(
                             str_detect(Strategy, "single") ~ "Single",
                             str_detect(Strategy, "crisis") ~ "Crisis",
                             str_detect(Strategy, "tranq") ~ "Tranquil",
                             str_detect(Strategy, "cvine") ~ "C-Vine",
                             TRUE ~ "Unknown"
                           )), 
                         aes(x = Bank, y = Weight, fill = RVine)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    geom_hline(aes(yintercept = 1/ncol(sim_returns), linetype = "Equal Weights")) +
    facet_wrap(~ Type) +
    scale_fill_brewer(palette = "Set1") +
    scale_linetype_manual(name = "", values = c("Equal Weights" = "dashed")) +
    labs(
      title = if (show_title) "Portfolio Weights by Strategy and RVine" else NULL,
      subtitle = if (show_title) "Based on 5000 RVine Simulations" else NULL,
      y = "Portfolio Weight", x = "", fill = "RVine"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(size = legend_text_size + 1),
      legend.text = element_text(size = legend_text_size), 
      strip.text = element_text(size = legend_text_size)
    ) + 
    ylim(0, 1)
  
  ggsave(paste0(plot_folder, data_type, "_weights5000.png"), plot = weights_plot,
         width = 8, height = 5, dpi = 200, bg = "transparent")
  
  message("weights5000 plot saved.")
  weights_plot
}

# VaR barchart
plot_var5000 <- function(var_5000, plot_folder, data_type, 
                         show_title = TRUE,
                         legend_text_size = 12) {
  message("Creating var5000 plot...")
  
  var_plot <- ggplot(-var_5000 %>%
                       pivot_longer(everything(), names_to = "Strategy_Regime", values_to = "VaR") %>%
                       mutate(
                         Strategy = case_when(
                           grepl("eq", Strategy_Regime) ~ "Equal Weights",
                           grepl("minvar", Strategy_Regime) ~ "MinVaR",
                           grepl("maxir", Strategy_Regime) ~ "MaxIR",
                           TRUE ~ "Other"
                         ),
                         RVine = case_when(
                           grepl("single", Strategy_Regime) ~ "Single",
                           grepl("crisis", Strategy_Regime) ~ "Crisis",
                           grepl("tranq", Strategy_Regime) ~ "Tranquil",
                           grepl("cvine", Strategy_Regime) ~ "C-Vine",
                           TRUE ~ "Unknown"
                         )
                       ), aes(x = Strategy, y = VaR, fill = RVine)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = if (show_title) "VaR by Strategy and Regime" else NULL,
      x = "", y = "VaR (1% quantile)", fill = "RVine"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(hjust = 1, size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = legend_text_size + 1),
      legend.text = element_text(size = legend_text_size), 
      strip.text = element_text(size = legend_text_size)
    )
  
  ggsave(paste0(plot_folder, data_type, "_var5000.png"), plot = var_plot,
         width = 8, height = 5, dpi = 200, bg = "transparent")
  
  message("var5000 plot saved.")
  var_plot
}


export_sim5000_table <- function(var_5000, es_5000, table_folder, data_type) {
  message("Creating and exporting simulation summary table...")
  
  sim5000_table <- cbind(var_5000 %>%
                           pivot_longer(everything(), names_to = "Strategy_Regime", values_to = "VaR"),
                         es_5000 %>%
                           pivot_longer(everything(), names_to = "Strategy_Regime", values_to = "ES") %>%
                           dplyr::select(ES)
  ) %>%
    mutate(
      Strategy = case_when(
        grepl("eq", Strategy_Regime) ~ "Equal Weights",
        grepl("minvar", Strategy_Regime) ~ "MinVar",
        grepl("maxir", Strategy_Regime) ~ "MaxIR",
        TRUE ~ "Other"
      ),
      RVine = case_when(
        grepl("single", Strategy_Regime) ~ "Full",
        grepl("crisis", Strategy_Regime) ~ "Crisis",
        grepl("tranq", Strategy_Regime) ~ "Tranquil",
        grepl("cvine", Strategy_Regime) ~ "C-Vine",
        TRUE ~ "Unknown"
      )
    ) %>%
    dplyr::select(RVine, Strategy, VaR, ES)
  
  sim5000_table_tex <- kable(sim5000_table, format = "latex", booktabs = TRUE, caption = "Simulation Results (5000)", digits = 3)
  writeLines(sim5000_table_tex, paste0(table_folder, data_type, "_sim5000_table.tex"))
  
  message("Simulation table exported to LaTeX.")
  sim5000_table
}

# U: Diagnostics ----------------------------------------------------------

qqplot_u <- function(bank, vine_model, seed = 7){
  set.seed(seed)
  u_sim_single <- RVineSim(N = nrow(u), rvine_single)
  qqplot(u[, bank], u_sim_single[, bank], 
         main = paste("QQ Plot: U of", bank), xlab = "Observed", ylab = "Simulated")
  abline(0, 1)
}

# CoVaR -------------------------------------------------------------------
# Compute simulated CoVaR from a vine modelxs
CoVaRcompute <- function(vine_model, price_columns, n_sim = 5e4, alpha = 0.05, eps = 0.005, seed = 7) {
  d <- ncol(u)
  
  # Simulate from R-vine
  set.seed(seed)
  sim_u <- RVineSim(n_sim, vine_model)
  
  # Create index pairs (i ≠ j)
  index_pairs <- which(!diag(d), arr.ind = TRUE)
  
  # Parallelized CoVaR computation
  covar_list <- pblapply(1:nrow(index_pairs), function(k) {
    i <- index_pairs[k, 1]
    j <- index_pairs[k, 2]
    
    u_i <- sim_u[, i]
    u_j <- sim_u[, j]
    
    # Distress condition: u_j ≤ alpha
    distress_sample <- u_i[u_j <= alpha]
    CoVaR_ij <- quantile(distress_sample, alpha, na.rm = TRUE)
    
    # Median condition: u_j ≈ 0.5 ± eps
    normal_sample <- u_i[abs(u_j - 0.5) < eps]
    CoVaR_median_ij <- quantile(normal_sample, alpha, na.rm = TRUE)
    
    list(i = i, j = j, CoVaR = CoVaR_ij, CoVaR_median_ij  = CoVaR_median_ij, 
         Delta_CoVaR = CoVaR_ij - CoVaR_median_ij)
  })
  
  # Fill result matrices
  CoVaR_mat <- matrix(NA, d, d, dimnames = list(price_columns, price_columns))
  Delta_mat <- matrix(NA, d, d, dimnames = list(price_columns, price_columns))
  CoVaR_med <- matrix(NA, d, d, dimnames = list(price_columns, price_columns))
  
  for (res in covar_list) {
    CoVaR_mat[res$i, res$j] <- res$CoVaR
    Delta_mat[res$i, res$j] <- res$Delta_CoVaR
    CoVaR_med[res$i, res$j] <- res$CoVaR_median_ij
  }
  
  list(CoVaR = CoVaR_mat, CoVaR_med = CoVaR_med, Delta_CoVaR = Delta_mat)
}

CoVaRlogret <- function(CoVaRmatrix, std_residuals, garch_fits){
  inv_cdfs <- lapply(std_residuals, function(x) {
    ecdf_x <- ecdf(x)
    x_vals <- sort(x)
    u_vals <- ecdf_x(x_vals)
    approxfun(u_vals, x_vals, rule = 2, ties = "mean")  # extrapolate at edges
  })
  inv_cdfs <- setNames(inv_cdfs, colnames(std_residuals))
  
  CoVaR_z <- as.data.frame(mapply(function(colname, u_col) {
    inv_cdfs[[colname]](u_col)
  }, colname = colnames(CoVaRmatrix), u_col = as.data.frame(CoVaRmatrix), SIMPLIFY = FALSE))
  
  garch_estimates <- lapply(garch_fits, function(fit) {
    list(sigma = tail(sigma(fit), 1), mean = tail(fitted(fit), 1))
  })
  
  sigmas <- sapply(garch_estimates, `[[`, "sigma")
  means <- sapply(garch_estimates, `[[`, "mean")
  
  CoVaR_logret <- sweep(CoVaR_z, 2, sigmas, `*`)
  CoVaR_logret <- sweep(CoVaR_logret, 2, means, `+`)
  
  CoVaR_logret <- as.matrix(CoVaR_logret)
  dimnames(CoVaR_logret) = list(price_columns, price_columns)
  CoVaR_logret
}

# Plot delta CoVaR
CoVaRplot <- function(CoVaR_matrix, subt = NULL, lower_limit = NULL, show_title = T, 
                      legend_text_size = 12){
  
  lower_limit <- if (is.null(lower_limit)) -0.25 else lower_limit
  CoVaR_matrix_long <- CoVaR_matrix %>%
    as.data.frame() %>% 
    rownames_to_column("Affected") %>%
    pivot_longer(-Affected, names_to = "Distressed", values_to = "DeltaCoVaR")
  
  ggplot(CoVaR_matrix_long, aes(x = Distressed, y = Affected, fill = DeltaCoVaR)) +
    geom_tile(color = "white") + 
    geom_text(aes(label = round(DeltaCoVaR, 2)), color = "black", size = 3) + 
    scale_fill_gradient(
      high = "#f7fbff", 
      low = "royalblue3",
      limits = c(lower_limit, 0),
      name = "ΔCoVaR"
    ) + 
    labs(title = if (show_title) expression(Delta*CoVaR ~ "Heatmap") else NULL,
         x = "Distressed Institution",
         y = "Affected Institution",
         subtitle = subt) +
    theme_minimal(base_size = 10) +
    theme(
      legend.title = element_text(size = legend_text_size + 3),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(2, "lines"),
      axis.text = element_text(size = 12), 
      axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plot delta in delta CoVaR
CoVaR_Delta_in_Delta <- function(CoVaR_matrix, lower_limit = NULL, 
                                 show_title = T, 
                                 legend_text_size = 12){
  lower_limit <- if (is.null(lower_limit)) -0.08 else lower_limit
  CoVaR_matrix_long <- CoVaR_matrix %>%
    as.data.frame() %>% 
    rownames_to_column("Affected") %>%
    pivot_longer(-Affected, names_to = "Distressed", values_to = "DeltaCoVaR")
  
  ggplot(CoVaR_matrix_long, aes(x = Distressed, y = Affected, fill = DeltaCoVaR)) +
    geom_tile(color = "white") + 
    geom_text(aes(label = round(DeltaCoVaR, 2)), color = "black", size = 3) + 
    scale_fill_gradient2(
      high = "firebrick3",
      mid = "#f7fbff",
      low = "royalblue3",
      limits = c(lower_limit, -lower_limit),
      name = "Crisis - Tranq"
    ) + 
    labs(title = if (show_title) expression("Change in " ~ Delta*CoVaR ~ "by Regime") else NULL,
         x = "Distressed Institution",
         y = "Affected Institution") +
    theme_minimal(base_size = 13) +
    theme(
      legend.title = element_text(size = legend_text_size + 5),
      legend.text = element_text(size = legend_text_size),
      axis.text = element_text(size = 12), 
      axis.text.x = element_text(angle = 45, hjust = 1))
}

# Save a CoVaR-based plot
save_covar_plot <- function(covar_data, title, filename, folder, lower_limit = NULL, width = 8, 
                            height = 5, dpi = 200, deltaindelta = FALSE, 
                            show_title = TRUE, 
                            legend_text_size = 12) {
  folder <- sub("[/\\\\]+$", "", folder)
  if (deltaindelta){
    plot <- CoVaR_Delta_in_Delta(covar_data, lower_limit = lower_limit, 
                                 show_title = show_title, legend_text_size = legend_text_size)
  } else {
    plot <- CoVaRplot(covar_data, subt = title, lower_limit = lower_limit,
                      show_title = show_title, legend_text_size = legend_text_size)
  }
  ggsave(filename = file.path(folder, filename), plot = plot,
         width = width, height = height, dpi = dpi, bg = "transparent")
  message("Saved CoVaR-based plot to: ", file.path(folder, filename))
  plot
}

CoVaR_table_summary <- function(CoVaR_matrix, CoVaR_matrix_logr){
  
  # Extract the conditional sigma and mean from each GARCH fit
  garch_params <- lapply(garch_fits, function(fit) {
    list(sigma = tail(sigma(fit), 1),
         mean = tail(fitted(fit), 1))
  })
  
  # Create a tidy table
  CoVaR_table <- expand.grid(
    Affected = rownames(CoVaR_matrix),
    Distressed = colnames(CoVaR_matrix),
    stringsAsFactors = FALSE
  )
  
  # Add ΔCoVaR values
  CoVaR_table <- CoVaR_table %>%
    filter(Affected != Distressed) %>%
    mutate(
      DeltaCoVaR_U = mapply(function(i, j) CoVaR_matrix[i, j], Affected, Distressed),
      DeltaCoVaR_Logret = mapply(function(i, j) CoVaR_matrix_logr[i, j], Affected, Distressed),
      Sigma_Affected = sapply(Affected, function(i) garch_params[[i]]$sigma),
      Mean_Affected = sapply(Affected, function(i) garch_params[[i]]$mean),
      Sigma_Distressed = sapply(Distressed, function(i) garch_params[[i]]$sigma),
      Mean_Distressed = sapply(Distressed, function(i) garch_params[[i]]$mean)
    )
  
  CoVaR_table <- CoVaR_table %>%
    arrange(desc(abs(DeltaCoVaR_U)))
  CoVaR_table
}

# Standard Correlation Techniques -----------------------------------------

# Quad programming for weights based on correlation
minvar_corr <- function(start_idx, returns, price_columns, rolling_window = NULL, rolling = T){
  if (rolling == TRUE){
    data <- returns[start_idx:(start_idx + rolling_window - 1), price_columns]
  } else {
    data <- returns[,  price_columns]
    message("Estimating covariance using full data (non-rolling)")
  }
  
  Sigma <- cov(data)
  
  # Objective function: 0.5 * w' Σ w
  Dmat <- 2 * Sigma
  n_assets <- length(price_columns)
  dvec <- rep(0, n_assets)
  
  # Constraints: sum(w) == 1 and w >= 0
  Amat <- cbind(rep(1, n_assets), diag(n_assets))  # (n_assets + 1) constraints
  bvec <- c(1, rep(0, n_assets))
  result <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  result$solution
}

# Rolling-window optimisation using correlations
minvar_corr_rolling <- function(returns, returns_test, rolling_window, 
                                price_columns, test_window){
  # Cutting the observations earlier than the first obs. for the rolling window
  returns <- returns %>% dplyr::select(-"regime")
  returns <- tail(returns, rolling_window)
  
  # Adding the test observations
  returns <- rbind(returns, returns_test)
  # Test_window-number starting indices, for test_window ahead estimation
  start_indices <- seq(1, test_window)
  helper <- function(x) minvar_corr(x, returns = returns, price_columns = price_columns, 
                                    rolling_window = rolling_window, rolling = T)
  pblapply(start_indices, helper)
}


# Big Summary -------------------------------------------------------------

# Store VaR and ES
var_es_backtest_table <- function(roll_dataframe, var = NULL, es = NULL){
  if (!is.null(var) && var ==TRUE){
    data.frame(Date = returns_backtest$Date, 
               Var_eq = unlist(t(sapply(roll_dataframe, "[[", 1))[, 'VaR_eq']),
               Var_mvar = unlist(t(sapply(roll_dataframe, "[[", 2))[, 'VaR_mvar']),
               Var_ir = unlist(t(sapply(roll_dataframe, "[[", 3))[, 'VaR_ir']))
  } else if (!is.null(es) && es == TRUE){
    data.frame(Date = returns_backtest$Date, 
               ES_eq = unlist(t(sapply(roll_dataframe, "[[", 1))[, 'ES_eq']),
               ES_mvar = unlist(t(sapply(roll_dataframe, "[[", 2))[, 'ES_mvar']),
               ES_ir = unlist(t(sapply(roll_dataframe, "[[", 3))[, 'ES_ir']))
  }
}

# Store weights
weights_backtest <- function(returns_backtest, roll_risk, roll_risk_regimes, 
                             roll_risk_rollvines, roll_risk_cvine, roll_minvar_corr) {
  
  extract_weights <- function(obj, idx, kind) {
    do.call(rbind, lapply(obj, function(x) x[[idx]][[kind]]))
  }
  
  weights_eq         <- extract_weights(roll_risk_regimes, 1, 'weights_eq')
  weights_mvar_single <- extract_weights(roll_risk, 2, 'weights_mvar')
  weights_ir_single   <- extract_weights(roll_risk, 3, 'weights_ir')
  
  weights_mvar_regimes <- extract_weights(roll_risk_regimes, 2, 'weights_mvar')
  weights_ir_regimes   <- extract_weights(roll_risk_regimes, 3, 'weights_ir')
  
  weights_mvar_roll <- extract_weights(roll_risk_rollvines, 2, 'weights_mvar')
  weights_ir_roll   <- extract_weights(roll_risk_rollvines, 3, 'weights_ir')
  
  weights_mvar_cvine <- extract_weights(roll_risk_cvine, 2, 'weights_mvar')
  weights_ir_cvine   <- extract_weights(roll_risk_cvine, 3, 'weights_ir')
  
  weights_minvar_corr <- do.call(rbind, roll_minvar_corr)
  portfolio_returns <- function(returns, weights) rowSums(returns * weights)
  
  data.frame(
    Date = returns_backtest$Date,
    Eq = portfolio_returns(returns_backtest[, -1], weights_eq),
    mvar_single = portfolio_returns(returns_backtest[, -1], weights_mvar_single),
    mvar_regimes = portfolio_returns(returns_backtest[, -1], weights_mvar_regimes),
    ir_single = portfolio_returns(returns_backtest[, -1], weights_ir_single),
    ir_regimes = portfolio_returns(returns_backtest[, -1], weights_ir_regimes),
    mvar_roll = portfolio_returns(returns_backtest[, -1], weights_mvar_roll),
    ir_roll = portfolio_returns(returns_backtest[, -1], weights_ir_roll),
    mvar_cvine = portfolio_returns(returns_backtest[, -1], weights_mvar_cvine),
    ir_cvine = portfolio_returns(returns_backtest[, -1], weights_ir_cvine),
    minvar_corr = portfolio_returns(returns_backtest[, -1], weights_minvar_corr)
  )
}


# VaR Tests ---------------------------------------------------------------

# Run UC and CC tests
run_vartests <- function(portfolio_returns, var_estimates, alpha = 0.01) {
  mapply(function(actual_ret, var_series) {
    VaRTest(alpha = alpha, actual = actual_ret, VaR = var_series)
  }, portfolio_returns, var_estimates, SIMPLIFY = FALSE)
}

save_var_test_table <- function(var_test_df, test_type, data_type, table_folder) {
  message(sprintf("Saving %s VaR test table...", test_type))
  
  var_test_table <- kable(var_test_df, format = "latex", booktabs = TRUE, 
                          caption = paste(test_type, "VaR Test"), digits = 3)
  
  file_path <- paste0(table_folder, data_type, "_vartest_", tolower(test_type), "_table.tex")
  writeLines(var_test_table, file_path)
  
  message(sprintf("%s VaR test table saved to: %s", test_type, file_path))
}

run_durtest <- function(test_obj) {
  data.frame(
    H0 = test_obj$H0,
    LRp = test_obj$LRp,
    decision = test_obj$Decision
  )
}


# Ploting -----------------------------------------------------------------

# Cumulative returns for subset of strategies
plot_cumulative_returns <- function(df, strategy_subset, filename_prefix, title = NULL, subtitle = NULL, 
                                    show_title = TRUE, legend_text_size = 12, legend.position = "right") {
  # Get labels and colors that match the final factor levels
  strategy_subset_labels <- strategy_labels[strategy_subset]
  strategy_subset_colors <- strategy_colors[strategy_subset]
  names(strategy_subset_colors) <- strategy_subset_labels
  
  plot_data <- df %>%
    pivot_longer(cols = -Date, names_to = "Strategy", values_to = "Return") %>%
    filter(Strategy %in% strategy_subset) %>%
    group_by(Strategy) %>%
    mutate(CumulativeReturn = cumsum(Return)) %>%
    ungroup() %>%
    mutate(Strategy = factor(
      Strategy,
      levels = strategy_subset,
      labels = strategy_subset_labels
    ))
  
  p <- ggplot(plot_data, aes(x = Date, y = CumulativeReturn, color = Strategy)) +
    geom_line(linewidth = 0.8) +
    theme_minimal() +
    scale_color_manual(values = strategy_subset_colors) +
    labs(
      title = if (show_title) title else NULL,
      subtitle = if (show_title) subtitle else NULL,
      x = NULL,
      y = "Cumulative Log Return",
      color = "Strategy"
    ) +
    theme(legend.position = legend.position, 
          legend.title = element_text(size = legend_text_size + 1),
          legend.text = element_text(size = legend_text_size),
          legend.key.size = unit(2, "lines"))
  
  ggsave(
    filename = paste0(plot_folder, data_type, "_", filename_prefix, ".png"),
    plot = p,
    width = 8, height = 5, dpi = 200, bg = "transparent"
  )
  message("Saved plot to: ", paste0(plot_folder, data_type, "_", filename_prefix, ".png"))
  p
}

# Plot VaR and ES for EW and all models
plot_var_es <- function(
    df_single, df_roll, df_regimes, df_cvine, portfolio_ret, 
    measure = c("VaR", "ES"), ylims = c(-0.1, 0.08),
    portfolio_type = "Eq", 
    plot_folder, data_type,
    font_size = 14,                  # NEW: default base size
    font_face = "plain",             # NEW: default face
    font_setting = NULL,              # NEW: presets like "bold_large", 
    legend.position = "right", 
    show_title = TRUE
) {
  measure <- match.arg(measure)
  prefix <- tolower(measure)
  
  # Apply font_setting presets
  if (!is.null(font_setting)) {
    if (font_setting == "regular") {
      font_size <- 14
      font_face <- "plain"
    } else if (font_setting == "bold_large") {
      font_size <- 18
      font_face <- "bold"
    } else if (font_setting == "presentation") {
      font_size <- 20
      font_face <- "bold"
    }
  }
  
  col_base <- paste0(ifelse(measure == "VaR", "Var", "ES"), "_", tolower(portfolio_type))
  
  df_plot <- df_single %>%
    dplyr::select(Date, all_of(col_base)) %>%
    mutate(
      Portfolio_ret = portfolio_ret[[portfolio_type]],
      !!paste0(col_base, "_roll") := df_roll[[col_base]],
      !!paste0(col_base, "_regimes") := df_regimes[[col_base]],
      !!paste0(col_base, "_cvine") := df_cvine[[col_base]]
    ) %>%
    pivot_longer(cols = -Date, names_to = "Return_VaR", values_to = "Value") %>%
    mutate(
      Return_VaR = factor(
        Return_VaR,
        levels = c(col_base, paste0(col_base, "_regimes"), paste0(col_base, "_roll"), paste0(col_base, "_cvine"), "Portfolio_ret")
      )
    ) %>%
    arrange(Date, Return_VaR)
  
  labels <- var_es_labels[levels(df_plot$Return_VaR)]
  colors <- var_es_colors[levels(df_plot$Return_VaR)]
  
  p <- ggplot(df_plot, aes(x = Date, y = Value, color = Return_VaR)) +
    geom_line(linewidth = 0.5) +
    theme_minimal() +
    scale_color_manual(values = colors, labels = labels) +
    labs(
      title = if (show_title) "Equal Weights Portfolio" else NULL,
      subtitle = if (show_title) paste0(measure, " 1% for Single, Regime-based, Rolling RVines and C-Vine") else NULL,
      x = "", y = "Log Return", color = "Strategy"
    ) +
    ylim(ylims[1], ylims[2]) +
    theme(
      plot.title = element_text(size = font_size + 8, face = font_face),
      plot.subtitle = element_text(size = font_size + 4, face = font_face),
      legend.title = element_blank(),
      legend.text = element_text(size = font_size, face = font_face),
      legend.key.size = unit(2, "lines"),
      axis.text = element_text(size = font_size),
      legend.position = legend.position 
    )
  
  ggsave(
    paste0(plot_folder, data_type, "_", tolower(portfolio_type), "_", tolower(measure), ".png"),
    plot = p, width = 8, height = 5, dpi = 200, bg = "transparent"
  )
  
  message("Saved ", measure, " plot to: ", plot_folder, data_type, "_", tolower(portfolio_type), "_", tolower(measure), ".png")
  p
}

# Plot VaR & ES for one MVaR strategy
plot_mvar_var_es <- function(
    VaR_df, ES_df, portfolio_ret, portfolio_name, strategy_name, 
    var_color, es_color, plot_folder, data_type, ylims = c(-0.1, 0.08),
    font_size = 14, font_face = "plain", font_setting = NULL, 
    show_title = TRUE, legend.position = "bottom") {
  
  # Apply font presets
  if (!is.null(font_setting)) {
    if (font_setting == "regular") {
      font_size <- 14
      font_face <- "plain"
    } else if (font_setting == "bold_large") {
      font_size <- 18
      font_face <- "bold"
    } else if (font_setting == "presentation") {
      font_size <- 20
      font_face <- "bold"
    }
  }
  
  # Compose variable names
  var_col <- "Var_mvar"
  get_es_col <- function(portf) {
    if (grepl("mvar", portf)) return("ES_mvar")
    if (grepl("ir", portf)) return("ES_ir")
    if (grepl("eq", portf)) return("ES_eq")
    stop("No matching ES column for portfolio_name")
  }
  
  es_col <- get_es_col(portfolio_name)
  portf_col <- paste0("Portf_ret_", portfolio_name)
  
  # Prepare data
  df_plot <- VaR_df %>% 
    dplyr::select(-c(Var_eq, Var_ir)) %>% 
    mutate(
      !!portf_col := portfolio_ret[[portfolio_name]],
      !!es_col := ES_df[[es_col]]
    ) %>% 
    pivot_longer(cols = -Date, names_to = "Return_VaR", values_to = "Value") %>% 
    mutate(
      Return_VaR = factor(Return_VaR, levels = c(var_col, es_col, portf_col))
    ) %>% 
    arrange(Date, Return_VaR)
  
  # Colors and labels for scale
  colors <- c(
    setNames("black", portf_col),
    Var_mvar = var_color,
    setNames(es_color, es_col)
  )
  labels <- c(
    setNames(paste(strategy_name, "Portfolio"), portf_col),
    Var_mvar = paste0("VaR 1% (", strategy_name, ")"),
    setNames(paste0("ES 1% (", strategy_name, ")"), es_col)
  )
  
  # Plot
  p <- ggplot(df_plot, aes(x = Date, y = Value, color = Return_VaR)) +
    geom_line(linewidth = 0.5) +
    theme_minimal() +
    theme(
      legend.position = legend.position,
      legend.title = element_blank(),
      legend.text = element_text(size = font_size, face = font_face),
      axis.text = element_text(size = font_size),
      plot.title = element_text(size = font_size + 6, face = font_face),
      plot.subtitle = element_text(size = font_size + 2, face = font_face),
      axis.title = element_text(size = font_size, face = font_face), 
      legend.key.size = unit(2, "lines")
    ) +
    scale_color_manual(values = colors, labels = labels) +
    labs(
      title = if (show_title) paste(strategy_name, "Portfolio") else NULL,
      subtitle = if (show_title) paste0("Weights optimised for VaR 1% (", strategy_name, ")") else NULL,
      x = "", y = "Log Return"
    ) +
    ylim(ylims[1],  ylims[2])
  
  # Save plot
  filename <- paste0(plot_folder, data_type, "_", portfolio_name, ".png")
  ggsave(filename, plot = p, width = 8, height = 5, dpi = 200, bg = "transparent")
  
  message("Saved plot to: ", filename)
  p
}

# Plot VaR & ES for one IR strategy
plot_ir_var_es <- function(
    VaR_df, ES_df, portfolio_ret, portfolio_name, strategy_name, 
    var_color, es_color, plot_folder, data_type, ylims = c(-0.1, 0.08),
    font_size = 14, font_face = "plain", font_setting = NULL, 
    legend.position = "bottom", show_title = TRUE
) {
  # Apply font presets
  if (!is.null(font_setting)) {
    if (font_setting == "regular") {
      font_size <- 14
      font_face <- "plain"
    } else if (font_setting == "bold_large") {
      font_size <- 18
      font_face <- "bold"
    } else if (font_setting == "presentation") {
      font_size <- 20
      font_face <- "bold"
    }
  }
  
  # Compose variable names
  var_col <- "Var_ir"
  get_es_col <- function(portf) {
    if (grepl("mvar", portf)) return("ES_mvar")
    if (grepl("ir", portf)) return("ES_ir")
    if (grepl("eq", portf)) return("ES_eq")
    stop("No matching ES column for portfolio_name")
  }
  
  es_col <- get_es_col(portfolio_name)
  portf_col <- paste0("Portf_ret_", portfolio_name)
  
  # Prepare data
  df_plot <- VaR_df %>% 
    dplyr::select(-c(Var_eq, Var_mvar)) %>% 
    mutate(
      !!portf_col := portfolio_ret[[portfolio_name]],
      !!es_col := ES_df[[es_col]]
    ) %>% 
    pivot_longer(cols = -Date, names_to = "Return_VaR", values_to = "Value") %>% 
    mutate(
      Return_VaR = factor(Return_VaR, levels = c(var_col, es_col, portf_col))
    ) %>% 
    arrange(Date, Return_VaR)
  
  # Colors and labels
  colors <- c(
    setNames("black", portf_col),
    setNames(var_color, var_col),
    setNames(es_color, es_col)
  )
  
  labels <- c(
    setNames(paste(strategy_name, "Portfolio"), portf_col),
    setNames(paste0("VaR 1% (", strategy_name, ")"), var_col),
    setNames(paste0("ES 1% (", strategy_name, ")"), es_col)
  )
  
  # Plot
  p <- ggplot(df_plot, aes(x = Date, y = Value, color = Return_VaR)) +
    geom_line(linewidth = 0.5) +
    theme_minimal() +
    theme(
      legend.position = legend.position,
      legend.title = element_blank(),
      legend.text = element_text(size = font_size, face = font_face),
      axis.text = element_text(size = font_size),
      plot.title = element_text(size = font_size + 6, face = font_face),
      plot.subtitle = element_text(size = font_size + 2, face = font_face),
      axis.title = element_text(size = font_size, face = font_face), 
      legend.key.size = unit(2, "lines")
    ) +
    scale_color_manual(values = colors, labels = labels) +
    labs(
      title = if (show_title) paste(strategy_name, "Portfolio") else NULL,
      subtitle = if (show_title) paste0("Weights optimised for IR 1% (", strategy_name, ")") else NULL,
      x = "", y = "Log Return"
    ) +
    ylim(ylims[1], ylims[2])
  
  # Save
  filename <- paste0(plot_folder, data_type, "_", portfolio_name, ".png")
  ggsave(filename, plot = p, width = 8, height = 5, dpi = 200, bg = "transparent")
  
  message("Saved plot to: ", filename)
  p
}

# Freq of positive returns
hit_ratio <- function(x) mean(x > 0)
