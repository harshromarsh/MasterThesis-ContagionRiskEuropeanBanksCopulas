### R-VINE COPULAS ON TIME SERIES OF EUROPEAN BANKS

# Downloading Data -------------------------------------
data_type <- "CDS"
bank_names <- readRDS('Data/bank_names.Rds')

# Backtest
test_size <- 500
var_alpha <- 0.01

source("Code/functions.R")
source("Code/config.R")

prices <- read.csv(prices_path)
colnames(prices) <- gsub(colname_replace_pattern, colname_replacement, colnames(prices))
colnames(prices)[-1] <- sapply(colnames(prices)[-1], bank_abbr)

data_setup(prices)

price_columns <- colnames(prices)[-1]  # exclude the 'Date' column
returns <- calculate_returns(prices = prices, price_columns = price_columns)

message("The CDS spread for ETE and EUROB become exactly the same starting from ", 
      as.character(prices$Date[which(prices$ETE==prices$EUROB)[1]]), "\nThis is likely because the news of EUROB restructuring led to Bloomberg imputing data with those of ETE.")

banks_to_exclude <- "EUROB"
returns <- returns %>% dplyr::select(!any_of(banks_to_exclude))
price_columns <- price_columns[!price_columns %in% banks_to_exclude]
bank_names <- bank_names[!bank_names %in% banks_to_exclude]

# Missing values ----------------------------------------------------------
### Handling missing values: drop some and interpolate values for others
returns_long <- pivot_longer(returns, -Date, names_to = "Ticker", values_to = "Return")
returns_long$Date <- as.Date(returns_long$Date)

# First non-NA values by column
as.Date(sapply(bank_names, function(x) as.Date(returns[which(!is.na(returns[, x]))[1], "Date"]))) %>% sort()

first_date <- as.Date(returns[which(!is.na(returns$DBK))[1], "Date"])
message("The first effective data point will thus be ", as.character(first_date))

missing_values_plot(returns_long) 

message("ACA has a low share of NAs, but almost all of them at the start of observation = excluded")
 
ret_with_missing <- pivot_returns_long(returns_long, first_date)
miss_values <- compute_missing_share(ret_with_missing, returns)
banks_to_exclude <- c(exclude_high_missing(miss_values, 0.1),  "ACA")
returns_clean <- update_returns(returns, ret_with_missing, banks_to_exclude, price_columns, bank_names)

# updated bank names, price columns after exclusion
bank_names <- returns_clean$bank_names
price_columns <- returns_clean$price_columns

returns <- interpolate_missing(returns_clean$returns, data_type, first_date, special_bank = ifelse(data_type == "Stock", "CABK", "DBK"))

#check NA
colSums(is.na(returns))

# Stationarity ------------------------------------------------------------
stationarity_results <- check_stationarity(returns)
save_stationarity_table(stationarity_results, long_caption, table_folder, data_type)

# Regime Detection --------------------------------------------------------
load("Saved Results/regimes.Rds") # loaded into `regimes`
load("Saved Results/vix.Rds") # loaded into `vix`
load("Saved Results/vix_backtest.Rds") # loaded into `vix_backtest`
load("Saved Results/Future_regimes.Rds") # loaded into `fut_states`

# Setting aside Data for Test ---------------------------------------------
split_data <- split_backtest_set(returns, returns, test_days = test_size) # returns as second input, as `vix` is taken from Stock

returns <- split_data$returns
returns_backtest <- split_data$returns_backtest

# Summary Statistics ------------------------------------------------------
summary_stats <- compute_descriptive_stats(returns, returns_backtest, price_columns, long_caption, table_folder, data_type)

# Volatility Clustering
banks_to_plot <- c('BNP', 'DBK', 'RBI', 'SAN')
vol_cl_data <- volatility_clustering(returns, returns_backtest, price_columns, banks_to_plot)

png(paste0(plot_folder, data_type, "_vol_clustering.png"), width = 1600, height = 1000, res = 200, bg = "transparent")
plot(vol_cl_data,
     main = NULL,
     multi.panel = TRUE, major.ticks = "years")
dev.off()
# Regime Fitting ----------------------------------------------------------

# Working out the any misalignments
vix$regime <- regimes
returns <- left_join(returns, vix %>% 
                   filter(Date %in% returns$Date) %>% 
                   dplyr::select(Date, regime), by = "Date") %>% fill(regime)

returns$regime <- as.factor(returns$regime)

fut_states <- left_join(returns_backtest[, 1:2], vix_backtest %>% 
            mutate(regimes = fut_states), by = "Date") %>% fill(regimes) %>% 
  dplyr::select(regimes) 
fut_states <- fut_states[[1]]

# GARCH -------------------------------------------------------------------

garch_result <- fit_garch_models(returns, price_columns, 
                                 armaorder = c(1, 1), 
                                 dist_grch = "sstd", 
                                 special_banks = c("RBI", "EBS", "SEB", "DBK"),
                                 special_armaorder = c(2, 2))

garch_fits <- garch_result$fits
std_residuals <- garch_result$residuals

resid_diagnostics <- check_residual_diagnostics(std_residuals)
resid_diagnostics %>% filter(pval_mean < 0.05 | pval_vol < 0.05)

sapply(std_residuals, sd) # check the standardisation

# GARCH Fit ---------------------------------------------------------------
emp_cdfs <- lapply(std_residuals, ecdf)
emp_cdfs <- setNames(emp_cdfs, colnames(std_residuals))

u <- pit(std_residuals) 

# Fit of GARCH ------------------------------------------------------------

plot_bank_qqplot(garch_fits, std_residuals, data_type, plot_folder, "BNP", qqplot_garch_sstd)
plot_bank_qqplot(garch_fits, std_residuals, data_type, plot_folder, "DBK", qqplot_garch_sstd)
plot_bank_qqplot(garch_fits, std_residuals, data_type, plot_folder, "RBI", qqplot_garch_sstd)
plot_bank_qqplot(garch_fits, std_residuals, data_type, plot_folder, "ETE", qqplot_garch_sstd)

ks_test_results <- pblapply(bank_names, function(bank) {
  kolmogorov_smirnov_test(bank, garch_fits, std_residuals, B = 2000)
})
names(ks_test_results) <- bank_names

save_ks_test_table(ks_test_results, data_type, table_folder)

save_ks_plot(ks_test_results, "SAN", data_type, plot_folder)
save_ks_plot(ks_test_results, "BNP", data_type, plot_folder)

# ECDF and PIT ------------------------------------------------------------
# Prob. Integral Transform
emp_cdfs <- lapply(std_residuals, ecdf)
emp_cdfs <- setNames(emp_cdfs, colnames(std_residuals))

u <- pit(std_residuals) 
save_pairplots(u, data_type = "Stock", plot_folder = plot_folder, banks_per_plot = 5)

# Single whole(train) sample R Vine Copula --------------------------------
rvine_single <- RVineStructureSelect(u, familyset = 0:6)

custom_rvine_plot(rvine_single, tree = 1, type = type, col = col)
save_rvine_tree_plot(rvine_single, tree = 1, data_type = data_type, type = type,
                     col = col, plot_folder = plot_folder, filename_prefix = "rvine_single")

# C-Vine (baseline) -------------------------------------------------------
cvine_single <- RVineStructureSelect(u, type = 1, familyset = 0:6) #type = 1

custom_rvine_plot(cvine_single, tree = 1, type = type, col = col)
save_rvine_tree_plot(cvine_single, tree = 1, data_type = data_type, type = type,
                     col = col, plot_folder = plot_folder, filename_prefix = "cvine_single")

vuong <- RVineVuongTest(u, rvine_single, cvine_single)

# Two-regime R Vine Copulas -----------------------------------------------

length(returns$regime) == nrow(u)
u_regime1 <- u[returns$regime == 1, ]
u_regime2 <- u[returns$regime == 2, ]

rvine_regime1 <- RVineStructureSelect(u_regime1, familyset = 0:6) 
rvine_regime2 <- RVineStructureSelect(u_regime2, familyset = 0:6) 

custom_rvine_plot(rvine_regime1, tree = 1, type = type, col = col)
save_rvine_tree_plot(rvine_regime1, tree = 1, data_type = data_type, type = type,
                     col = col, plot_folder = plot_folder, filename_prefix = "rvine_crisis")

custom_rvine_plot(rvine_regime2, tree = 1, type = type, col = col)
save_rvine_tree_plot(rvine_regime2, tree = 1, data_type = data_type, type = type,
                     col = col, plot_folder = plot_folder, filename_prefix = "rvine_tranq")
# Print-outs --------------------------------------------------------------

save_vine_tree_tables(rvine_single, data_type, "vine_single", table_folder, trees_to_export = 1:2, 
                      tree_to_print = 1)

save_vine_tree_tables(rvine_regime1, data_type, "vine_crisis", table_folder, trees_to_export = 1, 
                      tree_to_print = 1)

save_vine_tree_tables(rvine_regime2, data_type, "vine_tranq", table_folder, trees_to_export = 1, 
                      tree_to_print = 1)

save_vine_tree_tables(cvine_single, data_type, "vine_cvine", table_folder, trees_to_export = 1, 
                      tree_to_print = 1)

comparison_dep <- compare_tree_dependence(rvine_regime1, rvine_regime2)
comparison_dep

# Information Criteria ----------------------------------------------------

res_ic <- compute_ic_values(u, u_regime1, u_regime2,
                            rvine_single, rvine_regime1, rvine_regime2,
                            cvine_single, regimes, returns)
ic_df <- res_ic$ic_df
ic_df
save_ic_table(ic_df, paste0(table_folder, data_type, "_ic_table.tex"))

# Vuong Test ----------------------------------------------------

df1 <- res_ic$df$df1
df2 <- res_ic$df$df2
df_single <- res_ic$df$df_single

ll_single <- res_ic$ll$ll_single
ll_regimes <- numeric(length(ll_single))
ll_regimes[returns$regime == 1] <- res_ic$ll$ll_r1
ll_regimes[returns$regime == 2] <- res_ic$ll$ll_r2

vuong_result <- vuong_test(ll_regimes, ll_single, df1 + df2, df_single)
save_vuong_table(vuong_result, paste0(table_folder, data_type, "_vuong_test_table.tex"))

# Rolling R-Vines ----------------------------------------------------------

### SPECIFY ARMA ORDERS

roll_rvines <- rolling_rvines(u = u, returns = returns, rolling_window = 250, step_size = 250)
roll_rvines_backtest <- rolling_rvines(returns = returns_backtest, rolling_window = 250, step_size = 250, 
                                       compute_u = TRUE, price_columns = price_columns, 
                                       dist_grch = "sstd", armaorder = c(1, 1), 
                                       special_banks = c("RBI", "EBS", "SEB", "DBK"),
                                       special_armaorder = c(2, 2))

#xxx <- fit_garch_models(returns = returns_backtest[251:500, ], dist_grch = "sstd", armaorder = c(1, 1), 
#                 special_banks = c("RBI", "EBS", "SEB", "DBK"), log = T,
#                 special_armaorder = c(2, 2), price_columns = price_columns)

plot_roll_rvine(roll_rvines, indices = 1:4, 
                plot_folder = plot_folder, data_type = data_type, 
                tree = 1, type = type, col = col)

# Simulations -------------------------------------------------------------

simulated_logret(rolling = F, N = 1000, price_columns = price_columns, rvine = rvine_single, returns = returns, 
                 dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
                 special_armaorder = c(2, 2))

load("Saved Results/Future_regimes.Rds")

# Portfolios ---------------------------------------------------------

# weights
sim_returns <- simulated_logret(N = 5000, rvine = rvine_single, returns = returns, 
                                price_columns = price_columns, rolling = FALSE, seed = 7, 
                                dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
                                special_armaorder = c(2, 2))

sim_returns_crisis <- simulated_logret(N = 5000, rvine = rvine_regime1, returns = returns, 
                                       price_columns = price_columns, rolling = FALSE, seed = 7, 
                                       dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
                                       special_armaorder = c(2, 2))
sim_returns_tranquil <- simulated_logret(N = 5000, rvine = rvine_regime2, returns = returns, 
                                         price_columns = price_columns, rolling = FALSE, seed = 7, 
                                         dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
                                         special_armaorder = c(2, 2))
sim_returns_cvine <- simulated_logret(N = 5000, rvine = cvine_single, returns = returns, 
                                      price_columns = price_columns, rolling = FALSE, seed = 7, 
                                      dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
                                      special_armaorder = c(2, 2))

# Summary
vine_models <- list(
  single   = sim_returns,
  crisis   = sim_returns_crisis,
  tranquil = sim_returns_tranquil,
  cvine    = sim_returns_cvine
)

all_static_results <- lapply(vine_models, compute_weights_and_risk, alpha = var_alpha)

port_weights_5000 <- data.frame(
  Bank = bank_names,
  Weights_eq = rep(1 / ncol(sim_returns), ncol(sim_returns))
)

for (name in names(all_static_results)) {
  port_weights_5000[[paste0("Weights_minvar_", name)]] <- all_static_results[[name]]$weights_minvar
  port_weights_5000[[paste0("Weights_maxir_", name)]]  <- all_static_results[[name]]$weights_maxir
}

# Near-zero weights to exactly zero
weight_cols <- unclass(port_weights_5000)[-1] # get weight columns, exclude the Bank col
weight_cols <- lapply(weight_cols, norm_near_zero)

port_weights_5000 <- cbind(port_weights_5000$Bank, as.data.frame(weight_cols))
colnames(port_weights_5000)[1] <- "Bank"

var_5000 <- data.frame(
  Var_eq_single  = VaR(sim_returns, port_weights_5000$Weights_eq, var_alpha),
  Var_eq_crisis  = VaR(sim_returns_crisis, port_weights_5000$Weights_eq, var_alpha),
  Var_eq_tranq   = VaR(sim_returns_tranquil, port_weights_5000$Weights_eq, var_alpha),
  Var_eq_cvine   = VaR(sim_returns_cvine, port_weights_5000$Weights_eq, var_alpha)
)

es_5000 <- data.frame(
  ES_eq_single   = ES(sim_returns, port_weights_5000$Weights_eq, var_alpha),
  ES_eq_crisis   = ES(sim_returns_crisis, port_weights_5000$Weights_eq, var_alpha),
  ES_eq_tranq    = ES(sim_returns_tranquil, port_weights_5000$Weights_eq, var_alpha),
  ES_eq_cvine    = ES(sim_returns_cvine, port_weights_5000$Weights_eq, var_alpha)
)

# Append minvar and maxir VaR/ES
for (name in names(all_static_results)) {
  var_5000[[paste0("Var_minvar_", name)]] <- all_static_results[[name]]$VaR_minvar
  var_5000[[paste0("Var_maxir_", name)]]  <- all_static_results[[name]]$VaR_maxir
  es_5000[[paste0("ES_minvar_", name)]]   <- all_static_results[[name]]$ES_minvar
  es_5000[[paste0("ES_maxir_", name)]]    <- all_static_results[[name]]$ES_maxir
}

# Visualisations

plot_weights5000(port_weights_5000, sim_returns, plot_folder, data_type, 
                 show_title = F, legend_text_size = 14)
plot_var5000(var_5000, plot_folder, data_type, 
             show_title = F, legend_text_size = 14)
export_sim5000_table(var_5000, es_5000, table_folder, data_type)

# QQplots -----------------------------------------------------------------

qqplot_u("EBS", rvine_single)
qqplot_u("DBK", rvine_single)
qqplot_u("SAN", rvine_single)
qqplot_u("BNP", rvine_single)
qqplot_u("RBI", rvine_regime1)

# CoVaR -------------------------------------------------------------------

CoVaR_single <- CoVaRcompute(rvine_single, price_columns = price_columns)
CoVaR_crisis <- CoVaRcompute(rvine_regime1, price_columns = price_columns)
CoVaR_single_logret <- CoVaRlogret(CoVaR_single$CoVaR, std_residuals = std_residuals, garch_fits = garch_fits)
CoVaR_single_med_logret <- CoVaRlogret(CoVaR_single$CoVaR_med, std_residuals = std_residuals, garch_fits = garch_fits)
Delta_CoVaR_logret_single <- CoVaR_single_logret - CoVaR_single_med_logret

CoVaR_tranq <- CoVaRcompute(rvine_regime2, price_columns = price_columns)

# Saving Plots
save_covar_plot(CoVaR_single$Delta_CoVaR, NULL, 
                paste0(data_type, "_covar_single.png"), plot_folder, lower_limit = -0.5, 
                show_title = F)

save_covar_plot(CoVaR_crisis$Delta_CoVaR, "Crisis R-Vine", 
                paste0(data_type, "_covar_crisis.png"), plot_folder, lower_limit = -0.5)

save_covar_plot(CoVaR_tranq$Delta_CoVaR, "Tranq R-Vine", 
                paste0(data_type, "_covar_tranq.png"), plot_folder, lower_limit = -0.5)

save_covar_plot(CoVaR_crisis$Delta_CoVaR - CoVaR_tranq$Delta_CoVaR, 
                "ΔCoVaR: Crisis - Tranquil", 
                paste0(data_type, "_covar_deltaindelta.png"), plot_folder, deltaindelta = T, 
                lower_limit = -0.12)


head(CoVaR_table_summary(CoVaR_single$Delta_CoVaR, Delta_CoVaR_logret_single))
colMeans(CoVaR_single$Delta_CoVaR, na.rm = T) %>% sort()

# Backtest Results --------------------------------------------------------

roll_risk_regimes <- rolling_risk_regimes(rvine1 = rvine_regime1, rvine2 = rvine_regime2, returns = returns, 
                                          returns_test = returns_backtest, price_columns = price_columns, 
                                          vix = vix, vix_test = vix_backtest, test_window = test_size, 
                                          rolling_window = 1000, N = 1000, fut_states = fut_states, 
                                          dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
                                          special_armaorder = c(2, 2))

roll_risk_regimes_cleaned <- carry_forward_results(roll_risk_regimes) # 167
roll_risk_regimes <- roll_risk_regimes_cleaned$cleaned

save(roll_risk_regimes, file = "Saved Results/roll_risk_regimes_cds.RData")
#load("Saved Results/roll_risk_regimes_cds.RData")

roll_risk <- rolling_risk(rvine = rvine_single, returns = returns, returns_test = returns_backtest, 
                          price_columns, test_window = test_size, rolling_window = 1000, N = 1000, 
                          dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
                          special_armaorder = c(2, 2))

roll_risk_cleaned <- carry_forward_results(roll_risk) # 165
roll_risk <- roll_risk_cleaned$cleaned

save(roll_risk, file = "Saved Results/roll_risk_cds.RData")
#load("Saved Results/roll_risk_cds.RData")

cat("There are in total", length(roll_rvines), "RVines estimated in a rolling window manner.\nFor backtest of", test_size, "obs. we use the last RVine and an RVine estimated using first half")

if (test_size == 500) {
  rvine_list <- list(roll_rvines[[3]][[3]], roll_rvines_backtest[[4]][[3]])
} else if (test_size == 250){
  rvine_list <- list(roll_rvines[[5]][[3]])
}

roll_risk_rollvines <- rolling_risk_rollrvines(
  rvine_list = rvine_list, 
  returns = returns, returns_test = returns_backtest, price_columns = price_columns,
  test_window = test_size, rolling_window = 1000, N = 1000, 
  dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
  special_armaorder = c(2, 2))

roll_risk_rollvines_cleaned <- carry_forward_results(roll_risk_rollvines) # 165
roll_risk_rollvines <- roll_risk_rollvines_cleaned$cleaned

save(roll_risk_rollvines, file = "Saved Results/roll_risk_rollvines_cds.RData")
#load("Saved Results/roll_risk_rollvines_cds.RData")
roll_risk_cvine <- rolling_risk(rvine = cvine_single, returns = returns, returns_test = returns_backtest, 
                                price_columns, test_window = test_size, rolling_window = 1000, N = 1000, 
                                dist_grch = "sstd", armaorder = c(1, 1), special_banks = c("RBI", "EBS", "SEB", "DBK"),
                                special_armaorder = c(2, 2))

roll_risk_cvine_cleaned <- carry_forward_results(roll_risk_cvine) # 165
roll_risk_cvine <- roll_risk_cvine_cleaned$cleaned

save(roll_risk_cvine, file = "Saved Results/roll_risk_cvine_cds.RData")
#load("Saved Results/roll_risk_cvine_cds.RData")

# St. correlation things --------------------------------------------------

roll_minvar_corr <- minvar_corr_rolling(returns = returns, returns_test = returns_backtest, 
                                        price_columns = price_columns, test_window = test_size, rolling_window = 1000)

if (is.null(roll_risk[[1]])){
  roll_risk <- roll_risk[-1]
  roll_risk_cvine <- roll_risk_cvine[-1]
  roll_risk_rollvines <- roll_risk_rollvines[-1]
  roll_risk_regimes <- roll_risk_regimes[-1]
  roll_minvar_corr <- roll_minvar_corr[-1]
  returns_backtest <- returns_backtest[-1, ]
}
# Summarising -------------------------------------------------------------

# Portfolio returns, with estimated weights

VaR_backtest_single <- var_es_backtest_table(roll_risk, var = T)
ES_backtest_single <- var_es_backtest_table(roll_risk, es = T)

VaR_backtest_regimes <- var_es_backtest_table(roll_risk_regimes, var = T)
ES_backtest_regimes <- var_es_backtest_table(roll_risk_regimes, es = T)

VaR_backtest_rolling <- var_es_backtest_table(roll_risk_rollvines, var = T)
ES_backtest_rolling <- var_es_backtest_table(roll_risk_rollvines, es = T)

VaR_backtest_cvine <- var_es_backtest_table(roll_risk_cvine, var = T)
ES_backtest_cvine <- var_es_backtest_table(roll_risk_cvine, es = T)

# Portfolio weights (Backtest)

portfolio_ret <- weights_backtest(returns_backtest = returns_backtest, 
                                  roll_risk = roll_risk, 
                                  roll_risk_regimes = roll_risk_regimes, 
                                  roll_risk_rollvines = roll_risk_rollvines, 
                                  roll_risk_cvine = roll_risk_cvine, 
                                  roll_minvar_corr= roll_minvar_corr)

# VaR Tests -----------------------------------------------------------------

# UC & CC Tests -----------------------------------------------------------

return_series <- list(
  "EW, Single"     = portfolio_ret$Eq,
  "EW, Regime"     = portfolio_ret$Eq,
  "EW, Rolling"    = portfolio_ret$Eq,
  "EW, C-Vine"     = portfolio_ret$Eq,
  "MVaR, Single"   = portfolio_ret$mvar_single,
  "MVaR, Regime"   = portfolio_ret$mvar_regimes,
  "MVaR, Rolling"  = portfolio_ret$mvar_roll,
  "MVaR, C-Vine"   = portfolio_ret$mvar_cvine,
  "IR, Single"     = portfolio_ret$ir_single,
  "IR, Regime"     = portfolio_ret$ir_regimes,
  "IR, Rolling"    = portfolio_ret$ir_roll,
  "IR, C-Vine"     = portfolio_ret$ir_cvine, 
  "MinVar-Corr"    = portfolio_ret$minvar_corr
)

var_series <- list(
  "EW, Single"     = VaR_backtest_single$Var_eq,
  "EW, Regime"     = VaR_backtest_regimes$Var_eq,
  "EW, Rolling"    = VaR_backtest_rolling$Var_eq,
  "EW, C-Vine"     = VaR_backtest_cvine$Var_eq,
  "MVaR, Single"   = VaR_backtest_single$Var_mvar,
  "MVaR, Regime"   = VaR_backtest_regimes$Var_mvar,
  "MVaR, Rolling"  = VaR_backtest_rolling$Var_mvar,
  "MVaR, C-Vine"   = VaR_backtest_cvine$Var_mvar,
  "IR, Single"     = VaR_backtest_single$Var_ir,
  "IR, Regime"     = VaR_backtest_regimes$Var_ir,
  "IR, Rolling"    = VaR_backtest_rolling$Var_ir,
  "IR, C-Vine"     = VaR_backtest_cvine$Var_ir
)

# Run the VaR tests
vartests <- run_vartests(return_series[!names(return_series) %in% "MinVar-Corr"], var_series, var_alpha)

vartest_uc <- data.frame(
  Expected = sapply(vartests, function(x) x$expected.exceed),
  Actual   = sapply(vartests, function(x) x$actual.exceed),
  UC_H0    = sapply(vartests, function(x) x$uc.H0),
  P_value  = sapply(vartests, function(x) x$uc.LRp),
  UC_Decision = sapply(vartests, function(x) x$uc.Decision)
)

vartest_cc <- data.frame(
  Expected = sapply(vartests, function(x) x$expected.exceed),
  Actual   = sapply(vartests, function(x) x$actual.exceed),
  CC_H0    = sapply(vartests, function(x) x$cc.H0),
  P_value  = sapply(vartests, function(x) x$cc.LRp),
  CC_Decision = sapply(vartests, function(x) x$cc.Decision)
)
save_var_test_table(vartest_uc, "Unconditional", data_type, table_folder)
save_var_test_table(vartest_cc, "Conditional", data_type, table_folder)

# VaR Duration Test -------------------------------------------------------
# VaRDurTest, (alpha = 0.05)

# Run and collect all tests
vardurtest <- cbind(Strategy = rep(c("Equal Weight", "Min VaR", "Max IR"), each = 4), 
                    Model = rep(c("Single", "Regimes", "Rolling", "C-Vine"), by = 3),
                    rbind(
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$Eq, VaR_backtest_single$Var_eq))),
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$Eq, VaR_backtest_regimes$Var_eq))),
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$Eq, VaR_backtest_rolling$Var_eq))),
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$Eq, VaR_backtest_cvine$Var_eq))),
                      
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$mvar_single, VaR_backtest_single$Var_mvar))),
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$mvar_regimes, VaR_backtest_regimes$Var_mvar))),
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$mvar_roll, VaR_backtest_rolling$Var_mvar))),
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$mvar_cvine, VaR_backtest_cvine$Var_mvar))),
                      
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$ir_single, VaR_backtest_single$Var_ir))),
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$ir_regimes, VaR_backtest_regimes$Var_ir))),
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$ir_roll, VaR_backtest_rolling$Var_ir))), 
                      cbind(run_durtest(VaRDurTest(var_alpha, portfolio_ret$ir_cvine, VaR_backtest_cvine$Var_ir)))
                    )
)

save_var_test_table(vardurtest, "Duration", data_type, table_folder)

# DM test -----------------------------------------------------------------

# Run Diebold-Mariano test (power = 2 = squared loss, but here it's returns)
# Squared error loss
se_cvine <- (portfolio_ret$ir_cvine - VaR_backtest_cvine$Var_ir)^2
se_single <- (portfolio_ret$ir_single - VaR_backtest_single$Var_ir)^2

dm.test(se_cvine, se_single, h = 1, alternative = "two.sided")

# Backtest Plots ----------------------------------------------------------

# Define groups of strategies to plot
minvar_strategies <- c('mvar_single', 'mvar_regimes', 'mvar_roll', 'mvar_cvine', 
                       'minvar_corr')
ir_strategies <- c('Eq', 'ir_single', 'ir_regimes', 'ir_roll', 'ir_cvine')

# Cumulative Returns from Strategies
plot_cumulative_returns(
  df = portfolio_ret,
  strategy_subset = minvar_strategies,
  title = "Cumulative Portfolio Returns",
  subtitle = "Equal and Min VaR Weights",
  filename_prefix = "cumul_ret_backtest_minvar", 
  legend_text_size = 10, 
  show_title = F
)

plot_cumulative_returns(
  df = portfolio_ret,
  strategy_subset = ir_strategies,
  title = "Cumulative Portfolio Returns",
  subtitle = "Equal and Max IR Weights",
  filename_prefix = "cumul_ret_backtest_ir", 
  legend_text_size = 10, 
  show_title = F
)

# EW portfolio wih VaR and ES
plot_var_es(
  df_single = VaR_backtest_single,
  df_roll = VaR_backtest_rolling,
  df_regimes = VaR_backtest_regimes,
  df_cvine = VaR_backtest_cvine,
  portfolio_ret = portfolio_ret,
  measure = "VaR",
  portfolio_type = "Eq",
  plot_folder = plot_folder,
  data_type = data_type, 
  ylims = c(-0.2, 0.1), 
  font_size = 8, font_face = "bold",
  legend.position = "right"
)

plot_var_es(
  df_single = ES_backtest_single,
  df_roll = ES_backtest_rolling,
  df_regimes = ES_backtest_regimes,
  df_cvine = ES_backtest_cvine,
  portfolio_ret = portfolio_ret,
  measure = "ES",
  portfolio_type = "Eq",
  plot_folder = plot_folder,
  data_type = data_type, 
  ylims = c(-0.2, 0.1), 
  font_size = 8, font_face = "bold",
  legend.position = "bottom"
)

plot_mvar_var_es(
  VaR_df = VaR_backtest_single,
  ES_df = ES_backtest_single,
  portfolio_ret = portfolio_ret,
  portfolio_name = "mvar_single",
  strategy_name = "RVine Single",
  var_color = "#FF7F0E",
  es_color = "#FFBC80",
  plot_folder = plot_folder,
  data_type = data_type, 
  ylims = c(-0.3, 0.1), 
  font_setting = "regular"
)

plot_mvar_var_es(
  VaR_df = VaR_backtest_regimes,
  ES_df = ES_backtest_regimes,
  portfolio_ret = portfolio_ret,
  portfolio_name = "mvar_regimes",
  strategy_name = "RVine Regimes",
  var_color = "#2CA02C",
  es_color = "#A9DFBF",
  plot_folder = plot_folder,
  data_type = data_type,  
  ylims = c(-0.3, 0.1), 
  font_setting = "regular"
)

plot_mvar_var_es(
  VaR_df = VaR_backtest_rolling,
  ES_df = ES_backtest_rolling,
  portfolio_ret = portfolio_ret,
  portfolio_name = "mvar_roll",
  strategy_name = "RVine Rolling",
  var_color = "#8C564B",
  es_color = "#D7B8AE",
  plot_folder = plot_folder,
  data_type = data_type, 
  ylims = c(-0.3, 0.1), 
  font_setting = "regular"
)

plot_mvar_var_es(
  VaR_df = VaR_backtest_cvine,
  ES_df = ES_backtest_cvine,
  portfolio_ret = portfolio_ret,
  portfolio_name = "mvar_cvine",
  strategy_name = "CVine",
  var_color = "#17BECF",
  es_color = "#BADFF7",
  plot_folder = plot_folder,
  data_type = data_type,
  ylims = c(-0.3, 0.1), 
  font_setting = "regular"
)

# IR Portfolios

plot_ir_var_es(
  VaR_df = VaR_backtest_single,
  ES_df = ES_backtest_single,
  portfolio_ret = portfolio_ret,
  portfolio_name = "ir_single",
  strategy_name = "RVine Single",
  var_color = "#D62728",
  es_color = "#F4A6A6",
  plot_folder = plot_folder,
  data_type = data_type, 
  ylims = c(-0.3, 0.1), 
  font_setting = "regular"
)

plot_ir_var_es(
  VaR_df = VaR_backtest_regimes,
  ES_df = ES_backtest_regimes,
  portfolio_ret = portfolio_ret,
  portfolio_name = "ir_regimes",
  strategy_name = "RVine Regimes",
  var_color = "#9467BD",
  es_color = "#D4B9DA",
  plot_folder = plot_folder,
  data_type = data_type, 
  ylims = c(-0.3, 0.1), 
  font_setting = "regular"
)

plot_ir_var_es(
  VaR_df = VaR_backtest_rolling,
  ES_df = ES_backtest_rolling,
  portfolio_ret = portfolio_ret,
  portfolio_name = "ir_roll",
  strategy_name = "RVine Rolling",
  var_color = "#BCBD22",
  es_color = "#E4E6A1",
  plot_folder = plot_folder,
  data_type = data_type, 
  ylims = c(-0.3, 0.1), 
  font_setting = "regular"
)

plot_ir_var_es(
  VaR_df = VaR_backtest_cvine,
  ES_df = ES_backtest_cvine,
  portfolio_ret = portfolio_ret,
  portfolio_name = "ir_cvine",
  strategy_name = "CVine",
  var_color = "#7F7F7F",
  es_color = "lightgrey",
  plot_folder = plot_folder,
  data_type = data_type, 
  ylims = c(-0.3, 0.1), 
  font_setting = "regular"
)

# Portfolio Performance ---------------------------------------------------
xts_portfolio_returns <- xts(portfolio_ret[, -which(names(portfolio_ret) == "Date")],
                             order.by = portfolio_ret$Date)

# Combine all into one summary table
summary_table <- data.frame(
  AnnualizedReturn = apply(xts_portfolio_returns, 2, function(x) Return.annualized(x, scale = 252, geometric = FALSE)),
  AnnualizedVolatility = apply(xts_portfolio_returns, 2, function(x) StdDev.annualized(x, scale = 252, geometric = FALSE)),
  SharpeRatio = apply(xts_portfolio_returns, 2, function(x) SharpeRatio.annualized(x, Rf = 0, scale = 252, geometric = FALSE)),
  SortinoRatio = apply(xts_portfolio_returns, 2, function(x) SortinoRatio(x, MAR = 0)),
  MaxDrawdown = apply(xts_portfolio_returns, 2, function(x) maxDrawdown(x)),
  HitRatio = apply(xts_portfolio_returns, 2, hit_ratio),
  VaR = apply(xts_portfolio_returns, 2, function(x) quantile(x, probs = var_alpha, na.rm = TRUE))
)

summary_table
summ_strategy_table <- kable(summary_table, format = "latex", booktabs = TRUE, caption = "Portfolio Performance Metrics", 
                             digits = 3)
writeLines(summ_strategy_table, paste0(table_folder, data_type, "_summ_strategy_table.tex"))

# --- Stress Periods ---
stress_windows <- list(
  "Bond market stress" = "2023-07-01/2023-10-31",
  "Credit jitters"     = "2024-03-01/2024-03-31",
  "Geopolitical spike" = "2024-04-01/2024-05-15",
  "Stagflation shock"  = "2025-01-01/2025-02-29",
  "Yield surprise"     = "2025-05-01/2025-06-15"
)

for (name in names(stress_windows)) {
  range <- stress_windows[[name]]
  
  cat("\n---", name, "---\n")
  cut_xts_portfolio_returns <- xts_portfolio_returns[range]
  print(data.frame(
    AnnualizedReturn = apply(cut_xts_portfolio_returns, 2, function(x) Return.annualized(x, scale = 252, geometric = FALSE)),
    AnnualizedVolatility = apply(cut_xts_portfolio_returns, 2, function(x) StdDev.annualized(x, scale = 252, geometric = FALSE)),
    SharpeRatio = apply(cut_xts_portfolio_returns, 2, function(x) SharpeRatio.annualized(x, Rf = 0, scale = 252, geometric = FALSE)),
    MaxDrawdown = apply(cut_xts_portfolio_returns, 2, function(x) maxDrawdown(x)),
    HitRatio = apply(cut_xts_portfolio_returns, 2, hit_ratio),
    VaR = apply(cut_xts_portfolio_returns, 2, function(x) quantile(x, probs = var_alpha, na.rm = TRUE))
  ))
}

for (name in names(stress_windows)) {
  range <- stress_windows[[name]]
  
  r1 <- xts_portfolio_returns$mvar_regimes[range]
  r2 <- xts_portfolio_returns$mvar_single[range]
  VaR_backtest_regimes_xts <- xts(VaR_backtest_regimes[, -which(names(VaR_backtest_regimes) == "Date")],
                                  order.by = VaR_backtest_regimes$Date)
  VaR_backtest_single_xts <- xts(VaR_backtest_single[, -which(names(VaR_backtest_single) == "Date")],
                                 order.by = VaR_backtest_single$Date)
  v1 <- VaR_backtest_regimes_xts$Var_mvar[range]
  v2 <- VaR_backtest_single_xts$Var_mvar[range]
  
  se1 <- (r1 - v1)^2
  se2 <- (r2 - v2)^2
  
  cat("\n---", name, "---\n")
  print(dm.test(se1, se2, h = 1, alternative = "less"))
}

