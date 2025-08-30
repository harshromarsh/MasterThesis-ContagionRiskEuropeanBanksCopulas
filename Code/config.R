# config.R
if (data_type == "CDS") {
  prices_path <- 'Data/banks_cds.csv'
  colname_replace_pattern <- "\\."
  colname_replacement <- " "
  bank_abbr <- function(x) bank_names[names(bank_names) == x]
  long_caption <- "Log CDS Spread Returns"
  
  plot_folder <- "Plots/CDS/"
  table_folder <- "Tables/CDS/"
  
} else if (data_type == "Stock") {
  prices_path <- 'Data/banks_stock_prices.csv'
  colname_replace_pattern <- NULL
  colname_replacement <- NULL
  long_caption <- "Log Stock Returns"
  
  plot_folder <- "Plots/Stock/"
  table_folder <- "Tables/Stock/"
  
} else {
  stop("Unknown data_type")
}

# For the Plots
if (data_type == "Stock"){
  col <- "skyblue"
  type <- 1
} else if (data_type == "CDS"){
  col <- "darkseagreen"
  type <- 1
}

strategy_labels <- c(
  Eq = "Equal Weights",
  mvar_single = "Min VaR (Single RVine)",
  mvar_regimes = "Min VaR (Regime RVine)",
  ir_single = "Max IR (Single RVine)",
  ir_regimes = "Max IR (Regime RVine)", 
  mvar_roll = "Min VaR (Rolling RVine)", 
  ir_roll = "Max IR (Rolling RVine)", 
  mvar_cvine = "Min VaR (CVine)",   
  ir_cvine = "Max IR (CVine)", 
  minvar_corr = "Min Var (Corr)"
)

strategy_colors <- c(
  Eq = "#1F77B4",  
  mvar_single = "#FF7F0E", 
  mvar_regimes = "#2CA02C",
  ir_single = "#D62728", 
  ir_regimes = "#9467BD",
  mvar_roll = "#8C564B",
  ir_roll = "#BCBD22", 
  mvar_cvine = "#17BECF",   
  ir_cvine = "#7F7F7F", 
  minvar_corr = "#E377C2"
)

var_es_colors <- c(
  Portfolio_ret   = "black",
  Var_eq          = "#1F77D9",
  Var_eq_regimes  = "#191970",
  Var_eq_roll     = "#6A5ACD",
  Var_eq_cvine    = "#20B2AA",
  ES_eq           = "#6BA4F7",
  ES_eq_regimes   = "#1A2A8F",
  ES_eq_roll      = "#A89DF2",
  ES_eq_cvine     = "#77D6CF"
)

var_es_labels <- c(
  Portfolio_ret   = "Equal Weights Portfolio",
  Var_eq          = "VaR 1% (RVine Single)",
  Var_eq_regimes  = "VaR 1% (RVine Regimes)",
  Var_eq_roll     = "VaR 1% (RVine Rolling)",
  Var_eq_cvine    = "VaR 1% (CVine)",
  ES_eq           = "ES 1% (RVine Single)",
  ES_eq_regimes   = "ES 1% (RVine Regimes)",
  ES_eq_roll      = "ES 1% (RVine Rolling)",
  ES_eq_cvine     = "ES 1% (CVine)"
)
