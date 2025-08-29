### Downloading Stock Prices of Systemically Relevant Banks

# Libraries -------------------------------------------------------------
packages <- c('quantmod', 'dplyr', 'purrr', 'zoo')
instld <- sapply(packages, require, character.only = T)
if (sum(instld)<length(instld)) install.packages(packages[!instld])
invisible(lapply(packages, library, character.only = TRUE))
rm(packages, instld)

# Downloading Stock Prices ------------------------------------------------
# Parent banks and their tickers
parent_banks <- c(
  "Erste Group" = "EBS.VI",
  "Raiffeisen Bank Intl" = "RBI.VI",
  "UniCredit" = "UCG.MI",
  "BNP Paribas" = "BNP.PA",
  "KBC Group" = "KBC.BR",
  "ING Group" = "INGA.AS",
  "Deutsche Bank" = "DBK.DE",
  "Commerzbank" = "CBK.DE",
  "Santander" = "SAN.MC",
  "BBVA" = "BBVA.MC",
  "CaixaBank" = "CABK.MC",
  "Societe Generale" = "GLE.PA",
  "Credit Agricole" = "ACA.PA",
  "HSBC" = "HSBA.L",
  "Eurobank" = "EUROB.AT",
  "National Bank of Greece" = "ETE.AT",
  #"OTP Bank" = "OTP.BD",
  "Intesa Sanpaolo" = "ISP.MI",
  "Nordea" = "NDA-SE.ST",
  "SEB" = "SEB-A.ST",
  "Swedbank" = "SWED-A.ST",
  "Danske Bank" = "DANSKE.CO"#,
  #"Banca Transilvania" = "TLV.RO"
)

# Download adjusted close prices (last 5 years)
get_bank_stockprices <- function(ticker) {
  tryCatch({
    data <- getSymbols(ticker, src = "yahoo", from = "2000-01-01", 
                       to = "2025-05-24", auto.assign = FALSE)
    Ad(data) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column("Date") %>% 
      rename(!!ticker := 2) #!!: unquotes a variable to use its value instead of its name
  }, error = function(e) {
    message(sprintf("Error downloading %s: %s", ticker, e$message))
    return(NULL)
  })
}

stock_data <- lapply(parent_banks, get_bank_stockprices) %>% 
  purrr::reduce(full_join, by = 'Date') %>% 
  arrange(Date)

# Clean column names (remove exchange suffixes)
colnames(stock_data) <- gsub("\\..*$", "", colnames(stock_data))
colnames(stock_data) <- gsub("([^.\\-]+)[.\\-].*", "\\1", colnames(stock_data))
colnames(stock_data)[colnames(stock_data)=="HSBA"] <- "HSBC"

parent_banks <- gsub("([^.\\-]+)[.\\-].*", "\\1", parent_banks)
# STOXX index -------------------------------------------------------------
stoxx <- getSymbols("EXV1.DE", src = "yahoo", from = "2001-01-01", to = Sys.Date(), auto.assign = FALSE)
stoxx <- dailyReturn(Ad(stoxx), type = "log") #WHY?
colnames(stoxx) <- "stoxx" 

# VIX ---------------------------------------------------------------------
VIX <- getSymbols("^VIX", src = "yahoo", from = "2001-01-01", to = Sys.Date(), auto.assign = FALSE)
VIX <- Ad(VIX)
colnames(VIX) <- "VIX"


# HSBC --------------------------------------------------------------------
parent_banks['HSBC'] <- 'HSBC'

# Save to CSV, Rds
write.csv(stock_data, "Data/banks_stock_prices.csv", row.names = FALSE)
saveRDS(parent_banks, "Data/bank_names.Rds")
saveRDS(stoxx, "Data/stoxx.Rds")
saveRDS(VIX, "Data/VIX.Rds")
