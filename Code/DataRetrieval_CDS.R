### CDS Systemically Relevant Banks (cleaning up the excel from Bloomberg)

# Libraries -------------------------------------------------------------
packages <- c('here', 'readxl', 'dplyr')
instld <- sapply(packages, require, character.only = T)
if (sum(instld)<length(instld)) install.packages(packages[!instld])
invisible(lapply(packages, library, character.only = TRUE))
rm(packages, instld)


# Read the Excel ----------------------------------------------------------
cds_data <- read_excel("Data/CDS_data.xlsx", sheet = 1, col_names = FALSE)

### The odd columns are overlapping Date columns and the even ones - cds prices of european banks
headers <- cds_data[1, ]
cds_data <- cds_data[-1, ]
dates <- cds_data %>% 
  dplyr::select(seq(1, ncol(cds_data), by = 2)) %>%  # selecting odd columns
  mutate(across(where(is.character), ~ as.Date(as.numeric(.), origin = "1899-12-30")))

unique_dates <- dates%>% 
  pivot_longer(everything(), values_to = "Date") %>%  # Stack all dates
  distinct(Date) %>%  # Keep unique dates
  arrange(Date) %>%   # Sort chronologically
  mutate(Date = as.Date(Date))

# Get all PRICE columns (even-numbered)
prices <- cds_data %>% 
  dplyr::select(seq(2, ncol(cds_data), by = 2)) %>%   #selecting even columns
  mutate(across(everything(), function(x) as.numeric(x)))

names(prices) <- as.character(headers[seq(2, length(headers), by =2)])

merge_pr_date <- function(){
  helper <- function(i){
    data <- bind_cols(dates[, i], prices[, i])
    names(data)[1]<- "Date"
    data %>% 
      mutate(Date = as.Date(Date)) %>% # force simple Date instead of POSixt
      na.omit()
  }
  prices_list <- lapply(1:ncol(dates), helper)
  data <- left_join(unique_dates, prices_list[[1]], by = 'Date')
  for(i in 2:ncol(prices)){
    data <- left_join(data, prices_list[[i]], by = 'Date')
  }
  data
}

# Final data set with CDS prices 
prices <- merge_pr_date()
prices <- prices[-nrow(prices), ] #last row is empty

bank_names <- headers[seq(2, ncol(cds_data), by = 2)]

# Some checks and cleaning up the directory
cat("\nThe first date in sample is", format(prices[[1, 1]], "%Y-%m-%d"), "\nThe last date in sample is:", format(prices[[nrow(prices), 1]], "%Y-%m-%d"))

rm(list = setdiff(ls(), c("prices", "bank_names")))

# Save to CSV, Rds
write.csv(prices, file = "Data/banks_cds.csv", row.names = FALSE)
saveRDS(bank_names, "Data/bank_names_cds.Rds")

#names(bank_names)[which(names(bank_names) != bank_names_cds)]
#bank_names_cds[which(names(bank_names) != bank_names_cds)]
