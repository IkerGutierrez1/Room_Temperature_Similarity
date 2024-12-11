library(dplyr)
library(zoo)
library(ggplot2)
library(tidyr)
library(lubridate)

load_data <- function(path){
  df <- read.csv2(path)
  
  df <- df %>%
    filter(!is.na(timestamp))
  
  df <- interpolate_NA(df)
  
  return (df)
}


interpolate_NA <- function(df){
  df <- df %>%
    mutate(across(-timestamp, ~ as.numeric(as.character(.))))
  col_NA <- colSums(is.na(df))
  print(col_NA)
  for(col in names(col_NA)){
    if(col_NA[col] > 0){
      print(paste("Columna:", col, "- Interpolando", col_NA[col], "valores NA"))
      df <- df %>%
        mutate(!!col := na.approx(.[[col]], na.rm = FALSE))
    }
  }

  df <- df %>%
    filter(!rowSums(is.na(select(., -c(timestamp)))) == (ncol(df) - 2))
  return (df)
}
