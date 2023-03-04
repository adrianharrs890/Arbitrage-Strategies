rm(list = ls())
library(tidyverse)
library(rvest)
library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)
library(zoo)
library(quantmod)
library(tibble)
library(tseries)
require(gridExtra)

begin <-  "2022-01-01"
end <- as.character(as.Date(Sys.time()))


setwd('Desktop/stocks/')

spy_500 <- read_csv('constituents_csv.csv')

companies <- spy_500 %>%
  pull(Symbol)

dataList <- list()

# Step One
for(i in 1:length(companies)){
  
  if(i %in% c(43, 60,65,81,
              100,115,149,150,
              160,190,204,244,
              268,335,341,364, 
              474,496,499)){
    next
  }
  
  dataList[[i]] <- Ad(getSymbols(companies[i], from = begin, to = end,periodicity = "daily", auto.assign = FALSE)) %>%
    as.data.frame() %>%
    rownames_to_column(., "date") %>%
    rename(Price = 2) %>%
    mutate(Company = companies[i])
}

stocks <- do.call(rbind,dataList)
write_csv(stocks, 'stock_data.csv')

df <- stocks %>%
  rename(Symbol = Company) %>%
  left_join(
    spy_500
  )

left <- df %>%
  pull(Symbol) %>%
  unique()

# Step 2
vector <- list()
combs <- combn(left,2)
for(i in 1:117855){
  vector[[i]] <- combs[,i]
}

sample_list <- unique(lapply(vector , sort))

corData <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(corData ) <- c('pair', "corr")

# Step 3 
for(i in 1:length(sample_list)){
  one <- df %>%
    filter(Symbol %in% sample_list[[i]][1]) %>%
    pull(Price)
  two <- df %>%
    filter(Symbol %in% sample_list[[i]][2]) %>%
    pull(Price)
  
  if(length(one) !=  length(two)){
    next
  }
  model <- lm(log(one) ~ log(two) - 1)
  sprd<-residuals(model)
  corData[i, 'pair'] <- paste0(sample_list[[i]],collapse = "---")
  corData[i, 'corr'] <- cor(one, two)
  corData[i, 'beta'] <- as.numeric(coef(model)[1])
  corData[i, "pval"] <- adf.test(sprd, alternative="stationary", k=0)$p.value
}

write_csv(corData %>% 
  filter( pval<=0.05 & corr >0.95) %>%
  arrange(desc(corr)), 'pairs_data.csv')

data_prep <- function(inter){
  beta <- inter %>%
    pull(beta)
  
  pairs <- inter %>%
    separate(pair, c("one",'two')) %>%
    select(one,two) %>%
    gather(.,'var',"val") %>%
    pull(val) 
  
  one <- df %>%
    filter(Symbol ==  pairs[1]) %>%
    pull(Price)
  
  two <- df %>%
    filter(Symbol == pairs[2]) %>%
    pull(Price)
  
  tab <- df %>%
    filter(Symbol %in% pairs) %>%
    cbind(
      spread = one-beta*two
    ) %>%
    mutate(spread_z = (spread - mean(spread))/sd(spread),
           lower =  sd(spread)*-1,
           upper =   sd(spread),
           middle =  0,
           date = ymd(date),
           Price = as.numeric(Price),
           Sector = as.factor(Sector),
           x_y = case_when(
             pairs[1] == Symbol ~ "Y",
             TRUE ~ "X"
           )) %>%
    rowwise() %>%
    mutate(
      trade = case_when(
        spread_z  > upper ~ "Buy X Sell Y",
        spread_z  <  lower ~ "Buy Y Sell X",
        
        TRUE ~ "No Trade"
        
      )
    ) 
  
  return(tab)
}


# Step 4
inter <- corData %>% 
  filter( pval<=0.05 & corr >0.95) %>%
  arrange(pval) %>%
  slice(7)

dgp <- data_prep(inter)

library(tidyverse)
library(tidymodels)
library(modeltime)
library(timetk)
library(lubridate)

# Step 5
Y <- dgp %>%
  filter(x_y ==  "Y") %>%
  pull(Symbol) %>%
  unique()

X <- dgp %>%
  filter(x_y ==  "X") %>%
  pull(Symbol) %>%
  unique()

time_Y <- dgp %>%
  filter(Symbol ==  Y)  %>%
  mutate(date = ymd(date)) %>%
  select(date, Price) %>%
  set_names(c("date", "value")) 

time_X <- dgp %>%
  filter(Symbol ==  X)  %>%
  mutate(date = ymd(date)) %>%
  select(date, Price) %>%
  set_names(c("date", "value")) 

splits_Y <- time_Y  %>%
  time_series_split(assess = "3 months", cumulative = TRUE)

splits_X <- time_X  %>%
  time_series_split(assess = "3 months", cumulative = TRUE)

recipe_spec_Y <- recipe(value ~ date, training(splits_Y)) %>%
  step_timeseries_signature(date) %>%
  step_rm(contains("am.pm"), contains("hour"), contains("minute"),
          contains("second"), contains("xts")) %>%
  step_fourier(date, period = 365, K = 5) %>%
  step_dummy(all_nominal())

recipe_spec_X <- recipe(value ~ date, training(splits_X)) %>%
  step_timeseries_signature(date) %>%
  step_rm(contains("am.pm"), contains("hour"), contains("minute"),
          contains("second"), contains("xts")) %>%
  step_fourier(date, period = 365, K = 5) %>%
  step_dummy(all_nominal())


model_spec_prophet_boost_Y <- prophet_boost() %>%
  set_engine("prophet_xgboost", yearly.seasonality = TRUE) 

model_spec_prophet_boost_X <- prophet_boost() %>%
  set_engine("prophet_xgboost", yearly.seasonality = TRUE) 

workflow_fit_prophet_boost_Y <- workflow() %>%
  add_model(model_spec_prophet_boost_Y) %>%
  add_recipe(recipe_spec) %>%
  fit(training(splits_Y))

workflow_fit_prophet_boost_X <- workflow() %>%
  add_model(model_spec_prophet_boost_X) %>%
  add_recipe(recipe_spec) %>%
  fit(training(splits_X))


model_table_Y <- modeltime_table(
  workflow_fit_prophet_boost_Y
) 

model_table_X <- modeltime_table(
  workflow_fit_prophet_boost_X
) 

calibration_table_Y <- model_table %>%
  modeltime_calibrate(testing(splits_Y))

calibration_table_X <- model_table %>%
  modeltime_calibrate(testing(splits_X))

pred_data_Y <- calibration_table_Y  %>%
  modeltime_refit(time_Y) %>%
  modeltime_forecast(
    h = "2 months",
    actual_data = time_Y
  ) %>%
  as.data.frame() %>%
  filter(.index > max(time$date)) %>%
  rename(date =  .index, 
         Price = .value) %>%
  mutate(Symbol = Y ) %>%
  select(date,Price,Symbol)

pred_data_X <- calibration_table_X  %>%
  modeltime_refit(time_X) %>%
  modeltime_forecast(
    h = "2 months",
    actual_data = time_X
  ) %>%
  as.data.frame() %>%
  filter(.index > max(time$date)) %>%
  rename(date =  .index, 
         Price = .value) %>%
  mutate(Symbol = X ) %>%
  select(date,Price,Symbol)


forcasted_data <- dgp %>%
  bind_rows(
    pred_data_Y 
  ) %>%
  bind_rows(
    pred_data_X
  )

# Step 6

y <- forcasted_data %>%
  filter(Symbol == Y) %>%
  pull(Price)

x <- forcasted_data %>%
  filter(Symbol == X) %>%
  pull(Price)

a0_y <- y[1]
a0_x <- x[1]

P0_y  <- matrix(1)
P0_x  <- matrix(1)

dt_y  <- ct_y  <- matrix(0)
dt_x  <- ct_x  <- matrix(0)

Zt_y  <- Tt_y  <- matrix(1)
Zt_x  <- Tt_x  <- matrix(1)

theta_y  <- theta_var_y <- rep(NA, length(y) + 1)
theta_x  <- theta_var_x <- rep(NA, length(x) + 1)

theta_y[1] <- a0_y 
theta_x[1] <- a0_x

theta_var_y[1] <- P0_y 
theta_var_x[1] <- P0_x

# Random weights
w <- sample(seq(0.00, .4, length = 10000),6)

sigma_w_y <- sigma_w_x <- sqrt(w[1])
sigma_v_y <- sigma_v_x <- sqrt(w[2])

G_t_y <- Tt_y
G_t_x <- Tt_x

F_t_y <- Zt_y
F_t_x <- Zt_x


for (i in 1:length(y)) {
  
  theta_hat_y <- theta_y[i]
  theta_hat_x <- theta_x[i]
  
  e_t_y <- y[i] - theta_hat_y * G_t_y* F_t_y
  e_t_x <- x[i] - theta_hat_x * G_t_x* F_t_x
  
  R_t_y <- G_t_y * theta_var_y[i] * G_t_y + sigma_w_y ^ 2
  R_t_x <- G_t_x * theta_var_x[i] * G_t_x + sigma_w_x ^ 2

  theta_y[i + 1] <- G_t_y * theta_hat_y + R_t_y * F_t_y * (sigma_v_y ^ 2 + F_t_y ^ 2 * R_t_y) ^ (-1) * e_t_y
  theta_x[i + 1] <- G_t_x * theta_hat_x + R_t_x * F_t_x * (sigma_v_x ^ 2 + F_t_x ^ 2 * R_t_x) ^ (-1) * e_t_x
  
  theta_var_y[i + 1] <- R_t_y - R_t_y * F_t_y * (sigma_v_y ^ 2 + F_t_y ^ 2 * R_t_y) ^ (-1) * F_t_y * R_t_y
  theta_var_x[i + 1] <- R_t_x - R_t_x * F_t_x * (sigma_v_x ^ 2 + F_t_x ^ 2 * R_t_x) ^ (-1) * F_t_x * R_t_x
}


theta_y <- theta_y[-1]
theta_x <- theta_x[-1]

theta_var_y <- theta_var_y[-1]
theta_var_x <- theta_var_x[-1]


data <- tibble(y_algebra = theta_y,
               y = y,
               id = 1:length(y),
               x = x, 
               x_algebra = theta_x)

data %>%
  gather(.,'var',"val",- id) %>%
  ggplot(.) +
  aes(id, val, group = var, color = val) +
  geom_line() 


forcasted_data






# Old code 


# 'DISH','ILMN'
# Pairs
i = 5
for(i in 5:74){
  #inter <- corData %>% 
  # filter( pval<=0.05 & corr >0.95) %>%
   # arrange(desc(corr)) %>%
   # slice(i)
  
  inter <- corData %>% 
    filter( pval<=0.05 & corr >0.95) %>%
    arrange(pval) %>%
    slice(1)
  
  
  beta <- inter %>%
    pull(beta)
  
  pairs <- inter %>%
    separate(pair, c("one",'two')) %>%
    select(one,two) %>%
    gather(.,'var',"val") %>%
    pull(val) 
  
  one <- df %>%
    filter(Symbol ==  pairs[1]) %>%
    pull(Price)
  
  two <- df %>%
    filter(Symbol == pairs[2]) %>%
    pull(Price)
  
  plot_data <- df %>%
    filter(Symbol %in% pairs) %>%
    cbind(
      spread = one-beta*two
    ) %>%
    mutate(lower =  quantile(spread,.25),
           upper =   quantile(spread,.75),
           middle =  quantile(spread,.5),
           date = ymd(date),
           Price = as.numeric(Price),
           Sector = as.factor(Sector),
           x_y = case_when(
             pairs[1] == Symbol ~ "Y",
             TRUE ~ "X"
           )) %>%
    rowwise() %>%
    mutate(
      trade = case_when(
        spread > upper ~ "Buy X Sell Y",
        spread <  lower ~ "Buy Y Sell X",
        
        TRUE ~ "No Trade"
        
      )
    )
  
  
  
  p1 <- plot_data %>% 
    filter(date >= as.Date(Sys.time()) - 60) %>%
    ggplot(.) +
    aes(date,spread, color = trade, group = trade) +
    geom_point(size = 1) +
    geom_hline(aes(yintercept=upper), color = "red", linetype = "dotted") +
    geom_hline(aes(yintercept=middle), color = "blue", linetype = "dotted")+
    geom_hline(aes(yintercept=lower), color = "green", linetype = "dotted") +
    ggtitle(paste("Y:",pairs[1], "X:", pairs[2])) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  p2 <- plot_data %>% 
    filter(date >= as.Date(Sys.time()) - 60) %>%
    ggplot(.) +
    aes(date,Price, color =  Symbol, group = Symbol) +
    geom_line() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(pairs[1],"---",pairs[2])) 
  
  grid.arrange(p1, p2, ncol=2)
  
}
place <- c()
for(i in 5:74){
  #inter <- corData %>% 
   # filter( pval<=0.05 & corr >0.95) %>%
   # arrange(desc(corr)) %>%
   # slice(i)
  
  
  inter <- corData %>% 
    filter( pval<=0.05 & corr >0.95) %>%
    arrange(pval) %>%
    slice(1)
  
  beta <- inter %>%
    pull(beta)
  
  pairs <- inter %>%
    separate(pair, c("one",'two')) %>%
    select(one,two) %>%
    gather(.,'var',"val") %>%
    pull(val) 
  
  one <- df %>%
    filter(Symbol ==  pairs[1]) %>%
    pull(Price)
  
  two <- df %>%
    filter(Symbol == pairs[2]) %>%
    pull(Price)
  
  tab <- df %>%
    filter(Symbol %in% pairs) %>%
    cbind(
      spread = one-beta*two
    ) %>%
    mutate(spread_z = (spread - mean(spread))/sd(spread),
           lower =  sd(spread)*-1,
           upper =   sd(spread),
           middle =  0,
           date = ymd(date),
           Price = as.numeric(Price),
           Sector = as.factor(Sector),
           x_y = case_when(
             pairs[1] == Symbol ~ "Y",
             TRUE ~ "X"
           )) %>%
    rowwise() %>%
    mutate(
      trade = case_when(
        spread_z  > upper ~ "Buy X Sell Y",
        spread_z  <  lower ~ "Buy Y Sell X",
        
        TRUE ~ "No Trade"
        
      )
    ) %>%
    ungroup() %>%
    summarise(mean_trade = mean(trade %in% c("Buy X Sell Y","Buy Y Sell X" )))
  
  if(as.numeric(tab) != 0){
    place[i] <- i
  }
  
  

  
}


good_slice <- as.vector(na.omit(place))

# 1

inter <- corData %>% 
  filter( pval<=0.05 & corr >0.95) %>%
  arrange(desc(corr)) %>%
  slice(good_slice[12])


inter <- corData %>% 
  filter( pval<=0.05 & corr >0.95) %>%
  arrange(pval) %>%
  slice(1)


beta <- inter %>%
  pull(beta)

pairs <- inter %>%
  separate(pair, c("one",'two')) %>%
  select(one,two) %>%
  gather(.,'var',"val") %>%
  pull(val) 

one <- df %>%
  filter(Symbol ==  pairs[1]) %>%
  pull(Price)

two <- df %>%
  filter(Symbol == pairs[2]) %>%
  pull(Price)

df %>%
  filter(Symbol %in% pairs) %>%
  cbind(
    spread = one-beta*two
  ) %>%
  mutate(spread_z = (spread - mean(spread))/sd(spread),
         lower =  sd(spread)*-1,
         upper =   sd(spread),
         middle =  0,
         date = ymd(date),
         Price = as.numeric(Price),
         Sector = as.factor(Sector),
         x_y = case_when(
           pairs[1] == Symbol ~ "Y",
           TRUE ~ "X"
         )) %>%
  rowwise() %>%
  mutate(
    trade = case_when(
      spread_z  > upper ~ "Buy X Sell Y",
      spread_z  <  lower ~ "Buy Y Sell X",
      
      TRUE ~ "No Trade"
      
    )
  )  %>% 
  ggplot(.) +
  aes(date,spread_z, color = trade, group = trade) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=upper), color = "red", linetype = "dotted") +
  geom_hline(aes(yintercept=middle), color = "blue", linetype = "dotted")+
  geom_hline(aes(yintercept=lower), color = "green", linetype = "dotted") +
  ggtitle(paste("Y:",pairs[1], "X:", pairs[2])) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


df %>%
  filter(Symbol %in% pairs) %>%
  cbind(
    spread = one-beta*two
  ) %>%
  mutate(spread_z = (spread - mean(spread))/sd(spread),
         lower =  sd(spread)*-1,
         upper =   sd(spread),
         middle =  0,
         date = ymd(date),
         Price = as.numeric(Price),
         Sector = as.factor(Sector),
         x_y = case_when(
           pairs[1] == Symbol ~ "Y",
           TRUE ~ "X"
         )) %>%
  rowwise() %>%
  mutate(
    trade = case_when(
      spread_z  > upper ~ "Buy X Sell Y",
      spread_z  <  lower ~ "Buy Y Sell X",
      
      TRUE ~ "No Trade"
      
    )
  )   %>% 
  mutate(plot_symbol = paste(Symbol,x_y )) %>%
  ggplot(.) +
  aes(date,Price, color =  plot_symbol, group = plot_symbol) +
  geom_line() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste(pairs[1],"---",pairs[2])) 








plot_data <- df %>%
  filter(Symbol %in% pairs) %>%
  cbind(
    spread = one-beta*two
  ) %>%
  mutate(lower =  quantile(spread,.25),
         upper =   quantile(spread,.75),
         middle =  quantile(spread,.5),
         date = ymd(date),
         Price = as.numeric(Price),
         Sector = as.factor(Sector),
         x_y = case_when(
           pairs[1] == Symbol ~ "Y",
           TRUE ~ "X"
         )) %>%
  rowwise() %>%
  mutate(
    trade = case_when(
      spread > upper ~ "Buy X Sell Y",
      spread <  lower ~ "Buy Y Sell X",
      
      TRUE ~ "No Trade"
      
    )
  )


plot_data %>% 
  mutate(plot_symbol = paste(Symbol,x_y )) %>%
  ggplot(.) +
  aes(date,Price, color =  plot_symbol, group = plot_symbol) +
  geom_line() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste(pairs[1],"---",pairs[2])) 

df %>%
  filter(Symbol %in% pairs) %>%
  cbind(
    spread = one-beta*two
  ) %>%
  mutate(spread_z = (spread - mean(spread))/sd(spread),
         lower =  sd(spread)*-1,
         upper =   sd(spread),
         middle =  0,
         date = ymd(date),
         Price = as.numeric(Price),
         Sector = as.factor(Sector),
         x_y = case_when(
           pairs[1] == Symbol ~ "Y",
           TRUE ~ "X"
         )) %>%
  rowwise() %>%
  mutate(
    trade = case_when(
      spread_z  > upper ~ "Buy X Sell Y",
      spread_z  <  lower ~ "Buy Y Sell X",
      
      TRUE ~ "No Trade"
      
    )
  )  %>% 
  ggplot(.) +
  aes(date,spread_z, color = trade, group = trade) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=upper), color = "red", linetype = "dotted") +
  geom_hline(aes(yintercept=middle), color = "blue", linetype = "dotted")+
  geom_hline(aes(yintercept=lower), color = "green", linetype = "dotted") +
  ggtitle(paste("Y:",pairs[1], "X:", pairs[2])) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


p1 <- plot_data %>% 
  ggplot(.) +
  aes(date,spread, color = trade, group = trade) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=upper), color = "red", linetype = "dotted") +
  geom_hline(aes(yintercept=middle), color = "blue", linetype = "dotted")+
  geom_hline(aes(yintercept=lower), color = "green", linetype = "dotted") +
  ggtitle(paste("Y:",pairs[1], "X:", pairs[2])) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


p2 <- plot_data %>% 
  mutate(plot_symbol = paste(Symbol,x_y )) %>%
  ggplot(.) +
  aes(date,Price, color =  plot_symbol, group = plot_symbol) +
  geom_line() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste(pairs[1],"---",pairs[2])) 

grid.arrange(p1, p2, ncol=2)
