library(tidyverse)
library(xlsx)

###### fit weibull to median and IQR#####
# Function to estimate Weibull parameters
estimate_weibull_params <- function(data) {
  # Extract percentiles
  t_25 <- data$p25
  t_50 <- data$p50
  t_75 <- data$p75
  
  # Log of survival times
  log_t_25 <- log(t_25)
  log_t_50 <- log(t_50)
  log_t_75 <- log(t_75)
  
  # Log-log survival values based on percentiles
  loglog_S_25 <- log(-log(1 - 0.25))
  loglog_S_50 <- log(-log(1 - 0.50))
  loglog_S_75 <- log(-log(1 - 0.75))
  
  # Estimate beta using the differences between the log-log survival values and log times
  beta_1 <- (loglog_S_25 - loglog_S_50) / (log_t_25 - log_t_50)
  beta_2 <- (loglog_S_50 - loglog_S_75) / (log_t_50 - log_t_75)
  
  # Estimate beta and theta
  # Average the two beta estimates
  beta <- (beta_1 + beta_2) / 2
  
  # Estimate theta using the equation for the median (t_50)
  log_theta <- log_t_50 - (loglog_S_50 / beta)
  theta <- exp(log_theta)
  
  list(beta = beta, theta = theta)
}

# # Example data for the two transitions
# data <- data.frame(
#   Transition = c("HHD_Dth", "HHD_ICHD"),
#   p25 = c(396, 69),
#   p50 = c(883, 220),
#   p75 = c(1786, 597)
# )
# 
# # Apply the function to each group of transitions
# params_tibble <- data %>%
#   group_by(Transition) %>%
#   nest() %>%
#   mutate(weibull_params = map(data, ~ estimate_weibull_params(.x))) %>%
#   # Extract beta and theta, and remove weibull_params using transmute
#   transmute(
#     Transition,
#     p25 = map_dbl(data, ~ .x$p25),
#     p50 = map_dbl(data, ~ .x$p50),
#     p75 = map_dbl(data, ~ .x$p75),
#     beta = map_dbl(weibull_params, "beta"),
#     theta = map_dbl(weibull_params, "theta")
#   )
# 
# # Display the result
# params_tibble



TransitionTableData <- read.xlsx("C:/Users/Stacey.Croft/OneDrive - Midlands and Lancashire CSU/Desktop/1210 Renal Services previously 1149 and 1150/DES KRT/TransitionsDataset.xlsx",
                                 sheetName="Sheet1")
#View(TransitionTableData)

TransitionTableData<-TransitionTableData%>% 
  mutate(DaysDuration=as.numeric(DaysDuration))

z1 <- TransitionTableData %>%
  filter(DaysDuration != "x" & DaysDuration != 0) %>%
  dplyr::select(From, To, DaysDuration, Percent) %>%
  mutate(Transition = paste0(From, "_", To)) %>%
  group_by(Transition) %>%
  pivot_wider(names_from = Percent, values_from = DaysDuration) %>%
  rename(p25 = "25", p50 = "50", p75 = "75")
# Apply the function to each group of transitions
z1_params <- z1 %>%
  group_by(Transition) %>%
  nest() %>%
  mutate(weibull_params = map(data, ~ estimate_weibull_params(.x))) %>%
  mutate()%>%
  transmute(
    Transition,
    p25 = map_dbl(data, ~ .x$p25),
    p50 = map_dbl(data, ~ .x$p50),
    p75 = map_dbl(data, ~ .x$p75),
    shape = map_dbl(weibull_params, "beta"),
    scale = map_dbl(weibull_params, "theta"),
    mean =  map_dbl(weibull_params, "theta") * gamma(1 + 1 / map_dbl(weibull_params, "beta")),
    calc_p25 = map_dbl(weibull_params, "theta") * (log(4 / 3) ^ (1 / map_dbl(
      weibull_params, "beta"
    ))),
    median = map_dbl(weibull_params, "theta") * (log(2) ^ (1 / map_dbl(
      weibull_params, "beta"
    ))),
    calc_p75 = map_dbl(weibull_params, "theta") * (log(4) ^ (1 / map_dbl(
      weibull_params, "beta"
    )))
  )
z1_params <- data.frame(z1_params)

# Display the result
print(z1_params)

# Write to Excel file
write.xlsx(z1_params, file = "Weibull_params.xlsx", sheetName = "Sheet1", row.names = FALSE)

######transition probabilities#####
#compute the hazard rate for a Weibull distribution at time t using shape and scale, in a hazard function
# Function for weibull hazard function
hazard_weibull <- function(t, shape, scale) {
  (shape / scale) * (t / scale)^(shape - 1)
}

# Transition probabilities for remaining states in DES (excluding one or more states)
transition_probabilities <- function(t, shape_values, scale_values, exclude_states = NULL) {
  # Exclude specific states (if any)
  if (!is.null(exclude_states)) {
    shape_values <- shape_values[-exclude_states]
    scale_values <- scale_values[-exclude_states]
  }
  
  # Number of remaining states
  N <- length(shape_values)
  
  # Calculate the hazard rates for each remaining state at time t
  hazard_rates <- sapply(1:N, function(i) {
    hazard_weibull(t, shape_values[i], scale_values[i])
  })
  
  # Total hazard rate (sum of hazard rates)
  total_hazard_rate <- sum(hazard_rates)
  
  # Transition probabilities for remaining states (hazard rate / total hazard rate)
  probabilities <- hazard_rates / total_hazard_rate
  
  return(probabilities)  # These will sum to 1
}

# Example using ICHD
# Shape and scale parameters for 5 state transitions: 
# ICHD_ICHD
# ICHD_Dth
# ICHD_HHD
# ICHD_PD
# ICHD_Tx
shape_values <- c(1.11, 0.96, 0.63, 0.52, 1.07)# shape parameters
scale_values <- c(1494.1, 1591.5, 670.3, 463.0, 1702.8)# scale parameters

# Calculate transition probabilities at time t = 1825 (end of chart in draft paper)
t <- c(0,365, 730, 1095, 1460, 1825)

# Exclude the 2nd state (index = 2) from the transition probabilities
#exclude_states <- c(1,5)
exclude_states <- NULL

# Get transition probabilities for remaining states
probs <- transition_probabilities(t, shape_values, scale_values, exclude_states)
sum(probs)#check sum to 1
print(probs)

# # Choose next state based on the probabilities
# next_state <- sample(1:length(probs), size = 1, prob = probs)
# print(paste("Next state is:", next_state))
