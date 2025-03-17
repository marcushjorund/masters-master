set.seed(123)
n <- 1000
n_days <- 30

# Baseline variables
hospital <- sample(1:4, n, replace=TRUE)
age <- round(rnorm(n, 60, 10), 1)
gender <- rbinom(n, 1, 0.5)
prev_visits <- rpois(n, 3)
U <- rnorm(n)

# Time-varying discharge rates (instrument Z)
hospital_base_rates <- c(0.25, 0.35, 0.30, 0.40)
daily_discharge_rates <- matrix(nrow=4, ncol=n_days)
for(h in 1:4) {
  daily_discharge_rates[h,] <- plogis(
    qlogis(hospital_base_rates[h]) + rnorm(n_days, sd=0.3)
  )
}

# Simulate length of stay using baseline confounders and U
simulate_los <- function(n) {
  los <- numeric(n)
  
  for (i in 1:n) {
    h <- hospital[i]
    day <- 1
    base_stay <- 1
    additional_stay <- 0
    
    while (day <= n_days) {
      z <- daily_discharge_rates[h, day]
      
      # Logistic regression model for daily discharge probability
      logit_discharge <- -2.5 + 
        3*z - 
        0.7*U[i] + 
        0.03*age[i] - 
        0.4*gender[i] + 
        0.1*prev_visits[i]
      
      p_discharge <- plogis(logit_discharge)
      
      if (day > base_stay && rbinom(1, 1, p_discharge)) break
      
      day <- day + 1
      additional_stay <- additional_stay + 1
    }
    
    los[i] <- min(base_stay + additional_stay, n_days)
  }
  
  return(los)
}

# Mortality simulation (further adjusted for ~5% cumulative mortality)
simulate_mortality <- function(los) {
  death_day <- rep(NA, n)
  
  for(i in 1:n) {
    for(day in 1:n_days) {
      if(is.na(death_day[i])) {
        # In-hospital mortality (lower baseline hazard)
        if(day <= los[i]) {
          hazard <- exp(-8.5 + 
                          0.015*age[i] + 
                          0.12*gender[i] + 
                          0.04*prev_visits[i] + 
                          0.18*U[i] - 
                          0.1)
        } else { # Post-discharge mortality (risk decreases over time)
          days_post <- day - los[i]
          hazard <- exp(-6.5 + 
                          0.015*age[i] + 
                          0.12*gender[i] + 
                          0.04*prev_visits[i] + 
                          0.18*U[i] - 
                          0.08*days_post)
        }
        
        # Simulate mortality event
        if(rbinom(1, 1, min(1, hazard))) death_day[i] <- day
      }
    }
  }
  
  return(death_day)
}

# Longitudinal dataset creation with death truncation
create_long_data <- function(los, death_day) {
  long_data <- vector("list", n)
  
  for(i in 1:n) {
    h <- hospital[i]
    death_time <- ifelse(is.na(death_day[i]), n_days, death_day[i])
    days <- seq_len(death_time)
    A <- as.numeric(days <= los[i])
    Z <- daily_discharge_rates[h, days]
    event <- c(rep(0, death_time-1), as.numeric(!is.na(death_day[i])))
    
    long_data[[i]] <- data.frame(
      patient_id = i,
      hospital = factor(h),
      day = days,
      A = A,
      Z = Z,
      age = age[i],
      gender = factor(gender[i], levels = c(0,1), labels = c("Male", "Female")),
      prev_visits = prev_visits[i],
      event = event,
      days_post_discharge = pmax(days - los[i], 0)
    )
  }
  
  return(do.call(rbind, long_data))
}

# Execute simulation
los <- simulate_los(n)
death_day <- simulate_mortality(los)
sim_data <- create_long_data(los, death_day)

# Verify mortality rate
mortality_rate <- mean(tapply(sim_data$event, sim_data$patient_id, max))
print(paste("Cumulative mortality rate:", round(mortality_rate, 4)))
sim_data
