#true parameter values
psi <- c(0, 0.1, 0.2)
names(psi) <- c("psi21", "psi32", "psi31")
expit <- function(x){
  1/(1+exp(-x))
}
#sample size
n <- 1e5
generate_data <- function(n, psi, Y3L2coef, linear = TRUE, censoring = FALSE, test = FALSE){
  #generate baseline time-varying and time-invariant confounders
  U <- rnorm(n)
  C <- rnorm(n)
  #t = 1
  L1 <- C+U+rnorm(n)
  Y1 <- C+U+rnorm(n)
  X1.logit <- 0.35*L1 + 0.35*Y1 + 0.1*C
  X1 <- rbinom(n,1,expit(X1.logit))
  
  #t = 2
  L2 <- -1*L1+X1+C+U+rnorm(n)
  Y2 <- psi["psi21"]*X1 + C+ U + rnorm(n)
  if(linear){
    Y2 <- Y2 + Y1 + L1
  } else{
    Y2 <- Y2 + 0.1*exp(Y1)+0.1*sqrt(abs(L1))
  }
  X2.logit <- -0.1*X1 + 0.2*L2 + 0.2*Y2 + 0.1*L1 + 0.1*Y1 + 0.1*C
  X2 <- rbinom(n,1,expit(X2.logit))
  
  #t = 3
  Y3 <- psi["psi32"]*X2 + psi["psi31"]*X1 + C + U + rnorm(n) 
  if(linear){
    if(test){
      Y3 <- Y3 + Y2+ Y1 + L1+Y3L2coef*L2
    }
    else{
      Y3 <- Y3 + Y2 + Y1 
    }
  } else{
    Y3 <- Y3 + pnorm(Y2)*exp(Y1)+Y3L2coef*L2 + L1
  }
  # Censoring at t=2
  if(censoring){
    C2.logit <- -2 + 0.5*X1 + 0.3*L1 + 0.2*Y1 +C + U
    C2 <- rbinom(n, 1, expit(C2.logit))
    
    # Apply censoring at t=2
    L2[C2 == 1] <- NA
    Y2[C2 == 1] <- NA
    X2[C2 == 1] <- NA
    
    # t = 3 (only for uncensored individuals)
    Y3 <- rep(NA, n)
    uncensored <- which(C2 == 0)
    Y3[uncensored] <- psi["psi32"]*X2[uncensored] + psi["psi31"]*X1[uncensored] + 
      C[uncensored] + U[uncensored] + rnorm(length(uncensored))
    
    if(linear){
      if(test){
        Y3[uncensored] <- Y3[uncensored] + Y2[uncensored] + Y1[uncensored] + 
          L1[uncensored] + Y3L2coef*L2[uncensored]
      } else {
        Y3[uncensored] <- Y3[uncensored] + Y2[uncensored] + Y1[uncensored]
      }
    } else {
      Y3[uncensored] <- Y3[uncensored] + pnorm(Y2[uncensored])*exp(Y1[uncensored]) + Y3L2coef*L2[uncensored] + L1[uncensored]
    }
    
    # Censoring at t=3 (only for those not censored at t=2)
    C3 <- rep(1, n)
    C3.logit <--2+ 0.5*X2 + 0.3*L2 + 0.2*Y2+ C + U
    C3[uncensored] <- rbinom(length(uncensored), 1, expit(C3.logit[uncensored]))
    
    # Apply censoring at t=3
    Y3[C3 == 1] <- NA
    return(data.frame(C,L1,Y1,X1,L2,Y2,X2,Y3,U,C2,C3))
  }
  else{
    return(data.frame(C,L1,Y1,X1,L2,Y2,X2,Y3,U))
  }
}
set.seed(9000)
generated_data <- generate_data(n, psi, 0, censoring = TRUE, test = TRUE,linear = FALSE)
sum(generated_data$C2)/n
sum(generated_data$C3)
(sum(generated_data$C3)-sum(generated_data$C2))/sum(1-generated_data$C2)
write.csv(generated_data, "C:/Users/marcu/OneDrive/Dokumenter/Master/masters-master/generated_data.csv", row.names = FALSE)


generate_data_test <- function(n, psi, Y3L2coef, linear = TRUE, censoring = FALSE, test = FALSE) {

  # Unmeasured confounder
  U <- rnorm(n)
  
  # Baseline confounder
  C_base <- rnorm(n)
  
  # Time t = 1
  L.1 <- C_base + U + rnorm(n)
  Y_base <- C_base + U + rnorm(n)
  A1.logit <- 0.35 * L.1 + 0.35 * Y_base + 0.1 * C_base
  A.1 <- rbinom(n, 1, expit(A1.logit))
  
  # Time t = 2
  L.2 <- -1 * L.1 + A.1 +Y_base+ C_base + U + rnorm(n)
  Y.2 <- psi["psi21"] * A.1 + C_base + U + rnorm(n)
  if (linear) {
    Y.2 <- Y.2 + Y_base + L.1
  } else {
    Y.2 <- Y.2 + 0.1 * exp(Y_base) + 0.1 * sqrt(abs(L.1))
  }
  A2.logit <- -0.1 * A.1 + 0.2 * L.2 + 0.2 * Y.2 + 0.1 * L.1 + 0.1 * Y_base + 0.1 * C_base
  A.2 <- rbinom(n, 1, expit(A2.logit))
  
  # Time t = 3
  Y.3 <- psi["psi32"] * A.2 + psi["psi31"] * A.1 + C_base + U + rnorm(n)
  if (linear) {
    if (test) {
      Y.3 <- Y.3 + Y.2 + Y_base + L.1 + Y3L2coef * L.2
    } else {
      Y.3 <- Y.3 + Y.2 + Y_base
    }
  } else {
    Y.3 <- Y.3 + pnorm(Y.2) * exp(Y_base) + Y3L2coef * L.2 + L.1
  }
  
  # Censoring mechanism
  C.2 <- C.3 <- rep(0, n)
  
  if (censoring) {
    # Censoring at t = 2
    C2.logit <- -2 + 0.5 * A.1 + 0.3 * L.1 + 0.2 * Y_base + C_base + U
    C.2 <- rbinom(n, 1, expit(C2.logit))
    
    L.2[C.2 == 1] <- NA
    Y.2[C.2 == 1] <- NA
    A.2[C.2 == 1] <- NA
    
    # Censoring at t = 3
    C3.logit <- -2 + 0.5 * A.2 + 0.3 * L.2 + 0.2 * Y.2 + C_base + U
    C.3 <- rbinom(n, 1, expit(C3.logit))
    C.3[C.2 == 1] <- 1
    
    Y.3[C.3 == 1] <- NA
  }
  
  # Combine all variables in wide format
  id <- 1:n
  data_wide <- as.data.frame(cbind(id, A.1,A.2,L.1,L.2,Y.2,Y.3,C.2,C.3,Y_base,C_base)
  )
  print(C.2)
  print(C.3)
  # Convert data to long format
  data_long <- reshape(
    data_wide,
    direction = "long",
    varying = c("A.1", "A.2", "L.1", "L.2", "Y.2", "Y.3", "C.2", "C.3"),
    v.names = c("A", "L", "Y", "C"),
    times = c(1, 2),
    timevar = "time",
    idvar = "id"
  )


  
  # Order data by ID and time
  data_long <- data_long[order(data_long$id, data_long$time), ]
  data_long$A[seq(1,2*n, by = 2)] <- as.numeric(A.1)
  data_long$A[seq(2,2*n, by = 2)] <- as.numeric(A.2)
  data_long$C[seq(1,2*n, by = 2)] <- as.numeric(C.2)
  data_long$C[seq(2,2*n, by = 2)] <- as.numeric(C.3)
  # Return the long-format data
  return(data_long)
}

n<- 100
data <- generate_data_test(n, psi, Y3L2coef = 0, censoring = TRUE)
data
library(gesttools)
data <- FormatData(data, idvar = "id",An = "A", varying = c("A", "L", "Y"),timevar = "time", Cn = "C",
                   GenerateHistory = TRUE, GenerateHistoryMax = 1)
data$time
outcomemodels <- list("Y~A+L +Y_base + C_base", "Y~A+L+Lag1Y+Lag1A + Lag1L +Y_base + C_base")
propensitymodel <- c("A~Lag1A + Y_base + Y+L + Lag1L + C_base")
censoringmodel <- c("C~Lag1A + Y_base + Y+L + Lag1L + C_base")
gest <- gestMultiple(data, idvar = "id", timevar = "time", Yn = "Y", An = "A", Cn = "C", 
                     outcomemodels = outcomemodels, propensitymodel = propensitymodel, 
                     censoringmodel = censoringmodel, type = 3)
gest$psi
