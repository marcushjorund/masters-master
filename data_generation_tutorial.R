options(scipen=1, digits = 2)
#true parameter values
psi <- c(0.0, 0.1, 0.2)
names(psi) <- c("psi21", "psi32", "psi31")
expit <- function(x){
  1/(1+exp(-x))
}
#sample size
n <- 1e5
generate_data <- function(n, linear = TRUE){
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
    Y2 <- Y2 + 0.1*exp(Y1)
  }
  X2.logit <- -0.1*X1 + 0.2*L2 + 0.2*Y2 + 0.1*L1 + 0.1*Y1 + 0.1*C
  X2 <- rbinom(n,1,expit(X2.logit))
  
  #t = 3
  Y3 <- psi["psi32"]*X2 + psi["psi31"]*X1 + C + U + rnorm(n) 
  if(linear){
    Y3 <- Y3 + Y2+Y1
  } else{
    Y3 <- Y3 + pnorm(Y2)*exp(Y1)
  }
  return(data.frame(C,L1,Y1,X1,L2,Y2,X2,Y3))
}
set.seed(42)
generated_data <- generate_data(n)
head(generated_data)
write.csv(generated_data, "C:/Users/marcu/OneDrive/Dokumenter/Master/masters-master/generated_data.csv", row.names = FALSE)


