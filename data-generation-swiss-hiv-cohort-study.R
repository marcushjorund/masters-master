set.seed(420)
gamma <- c(-3,0.05,-1.5,0.1)
theta <- c(-0.405,0.0205,-0.00405)
n <- 1000
T <- 40
k <- 5

expit <- function(x){
  1/(1+exp(-x))
}
u_0 <- runif(n)
#matrix with T+1 rows, to account for initial conditions
u <- matrix(c(u_0,rep(0,n*T)), nrow = T+1, ncol = n)
epsilon_0 <- rnorm(n,mean=0, sd = sqrt(20))
epsilon <- matrix(data = c(epsilon_0,rep(0,n*T)), ncol = n, nrow = T+1,byrow = TRUE)
L_0 <- qgamma(u_0, shape = 3, scale = 154, lower.tail = FALSE)+epsilon[1,]

L <- matrix(c(L_0,rep(0,n*T)), nrow = T+1, ncol = n, byrow = TRUE)
A_start <- rep(0,n)
A_0 <- rbinom(n,size = 1,prob = expit(theta[1]+theta[3]*(L_0-500)))
A <- matrix(c(A_start,A_0, rep(NA,n*T)), nrow = T+2, ncol = n, byrow = TRUE)
lambda_0 <- expit(gamma[1]+gamma[3]*A_0)
lambda <- matrix(c(lambda_0,rep(0,n*T)), nrow = T+1, ncol = n, byrow = TRUE)
Y <- matrix(data = c(as.integer(lambda_0>u_0),rep(NA,n*T)), nrow = T+1, ncol = n, byrow = TRUE)
Y[,Y[1,]==1]<-1
T_star <- rep(NA,n)
T_star[A_0==1] <-0
delta <- matrix(data = rnorm(n*T, mean = 0 , sd = sqrt(0.05)), nrow = T, ncol = n, byrow = TRUE)
data_long <- matrix(data = NA, ncol = 6)
data_long
colnames(data_long) <- c("id", "U", "time", "Y", "A", "L")
counter <- 1
for (i in 1:n){
  t <- 1
  while((Y[t,i] == 0)&(t <(T+1))){
    u[t+1,i] <- min(1, max(0,u[t,i]+delta[t,i]))
    if (t%%k != 0){
      L[t+1,i] <- L[t,i]
      A[t+2,i] <- A[t+1,i]
      data_long[,c("A","L")]
    }
    else{
      epsilon[t+1,i] <- rnorm(1,mean = 100*(u[t+1,i]-2), sd=sqrt(50))
      L[t+1,i] <- max(0, L[t,i]+150*A[t+2-k,i]*(1-A[t-k+1,i])+epsilon[t+1,i])
      if(A[t+1,i] == 0){
        A[t+2,i] <- rbinom(1, size = 1, prob = expit(theta[1] + theta[2]*t+theta[3]*(L[t+1,i]-500)))
      }
      else{
        A[(t+2):(T+2),i] <- 1
      }
      if ((A[t+2,i]==1)&(A[t+2-k,i]==0)){
        T_star[i] <- t
      }
    }
    lambda[t+1,i] <- expit(gamma[1]+gamma[2]*((1-A[t+2,i])*t+A[t+2,i]*T_star[i])
                    +gamma[3]*A[t+2,i]+gamma[4]*A[t+2,i]*(t-T_star[i]))
    if((1-prod(1-lambda[rep(1,t),i]))>=u_0[i]){
      Y[(t+1):(T+1),i] <- 1
    }
    else{
      Y[t+1,i] <- 0
    }
  id <- i
  U <- u[t+1,i]
  time <- t
  y <- Y[t+1,i]
  a <- A[t+2,i]
  l <- L[t+1,i]
  data_long <- rbind(data_long,c(id,U,time,y,a,l))
  t <- t+1
  }
}
data_long <- data_long[2:dim(data_long)[1],]

L[,9]
Y[,9]
A[,9]
#Testing with g-estimation
library(gesttools)
data <- dataexamples(n = 1000, seed = 420, Censoring = TRUE)
data$datagestmult
?dataexamples
