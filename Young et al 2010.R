sim <- 1000
n <- 2500
K <- 9
#Simulate the counterfactual survival time under no treatment T_0
T_0 <- rexp(n,0.01)
V_start <- A_start<- Y_0 <- 0
#id,U,time,y,a,l needed for gesttools
snaftm_data <- data.frame()