source("data_generation_tutorial.R")
#Linear outcome model, not propensity score adjusted
library(ipw)
library(matrixStats)

psi <- c(0,0.1,0.2)
names(psi) <- c("psi21","psi32","psi31")
Y3L2coef <-0
B <- 200
n <- 1e4
names(psi) <- c("psi21", "psi32", "psi31")

#Bootstrapping linear outcome, no propensity score

#3 columns for the three effect estimates, and one column for the presumed biased estimate for psi31
psi_matrix <- matrix(0, nrow = B, ncol = 4)
for (b in 1:B){
  generated_data <- generate_data(n, psi, Y3L2coef, linear = TRUE, test = TRUE)
  #effect of X2 on Y3
  model_y3x2 <- lm(Y3 ~ X2+Y2+X1+L1+L2+Y1+L2+C, data = generated_data)
  psi32 <- coef(model_y3x2)["X2"]
  #effect of X1 on Y2
  model_y2x1 <- lm(Y2~X1+L1+Y1+C, data = generated_data)
  psi21 <- coef(model_y2x1)["X1"]
  #residualise 
  generated_data[,"R32"] <- generated_data[,"Y3"]-psi32*generated_data[,"X2"]
  model_y3x1 <- lm(R32 ~ L1+X1+Y1+C, na.action = "na.fail",data = generated_data)
  psi31 <- coef(model_y3x1)["X1"]
  psi_est <- c(psi21, psi32, psi31)
  names(psi_est) <- c("psi21", "psi32", "psi31")
  psi_matrix[b,] <- c(psi_est,coef(model_y3x2)["X1"])
}
estimates <- colMeans(psi_matrix)
biases <- psi-estimates[1:3]
SD <- sqrt(diag(cov(psi_matrix)))
CI <- colQuantiles(psi_matrix, probs = c(0.025, 0.975))
estimates
biases
CI
SD

#Bootstrapping propensity score adjusted linear model

#3 columns for the three effect estimates, and one column for the presumed biased estimate for psi31
psi_matrix_propscore_adj <- matrix(0, nrow = B, ncol = 4)
for (b in 1:B){
  generated_data <- generate_data(n, psi, Y3L2coef, linear = TRUE, test = TRUE)
  p2_model <- glm(X2~X1+Y1+Y2+L1+L2+C, family = binomial("logit"), data = generated_data)
  p1_model <- glm(X1~Y1 + L1 + C, family = binomial("logit"), data = generated_data)
  generated_data[,"P2"] <- predict.glm(p2_model, type = "response")
  generated_data[,"P1"] <- predict.glm(p1_model, type = "response")
  #effect of X2 on Y3
  model_y3x2 <- lm(Y3 ~ X2+X1+Y2+Y1+L1+L2+C+P2, data = generated_data)
  psi32 <- coef(model_y3x2)["X2"]
  
  #residualise 
  generated_data[,"R32"] <- generated_data[,"Y3"]-psi32*generated_data[,"X2"]
  #effect of X1 on Y3
  model_y3x1 <- lm(R32 ~ X1+Y1+L1+C+P1, data = generated_data)
  psi31 <- coef(model_y3x1)["X1"]
  #effect of X1 on Y2
  
  model_y2x1 <- lm(Y2~X1+L1+Y1+C+P1, data = generated_data)
  psi21 <- coef(model_y2x1)["X1"]
  
  psi_est <- c(psi21, psi32, psi31)
  names(psi_est) <- c("psi21", "psi32", "psi31")
  psi_matrix_propscore_adj[b,] <- c(psi_est, coef(model_y3x2)["X1"])
}

estimates_propscore_adj <- colMeans(psi_matrix_propscore_adj)
biases_propscore_adj <- psi-estimates_propscore_adj[1:3]
SD_propscore_adj <- sqrt(diag(cov(psi_matrix_propscore_adj)))
CI_propscore_adj <- colQuantiles(psi_matrix_propscore_adj, probs = c(0.025, 0.975))
#Note that propensity score adjusted increases estimated variances of estimates

#Non linear model, not propensity score adjusted

#3 columns for the three effect estimates, and one column for the presumed biased estimate for psi31
psi_matrix_nonlinear <- matrix(0, nrow = B, ncol = 4)
for (b in 1:B){
  generated_data <- generate_data(n, psi, Y3L2coef, linear = FALSE, test = TRUE)
  #effect of X2 on Y3
  model_y3x2 <- lm(Y3 ~ X2+X1+Y2+Y1+L1+L2+C, data = generated_data)
  psi32 <- coef(model_y3x2)["X2"]
  
  #residualise 
  generated_data[,"R32"] <- generated_data[,"Y3"]-psi32*generated_data[,"X2"]
  #effect of X1 on Y3
  model_y3x1 <- lm(R32 ~ X1+Y1+L1+C, data = generated_data)
  psi31 <- coef(model_y3x1)["X1"]
  #effect of X1 on Y2
  
  model_y2x1 <- lm(Y2~X1+L1+Y1+C, data = generated_data)
  psi21 <- coef(model_y2x1)["X1"]
  
  psi_est <- c(psi21, psi32, psi31)
  names(psi_est) <- c("psi21", "psi32", "psi31")
  psi_matrix_nonlinear[b,] <- c(psi_est, coef(model_y3x2)["X1"])
}

estimates_nonlinear <- colMeans(psi_matrix_nonlinear)
biases_nonlinear <- psi-estimates_nonlinear[1:3]
SD_nonlinear <- sqrt(diag(cov(psi_matrix_nonlinear)))
CI_nonlinear <- colQuantiles(psi_matrix_nonlinear, probs = c(0.025, 0.975))

#Propensity score adjusted non-linear model

#3 columns for the three effect estimates, and one column for the presumed biased estimate for psi31
psi_matrix_nonlinear_propscore_adj <- matrix(0, nrow = B, ncol = 4)
for (b in 1:B){
  generated_data <- generate_data(n, psi, Y3L2coef, linear = FALSE, test = TRUE)
  p2_model <- glm(X2~X1+Y1+Y2+L1+L2+C, family = binomial("logit"), data = generated_data)
  p1_model <- glm(X1~Y1 + L1 + C, family = binomial("logit"), data = generated_data)
  generated_data[,"P2"] <- predict.glm(p2_model, type = "response")
  generated_data[,"P1"] <- predict.glm(p1_model, type = "response")
  #effect of X2 on Y3
  model_y3x2 <- lm(Y3 ~ X2+X1+Y2+Y1+L1+L2+C+P2, data = generated_data)
  psi32 <- coef(model_y3x2)["X2"]
  
  #residualise 
  generated_data[,"R32"] <- generated_data[,"Y3"]-psi32*generated_data[,"X2"]
  #effect of X1 on Y3
  model_y3x1 <- lm(R32 ~ X1+Y1+L1+C+P1, data = generated_data)
  psi31 <- coef(model_y3x1)["X1"]
  #effect of X1 on Y2
  
  model_y2x1 <- lm(Y2~X1+L1+Y1+C+P1, data = generated_data)
  psi21 <- coef(model_y2x1)["X1"]
  
  psi_est <- c(psi21, psi32, psi31)
  names(psi_est) <- c("psi21", "psi32", "psi31")
  psi_matrix_nonlinear_propscore_adj[b,] <- c(psi_est, coef(model_y3x2)["X1"])
}
estimates_nonlinear_propscore_adj <- colMeans(psi_matrix_nonlinear_propscore_adj)
biases_nonlinear_propscore_adj <- psi-estimates_nonlinear_propscore_adj[1:3]
SD_nonlinear_propscore_adj <- sqrt(diag(cov(psi_matrix_nonlinear_propscore_adj)))
CI_nonlinear_propscore_adj <- colQuantiles(psi_matrix_nonlinear_propscore_adj, probs = c(0.025, 0.975))
estimates
estimates_nonlinear
estimates_propscore_adj
estimates_nonlinear_propscore_adj
biases
biases_nonlinear
biases_propscore_adj
biases_nonlinear_propscore_adj

SD
SD_nonlinear
SD_propscore_adj
SD_nonlinear_propscore_adj

CI
CI_nonlinear
CI_propscore_adj
CI_nonlinear_propscore_adj

#Including boxplot for all estimates for all models
library(GGally)
library(reshape2)
library(tidyr)
library(dplyr)
estimates_censoring
psi_dataframe <- as.data.frame(psi_matrix[,1:3])
psi_dataframe_nonlinear <- as.data.frame(psi_matrix_nonlinear[,1:3])
psi_dataframe_propscore_adj <- as.data.frame(psi_matrix_propscore_adj[,1:3])
psi_dataframe_nonlinear_propscore_adj <- as.data.frame(psi_matrix_nonlinear_propscore_adj[,1:3])
colnames(psi_dataframe) <- colnames(psi_dataframe_nonlinear) <- colnames(psi_dataframe_propscore_adj) <- colnames(psi_dataframe_nonlinear_propscore_adj) <- c("psi21", "psi32", "psi31")
psi_dataframe <- psi_dataframe %>%
  mutate(matrix = "Linear outcome not propensity score adjusted")
psi_dataframe_nonlinear <- psi_dataframe_nonlinear %>%
  mutate(matrix = "Non-linear outcome not propensity score adjusted")
psi_dataframe_nonlinear_propscore_adj <- psi_dataframe_nonlinear_propscore_adj %>%
  mutate(matrix = "Non-linear outcome propensity score adjusted")
psi_dataframe_propscore_adj <- psi_dataframe_propscore_adj %>%
  mutate(matrix = "Linear outcome propensity score adjusted")
combined_df <- bind_rows(psi_dataframe, psi_dataframe_nonlinear, psi_dataframe_propscore_adj, psi_dataframe_nonlinear_propscore_adj)
long_df <- combined_df %>% 
  pivot_longer(cols = c("psi21", "psi32", "psi31"),
               names_to = "Parameter",
               values_to = "Estimate")
ggplot(long_df, aes(x = Parameter, y = Estimate, fill = matrix)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Linear outcome not propensity score adjusted" = "red", "Non-linear outcome not propensity score adjusted" = "magenta", "Linear outcome propensity score adjusted" = "purple", "Non-linear outcome propensity score adjusted" = "yellow")) +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.3, 0.4))+
  labs(title = "Comparison of Bootstrapped Estimates",
       x = "Parameter",
       y = "Estimate",
       fill = "Simulation setting and estimation model")+
  geom_hline(aes(yintercept = 0, color = "psi21"), show.legend = TRUE, linetype = 2) + 
  geom_hline(aes(yintercept = 0.1, color = "psi32"), show.legend = TRUE, linetype = 2) +
  geom_hline(aes(yintercept = 0.2, color = "psi31"), show.legend = TRUE, linetype = 2) +
  scale_color_manual(values = c("psi21" = "green", "psi32" = "blue", "psi31" = "orange")) +
  guides(
    color = guide_legend(title = "True parameter values"),
    linetype = "none"
  )

#Censoring with IPCW estimation
n <- 1e5

#3 columns for the three effect estimates, and one column for the presumed biased estimate for psi31
psi_matrix_censoring <- matrix(0, nrow = B, ncol = 4)
for (b in 1:B){
  n <- 1e5
  generated_data <- generate_data(n, psi, Y3L2coef, linear = FALSE, censoring = TRUE,test = TRUE)
  #effect of X2 on Y3
  ipw_c2 <- ipwpoint(exposure = C2, family = "binomial", link = "logit",
                     denominator = ~ X1 + L1 + Y1+C,
                     numerator = ~1,
                     data = generated_data)
  ipw_c3 <- ipwpoint(exposure = C3, family = "binomial", link = "logit",
                     denominator = ~ X2 + L2 + Y2 + X1 + L1 + Y1+C,
                     numerator = ~X1+Y1+L1+C,
                     data = generated_data[generated_data$C2 == 0, ])
  # Combine weights
  
  generated_data$ipw <- 1
  generated_data$ipw[generated_data$C2 == 0] <- ipw_c2$ipw.weights[generated_data$C2==0]
  generated_data$ipw[generated_data$C3==0] <- ipw_c3$ipw.weights[generated_data[generated_data$C2 ==0,]$C3 == 0]
  #Propensity scores
  uncensored_Y2 <- which(!is.na(generated_data$Y2))
  uncensored_Y3 <- which(!is.na(generated_data$Y3))
  p2_model <- glm(X2~X1+Y1+Y2+L1+L2+C, family = binomial("logit"),subset = uncensored_Y3,data = generated_data, na.action = na.exclude)
  p1_model <- glm(X1~Y1 + L1 + C, family = binomial("logit"), subset = uncensored_Y2,data = generated_data, na.action = na.exclude)
  generated_data[uncensored_Y3,"P2"] <- predict.glm(p2_model, type = "response")
  generated_data[uncensored_Y2,"P1"] <- predict.glm(p1_model, type = "response")
  # Effect of X2 on Y3
  model_y3x2 <- lm(Y3 ~ X2 + Y2 + X1 + L1 + L2 + Y1 + C + P2 + P1, 
                   data = generated_data, 
                   weights = ipw, 
                   na.action = na.exclude)
  psi32 <- coef(model_y3x2)["X2"]
  
  # Effect of X1 on Y2
  model_y2x1 <- lm(Y2 ~ X1 + L1 + Y1 + C + P1, 
                   data = generated_data, 
                   weights = ipw, 
                   na.action = na.exclude)
  psi21 <- coef(model_y2x1)["X1"]
  
  # Residualize
  generated_data$R32 <- generated_data$Y3 - psi32 * generated_data$X2
  
  # Final model
  model_y3x1 <- lm(R32 ~ X1 + L1 + Y1 + C + P1, 
                   data = generated_data, 
                   weights = ipw, 
                   na.action = na.exclude)
  psi31 <- coef(model_y3x1)["X1"]
  
  # Combine results
  psi_est <- c(psi21, psi32, psi31)
  names(psi_est) <- c("psi21", "psi32", "psi31")
  psi_est
  psi_matrix_censoring[b,] <- c(psi_est, coef(model_y3x2)["X1"])
}
estimates_censoring <- colMeans(psi_matrix_censoring)
bias_censoring <- psi - estimates_censoring[1:3]
SD_censoring <- sqrt(diag(cov(psi_matrix_censoring)))
CI_censoring <- colQuantiles(psi_matrix_censoring, probs = c(0.025, 0.975))
estimates_censoring
bias_censoring
SD_censoring
CI_censoring

#Censoring naive implementation without IPCW estimation
n <- 1e5

#3 columns for the three effect estimates, and one column for the presumed biased estimate for psi31
psi_matrix_censoring_naive <- matrix(0, nrow = B, ncol = 4)
for (b in 1:B){
  generated_data <- generate_data(n, psi, Y3L2coef, linear = FALSE, censoring = TRUE,test = TRUE)
  #Propensity scores
  uncensored_Y2 <- which(!is.na(generated_data$Y2))
  uncensored_Y3 <- which(!is.na(generated_data$Y3))
  p2_model <- glm(X2~X1+Y1+Y2+L1+L2+C, family = binomial("logit"),subset = uncensored_Y3,data = generated_data, na.action = na.exclude)
  p1_model <- glm(X1~Y1 + L1 + C, family = binomial("logit"), subset = uncensored_Y2,data = generated_data, na.action = na.exclude)
  generated_data[uncensored_Y3,"P2"] <- predict.glm(p2_model, type = "response")
  generated_data[uncensored_Y2,"P1"] <- predict.glm(p1_model, type = "response")
  # Effect of X2 on Y3
  model_y3x2 <- lm(Y3 ~ X2 + Y2 + X1 + L1 + L2 + Y1 + C + P2 + P1, 
                   data = generated_data, 
                   na.action = na.exclude)
  psi32 <- coef(model_y3x2)["X2"]
  
  # Effect of X1 on Y2
  model_y2x1 <- lm(Y2 ~ X1 + L1 + Y1 + C + P1, 
                   data = generated_data, 
                   na.action = na.exclude)
  psi21 <- coef(model_y2x1)["X1"]
  
  # Residualize
  generated_data$R32 <- generated_data$Y3 - psi32 * generated_data$X2
  
  # Final model
  model_y3x1 <- lm(R32 ~ X1 + L1 + Y1 + C + P1, 
                   data = generated_data,
                   na.action = na.exclude)
  psi31 <- coef(model_y3x1)["X1"]
  
  # Combine results
  psi_est <- c(psi21, psi32, psi31)
  names(psi_est) <- c("psi21", "psi32", "psi31")
  psi_est
  psi_matrix_censoring_naive[b,] <- c(psi_est, coef(model_y3x2)["X1"])
}
estimates_censoring_naive <- colMeans(psi_matrix_censoring_naive)
bias_censoring_naive <- psi - estimates_censoring_naive[1:3]
SD_censoring_naive <- sqrt(diag(cov(psi_matrix_censoring_naive)))
CI_censoring_naive <- colQuantiles(psi_matrix_censoring_naive, probs = c(0.025, 0.975))

library(GGally)
library(reshape2)
library(tidyr)
library(dplyr)
estimates_censoring
psi_dataframe_censoring <- as.data.frame(psi_matrix_censoring[,1:3])
psi_dataframe_censoring_naive <- as.data.frame(psi_matrix_censoring_naive[,1:3])

colnames(psi_dataframe_censoring) <- c("psi21", "psi32", "psi31")
colnames(psi_dataframe_censoring_naive) <- c("psi21", "psi32", "psi31")
psi_dataframe_censoring <- psi_dataframe_censoring %>%
  mutate(matrix = "IPCW censoring estimation")
psi_dataframe_censoring_naive <- psi_dataframe_censoring_naive %>%
  mutate(matrix = "Naive censoring estimation")
combined_df <- bind_rows(psi_dataframe_censoring, psi_dataframe_censoring_naive)
long_df <- combined_df %>% 
  pivot_longer(cols = c("psi21", "psi32", "psi31"),
               names_to = "Parameter",
               values_to = "Estimate")
ggplot(long_df, aes(x = Parameter, y = Estimate, fill = matrix)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("IPCW censoring estimation" = "magenta", "Naive censoring estimation" = "red")) +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.125, 0.3)) +
  labs(title = "Comparison of Bootstrapped Estimates",
       x = "Parameter",
       y = "Estimate",
       fill = "Censoring adjustment") +
  geom_hline(aes(yintercept = 0, color = "psi21"), show.legend = TRUE, linetype = 2) + 
  geom_hline(aes(yintercept = 0.1, color = "psi32"), show.legend = TRUE, linetype = 2) +
  geom_hline(aes(yintercept = 0.2, color = "psi31"), show.legend = TRUE, linetype = 2) +
  scale_color_manual(values = c("psi21" = "green", "psi32" = "blue", "psi31" = "orange")) +
  guides(
    color = guide_legend(title = "True parameter values"),
    linetype = "none"
  )
estimates_censoring
estimates_censoring_naive
bias_censoring
bias_censoring_naive
SD_censoring
SD_censoring_naive
CI_censoring
CI_censoring_naive

#testing
n <- 1e5
source("data_generation_tutorial.R")
generated_data_test <- generate_data_test(100, psi, Y3L2coef, linear = FALSE, censoring = TRUE,test = TRUE)
data <- FormatData(generated_data_test, idvar = "id",An = "A", varying = c("A", "L", "Y"),timevar = "time", Cn = "C",
                   GenerateHistory = TRUE, GenerateHistoryMax = 1)
ipw <- ipwtm(exposure = C, family = "binomial",data = generated_data_test[generated_data_test], id = id, type = "all",timevar = time,link = "logit", denominator = ~A+L+Y+Y_base + C_base)
#effect of X2 on Y3
ipw_c2 <- ipwpoint(exposure = C2, family = "binomial", link = "logit",
                   denominator = ~ X1 + L1 + Y1+C,
                   numerator = ~1,
                   data = generated_data)
ipw_c3 <- ipwpoint(exposure = C3, family = "binomial", link = "logit",
                   denominator = ~ X2 + L2 + Y2 + X1 + L1 + Y1+C,
                   numerator = ~1,
                   data = generated_data[generated_data$C2 == 0, ])
# Combine weights

generated_data$ipw2 <- 1
generated_data$ipw3 <- 1
generated_data$ipw2[generated_data$C2 == 0] <- ipw_c2$ipw.weights[generated_data$C2==0]
generated_data$ipw3[generated_data$C3==0] <- ipw_c3$ipw.weights[generated_data[generated_data$C2 ==0,]$C3 == 0]
#Propensity scores
uncensored_Y2 <- which(!is.na(generated_data$Y2))
uncensored_Y3 <- which(!is.na(generated_data$Y3))
p2_model <- glm(X2~X1+Y1+Y2+L1+L2+C, family = binomial("logit"),subset = uncensored_Y3,data = generated_data, na.action = na.exclude)
p1_model <- glm(X1~Y1 + L1 + C, family = binomial("logit"), subset = uncensored_Y2,data = generated_data, na.action = na.exclude)
generated_data[uncensored_Y3,"P2"] <- predict.glm(p2_model, type = "response")
generated_data[uncensored_Y2,"P1"] <- predict.glm(p1_model, type = "response")
# Effect of X2 on Y3
model_y3x2 <- lm(Y3 ~ X2 + Y2 + X1 + L1 + L2 + Y1 + C + P2 + P1, 
                 data = generated_data, 
                 weights = ipw3, 
                 na.action = na.exclude)
psi32 <- coef(model_y3x2)["X2"]

# Effect of X1 on Y2
model_y2x1 <- lm(Y2 ~ X1 + L1 + Y1 + C + P1, 
                 data = generated_data, 
                 weights = ipw2, 
                 na.action = na.exclude)
psi21 <- coef(model_y2x1)["X1"]

# Residualize
generated_data$R32 <- generated_data$Y3 - psi32 * generated_data$X2

# Final model
model_y3x1 <- lm(R32 ~ X1 + L1 + Y1 + C + P1, 
                 data = generated_data, 
                 weights = ipw3, 
                 na.action = na.exclude)
psi31 <- coef(model_y3x1)["X1"]

# Combine results
psi_est <- c(psi21, psi32, psi31)
names(psi_est) <- c("psi21", "psi32", "psi31")
psi_est
psi_matrix_censoring[b,] <- c(psi_est, coef(model_y3x2)["X1"])
