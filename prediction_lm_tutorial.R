setwd("C:/Users/marcu/OneDrive/Dokumenter/Master/masters-master/")
generated_data <- read.csv(file = "generated_data.csv", header = TRUE)
#effect of X2 on Y3
model_y3x2 <- lm(Y3 ~ X2+X1+Y2+Y1+L1+L2+C, data = generated_data)
psi32 <- coef(model_y3x2)["X2"]
#residualise 
generated_data[,"R32"] <- generated_data[,"Y3"]-psi32*generated_data[,"X2"]
#effect of X1 on Y3
model_y3x1 <- lm(R32 ~ C+Y1+L1+X1, data = generated_data)
psi31 <- coef(model_y3x1)["X1"]

#effect of X1 on Y2

model_y2x1 <- lm(Y2~X1+L1+Y1+C, data = generated_data)
psi21 <- coef(model_y2x1)["X1"]

psi <- c(psi21, psi32, psi31)
names(psi) <- c("psi21", "psi32", "psi31")
psi
