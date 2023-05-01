stressdata <- read.csv("/home/cc/homework/biostat/stressEcho.csv")
layout(matrix(1:10, nrow = 2))

summary(stressdata)

# Bivariate analysis: Correlation matrix
cor_matrix <- cor(stressdata[, sapply(stressdata, is.numeric)])
print(cor_matrix)

# argumentnames <- c("bhr", "basebp", "basedp", "pkhr", "sbp", "dp", "dose", "maxhr", "pctMphr", "mbp", "dpmaxdo", "dobdose", "age", "gender", "baseEF", "dobEF", "chestpain", "restwma", "posSE", "hxofHT", "hxofDM", "hxofCig", "hxofMI", "hxofPTCA", "hxofCABG", "ecg")

# Univariable logistic regression for each predictor
univariable_models <- lapply(argumentnames, function(predictor) {
  formula <- paste("any.event ~", predictor)
  model <- glm(formula, data = stressdata, family = "binomial")
  return(model)
})

# Print summary for each model
lapply(univariable_models, summary)

# Extract p-values for each coefficient (excluding the intercept) in each model
p_values_list <- lapply(univariable_models, function(model) {
  coefficients <- summary(model)$coefficients[-1, "Pr(>|z|)"] # Exclude the intercept
  return(coefficients)
})

# Combine variable names and p-values into a single data frame
p_values_df <- data.frame()
for (i in seq_along(argumentnames)) {
  predictor_name <- ifelse(is.null(names(p_values_list[[i]])), argumentnames[i], names(p_values_list[[i]]))
  predictor_p_values <- data.frame(variable = predictor_name, p_value = p_values_list[[i]])
  p_values_df <- rbind(p_values_df, predictor_p_values)
}

# Filter data frame to get variables with p-value < 0.1
significant_predictors <- p_values_df[p_values_df$p_value < 0.05, ]

print(significant_predictors)






# Load the necessary library
library(MASS)

# Define the significant predictor variables (replace 'ecgMI' with 'ecg')
significant_predictors <- c("dp", "dpmaxdo", "baseEF", "dobEF", "restwma", "posSE", "hxofHT", "hxofDM", "hxofMI", "ecg")

# Perform multivariable logistic regression using significant predictors
formula <- paste("any.event ~", paste(significant_predictors, collapse = " + "))
multivariable_model <- glm(formula, data = stressdata, family = "binomial")

# Print the summary of the multivariable model
summary(multivariable_model)

# Model selection using stepwise selection
stepwise_model <- stepAIC(multivariable_model, direction = "both", trace = FALSE)

# Print the summary of the final model
summary(stepwise_model)








hist(stressdata$dp, main = "Age Distribution", xlab = "Age", col = "lightblue")
boxplot(stressdata$dp, main = "Age Boxplot", ylab = "Age", col = "lightblue")

hist(stressdata$pctMphr, main = "pctMphr Distribution", xlab = "pctMphr", col = "lightblue")
boxplot(stressdata$pctMphr, main = "pctMphr Boxplot", ylab = "pctMphr", col = "lightblue")

hist(stressdata$dobEF, main = "dobEF Distribution", xlab = "dobEF", col = "lightblue")
boxplot(stressdata$hxofDM, main = "dobEF Boxplot", ylab = "dobEF", col = "lightblue")

plot(stressdata$age, stressdata$any.event, main = "Age vs. any.event", xlab = "Age", ylab = "any.event", col = "blue")
boxplot(ecg_MI ~ any.event, data = stressdata, main = "Age Distribution by Any Event", xlab = "Any Event", ylab = "Age", col = c("lightblue", "lightgreen"))

boxplot(dp ~ any.event, data = stressdata, main = "Age Distribution by Any Event", xlab = "Any Event", ylab = "Age", col = c("lightblue", "lightgreen"))

boxplot(hxofDM ~ any.event, data = stressdata, main = "Age Distribution by Any Event", xlab = "Any Event", ylab = "Age", col = c("lightblue", "lightgreen"))







# fit a logistic regression model without newMI, newPTCA, newCABG, and death
model <- glm(any.event ~ . - newMI - newPTCA - newCABG - death -X, data = stressdata, family = binomial())

# print the model summary
summary(model)


# subset the stressdata into two groups based on the values of 'any.event'
group1 <- subset(stressdata, any.event == 0)$newMI
group2 <- subset(stressdata, any.event == 1)$newMI

# perform two-sample z-test
ztest <- t.test(group1, group2, var.equal = TRUE)
print(ztest)

# extract the test statistic and p-value
test_statistic <- ztest$statistic
p_value <- ztest$p.value