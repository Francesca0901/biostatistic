stressdata <- read.csv("/home/cc/homework/biostat/stressEcho.csv")
layout(matrix(1, nrow = 1))

summary(stressdata)

# Bivariate analysis: Correlation matrix
cor_matrix <- cor(stressdata[, sapply(stressdata, is.numeric)])
print(cor_matrix)

argumentnames <- c("bhr", "basebp", "basedp", "pkhr", "sbp", "dp", "dose", "maxhr", "pctMphr", "mbp", "dpmaxdo", "dobdose", "age", "gender", "baseEF", "dobEF", "chestpain", "restwma", "posSE", "hxofHT", "hxofDM", "hxofCig", "hxofMI", "hxofPTCA", "hxofCABG", "ecg")

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

stressdata$ecgMI <- as.integer(stressdata$ecg == "MI")

# Define the significant predictor variables
significant_predictors <- c("dp", "dpmaxdo", "baseEF", "dobEF", "restwma", "posSE", "hxofHT", "hxofDM", "hxofMI", "ecgMI")

# Perform multivariable logistic regression using significant predictors and selected interaction terms
formula <- paste("any.event ~", paste(significant_predictors, collapse = " + "), "+ restwma:ecgMI + restwma:ecgMI:posSE")
multivariable_model_interactions <- glm(formula, data = stressdata, family = "binomial")

# Print the summary of the multivariable model with selected interaction terms
summary(multivariable_model_interactions)

# Model selection using stepwise selection
stepwise_model_interactions <- stepAIC(multivariable_model_interactions, direction = "both", trace = FALSE)

# Print the summary of the final model with selected interaction terms
summary(stepwise_model_interactions)



library(ResourceSelection)
# Assuming your stepwise model is stored in a variable called 'stepwise_model'
hoslem_test <- hoslem.test(stepwise_model_interactions$y, fitted(stepwise_model_interactions), g = 10)

# Print the results
print(hoslem_test)

library(pROC)
# Get the predicted probabilities from the stepwise model
predicted_probabilities <- predict(stepwise_model_interactions, type = "response")

# Create the ROC curve
roc_obj <- roc(stressdata$any.event, predicted_probabilities)

# Calculate the AUC
auc(roc_obj)

library(dplyr)
library(tidyr)
# Define the predictor variable names
final_predictors <- c("dp", "dpmaxdo", "baseEF", "dobEF", "restwma", "posSE", "hxofHT", "hxofDM", "hxofMI", "ecgMI")

# Select only significant predictor variables
mydata <- stressdata[, c(significant_predictors, "any.event")]
mydata <- mydata %>%
  mutate(logit = log(predicted_prob/(1-predicted_prob))) %>%
  gather(key = "predictors", value = "predictor.value", significant_predictors)

ggplot(mydata, aes(logit, predictor.value))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() + 
  facet_wrap(~predictors, scales = "free_y")




# Calculate the predicted probabilities
predicted_prob <- predict(stepwise_model_interactions, type = "response")

# Calculate the logit values
logit_predicted_prob <- log(predicted_prob / (1 - predicted_prob))
# Set the layout for the scatterplots
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))

predictors <- c("baseEF", "dobEF", "dp", "dpmaxdo")
# Create scatterplots of logit_predicted_prob vs. predictors
for (predictor in significant_predictors) {
  plot(logit_predicted_prob ~ stressdata[[predictor]], xlab = predictor, ylab = "Logit of Predicted Probability")
}

library(ggplot2)
# Create the scatter plots:
ggplot(mydata, aes(logit, predictor.value))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() + 
  facet_wrap(~predictors, scales = "free_y")





# Calculate Cook's distance
cooks_dist <- cooks.distance(stepwise_model_interactions)

# Calculate standardized residuals
std_resid <- rstandard(stepwise_model_interactions)

# Plot Cook's distance
plot(cooks_dist, type = "h", main = "Cook's Distance", xlab = "Observation", ylab = "Cook's Distance")

# Plot standardized residuals
plot(std_resid, type = "h", main = "Standardized Residuals", xlab = "Observation", ylab = "Standardized Residuals")


library(broom)
# Extract model results
model.data <- augment(stepwise_model_interactions) %>% 
  mutate(index = 1:n()) 

model.data %>% top_n(3, .cooksd)


ggplot(model.data, aes(index, .std.resid)) + 
  geom_point(aes(color = any.event), alpha = .5) +
  geom_point(data = model.data %>% filter(abs(.std.resid) > 3), aes(color = "red", shape = "x"), size = 3) +
  theme_bw()




avg_age <- mean(stressdata$dobEF)

# calculate number and percentage of people above average age with any.event
above_avg_age <- stressdata$dobEF > avg_age
below_avg_age <- stressdata$dobEF < avg_age
n_above = sum(above_avg_age)
n_below = sum(below_avg_age)
n_above_avg_with_event <- sum(stressdata$any.event == 1 & above_avg_age)
n_below_avg_with_event <- sum(stressdata$any.event == 1 & below_avg_age)
n_above_avg_without_event <- sum(stressdata$any.event == 0 & above_avg_age)
n_below_avg_without_event <- sum(stressdata$any.event == 0 & below_avg_age)
pct_above_avg_with_event <- n_above_avg_with_event / sum(above_avg_age) * 100
pct_below_avg_with_event <- n_below_avg_with_event / sum(below_avg_age) * 100
pct_above_avg_without_event <- n_above_avg_without_event / sum(above_avg_age) * 100
pct_below_avg_without_event <- n_below_avg_without_event / sum(below_avg_age) * 100





# hist(stressdata$dp, main = "Age Distribution", xlab = "Age", col = "lightblue")
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