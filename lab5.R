data("plasma", package = "HSAUR3")
layout(matrix(1:2, nrow = 2))
cdplot(ESR ~ fibrinogen, data = plasma)
cdplot(ESR ~ globulin, data = plasma)

plasma.glm.1 <- glm(ESR ~ fibrinogen, data = plasma, family = binomial())

summary(plasma.glm.1)
confint(plasma.glm.1, parm = "fibrinogen")

exp(coef(plasma.glm.1)["fibrinogen"])
exp(confint(plasma.glm.1, parm = "fibrinogen"))


plasma.glm.2 <- glm(ESR ~ fibrinogen + globulin, data = plasma, family = binomial())
summary(plasma.glm.2)
confint(plasma.glm.2, parm = "fibrinogen")

exp(coef(plasma.glm.2)["fibrinogen"])
exp(confint(plasma.glm.2, parm = "fibrinogen"))

anova(plasma.glm.1, plasma.glm.2, test = "Chisq")

prob <- predict(plasma.glm.2, type = "response")

plot(globulin ~ fibrinogen, data = plasma, xlim = c(2,6), ylim = c(25,55), pch = "*")
symbols(plasma$fibrinogen, plasma$globulin, circles = prob, add = TRUE)

