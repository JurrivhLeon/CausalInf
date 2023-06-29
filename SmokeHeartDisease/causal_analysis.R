# Import preliminary packages.
options(warn = -1)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(grid)
library(MatchIt)
library(tableone)
library(survey)

# Read in the data.
setwd('D:\\Users\\Jurrivh Liao\\Documents\\All of statistics\\Causal Inference\\Final')
smht <- read.csv('smoke_heart.csv')

# Compare the covariates.
bin.total <-  data.frame(
  avg = apply(smht[, c(2, 4, 5, 6, 8, 13)], 2, sum)
)
cont.total <- data.frame(
  avg = apply(smht[, c(3, 9, 10, 11)], 2, mean),
  std = apply(smht[, c(3, 9, 10, 11)], 2, sd)
)
table(smht$work_type)

smht.smoke <- smht[smht$smoking_status == 1,]
bin.smoke <-  data.frame(
  avg = apply(smht.smoke[, c(2, 4, 5, 6, 8, 13)], 2, sum)
)
cont.smoke <- data.frame(
  avg = apply(smht.smoke[, c(3, 9, 10, 11)], 2, mean),
  std = apply(smht.smoke[, c(3, 9, 10, 11)], 2, sd)
)
table(smht.smoke$work_type)

smht.nonsmoke <- smht[smht$smoking_status == 0,]
bin.nonsmoke <-  data.frame(
  avg = apply(smht.nonsmoke[, c(2, 4, 5, 6, 8, 13)], 2, sum)
)
cont.nonsmoke <- data.frame(
  avg = apply(smht.nonsmoke[, c(3, 9, 10, 11)], 2, mean),
  std = apply(smht.nonsmoke[, c(3, 9, 10, 11)], 2, sd)
)
table(smht.nonsmoke$work_type)

# Compute point estimate and confidence interval of OR.
or.est <- function(model){
  point <- model$coefficients['smoking_status']
  sd <- summary(model)$coefficients[,2]['smoking_status']
  or <- exp(point)
  ci <- exp(c(point - qnorm(0.975) * sd, point + qnorm(0.975) * sd))
  return(c(or=or, ci=ci))
}

# Regression.
crude.model <- glm(heart_disease ~ smoking_status, family = binomial, data = smht)
summary(crude.model)
or.est(crude.model)
adj1.model <- glm(heart_disease ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi, 
                family = binomial, data = smht)
summary(adj1.model)
or.est(adj1.model)
adj2.model <- glm(heart_disease ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi + hypertension + stroke + avg_glucose_level + cholesterol, 
                 family = binomial, data = smht)
summary(adj2.model)
or.est(adj2.model)

# Effect modification.
smht.male <- smht[smht$gender == 1,]
smht.female <- smht[smht$gender == 0,]

crude.model.male <- glm(heart_disease ~ smoking_status, family = binomial, data = smht.male)
summary(crude.model.male)
or.est(crude.model.male)
adj1.model.male <- glm(heart_disease ~ smoking_status + age + ever_married + as.factor(work_type) + Residence_type + bmi, 
                  family = binomial, data = smht.male)
summary(adj1.model.male)
or.est(adj1.model.male)
adj2.model.male <- glm(heart_disease ~ smoking_status + age + ever_married + as.factor(work_type) + Residence_type + bmi + hypertension + stroke + avg_glucose_level + cholesterol, 
                  family = binomial, data = smht.male)
summary(adj2.model.male)
or.est(adj2.model.male)

crude.model.female <- glm(heart_disease ~ smoking_status, family = binomial, data = smht.female)
summary(crude.model.female)
or.est(crude.model.female)
adj1.model.female <- glm(heart_disease ~ smoking_status + age + ever_married + as.factor(work_type) + Residence_type + bmi, 
                  family = binomial, data = smht.female)
summary(adj1.model.female)
or.est(adj1.model.female)
adj2.model.female <- glm(heart_disease ~ smoking_status + age + ever_married + as.factor(work_type) + Residence_type + bmi + hypertension + stroke + avg_glucose_level + cholesterol, 
                  family = binomial, data = smht.female)
summary(adj2.model.female)
or.est(adj2.model.female)

# Check interaction.

inter1.model <- glm(heart_disease ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi + gender * smoking_status, 
                   family = binomial, data = smht)
summary(inter1.model)
inter2.model <- glm(heart_disease ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi + + hypertension + stroke + avg_glucose_level + cholesterol + gender * smoking_status, 
                    family = binomial, data = smht)
summary(inter2.model)

# Propensity score model.
ps.model <- glm(smoking_status ~ gender + age + ever_married + as.factor(work_type) + Residence_type + bmi, 
                family = binomial, data = smht)
summary(ps.model)

plt <- ggplot() + 
  geom_histogram(data = data.frame(ps = ps.model$fitted.values[smht$smoking_status == 1]), 
                 bins = 50, aes(x = ps), fill='#0057B7', color="white")+
  geom_histogram(data = data.frame(ps = ps.model$fitted.values[smht$smoking_status == 0]), 
                 bins = 50, aes(x = ps), fill='#FFD700', color="white", alpha=0.5)+
  labs(title = 'Estimated PS distribution', x = 'Propensity Score', y = 'Count') + 
  theme(plot.title = element_text(size=13.5, hjust=0.5, margin=unit(c(0,0,.45,0), 'cm')),
        axis.title.x = element_text(size=12, margin=unit(c(.25,0,0,0), 'cm')), 
        axis.title.y = element_text(size=12, margin=unit(c(0,.45,0,0), 'cm')), 
        axis.text = element_text(size=10.5))
print(plt)

match.model <- matchit(smoking_status ~ gender + age + ever_married + as.factor(work_type) + Residence_type + bmi,
                       method = 'nearest', data=smht, distance = 'glm',
                       ratio=1, replace=FALSE)
smht.matched <- match.data(match.model)
# Draw a jitter plot.
plot(match.model, type='jitter')

# Check the covariate balance.
plotLowess <- function(data, variable, ylabel, show.legend=FALSE) {
  data$variable <- data[, variable]
  data$smoking_status <- as.factor(data$smoking_status)
  support <- c(min(data$variable), max(data$variable))
  ggplot(data, aes(x = distance, y = variable, color = smoking_status)) +
    geom_point(alpha = 0.2, size = 1.2, show.legend=FALSE) +
    geom_smooth(method = "loess", se = F, show.legend=FALSE) +
    scale_color_manual(values = c("#0057B7", "#FFD700")) +
    labs(x = "Propensity score", y = ylabel, col = 'Smoking status') +
    theme_bw() +
    ylim(support)
}

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))
print(plotLowess(smht.matched, 'gender', 'Gender'), 
      vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(plotLowess(smht.matched, 'ever_married', 'Marital status'), 
      vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(plotLowess(smht.matched, 'Residence_type', 'Residence type'), 
      vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(plotLowess(smht.matched, 'age', 'Age'), 
      vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(plotLowess(smht.matched, 'bmi', 'BMI', show.legend=TRUE), 
      vp=viewport(layout.pos.row=2, layout.pos.col=2))

table1 <- CreateTableOne(var=c('gender', 'age', 'bmi', 'ever_married', 'Residence_type', 'work_type'), 
                         strata='smoking_status', data=smht, test = TRUE)
table2 <- CreateTableOne(var=c('gender', 'age', 'bmi', 'ever_married', 'Residence_type', 'work_type'), 
                         strata='smoking_status', data=smht.matched, test=TRUE)
kableone(table1, smd=TRUE)  # Before Matching.
kableone(table2, smd=TRUE)  # After Matching.

match.reg <- glm(heart_disease ~ smoking_status, 
                 family = binomial, data = smht.matched)
summary(match.reg)
or.est(match.reg)

match.adjust <- glm(heart_disease ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi, 
                 family = binomial, data = smht.matched)
summary(match.adjust)
or.est(match.adjust)

# Mediation analysis.
# Mediator detection.
med1.lr <- glm(hypertension ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi,
               family = binomial, data = smht)
summary(med1.lr)

med2.lr <- glm(stroke ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi,
               family = binomial, data = smht)
summary(med2.lr)

med3.lr <- lm(avg_glucose_level ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi,
              data = smht)
summary(med3.lr)

med4.lr <- lm(cholesterol ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi,
              data = smht)
summary(med4.lr)

# Test the normality of mediator.
cholesterol.standardized <- with(smht, (cholesterol - mean(cholesterol))/ sd(cholesterol))
chol.data <- data.frame(smoke = as.factor(smht$smoking_status), chol=smht$cholesterol)
plt.qq <- ggqqplot(chol.data, x = 'chol', color = 'smoke',
                   palette = c("#0057B7", "#FFDD00"), shape='o') +
  grids(linetype='solid') +
  labs(title = 'QQ plot of Cholesterol by Smoking Status', x = 'Normal quantile', 
       y = 'Cholsterol (mmol/L)') +
  theme(plot.title = element_text(size=13.5, hjust=0.5, margin=unit(c(0,0,.45,0), 'cm')),
        axis.title.x = element_text(size=12, margin=unit(c(.25,0,0,0), 'cm')), 
        axis.title.y = element_text(size=12, margin=unit(c(0,.45,0,0), 'cm')), 
        axis.text = element_text(size=10.5))

plt.chol <- ggplot(data=chol.data) +
  geom_density(aes(x=chol, color=smoke), show.legend=FALSE) + 
  stat_density(aes(x=chol, color=smoke), geom="line",position="identity") +
  scale_color_manual(values = c("#0057B7", "#FFD700")) +
  labs(title = 'Density of Cholesterol by Smoking Status', x = 'Cholsterol (mmol/L)', 
       y = 'Density') +
  theme(plot.title = element_text(size=13.5, hjust=0.5, margin=unit(c(0,0,.45,0), 'cm')),
        axis.title.x = element_text(size=12, margin=unit(c(.25,0,0,0), 'cm')), 
        axis.title.y = element_text(size=12, margin=unit(c(0,.45,0,0), 'cm')), 
        axis.text = element_text(size=10.5))

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(plt.qq, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(plt.chol, vp=viewport(layout.pos.row=1, layout.pos.col=2))

# With interaction between mediator & exposure.
med.model <- glm(heart_disease ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi + cholesterol + cholesterol*smoking_status, 
                 family = binomial, data = smht)
summary(med.model)

# Estimate NDE and NIE.
med.est <- function(model_lr, model_lg, interaction=TRUE){
  theta1 <- model_lg$coefficients['smoking_status']
  theta2 <- model_lg$coefficients['cholesterol']
  if(interaction) theta3 <- model_lg$coefficients['smoking_status:cholesterol']
  else theta3 <- 0
  sigmasq <- summary(model_lr)$sigma ^ 2
  beta1 <- model_lr$coefficients['smoking_status']
  fitv <- mean(model_lr$fitted.values) - beta1 * mean(c(model_lr$model$smoking_status))
  lognde <- theta1 + theta3 * (fitv + theta2 * sigmasq) + 0.5 * theta3^2 * sigmasq
  nde <- exp(lognde)
  lognie <- (theta2 + theta3) * beta1
  nie <- exp(lognie)
  te <- nde * nie
  pm <- (nde * (nie - 1)) / (nde * nie - 1)
  return(c(nde=nde, nie=nie, te=te, pm=pm))
}

# Construct a confidence interval for NDE and NIE by bootstrapping.
med.conf.int <- function(dataset, b=1000){
  n <- dim(dataset)[1]
  nde.all <- rep(NaN, b)
  nie.all <- rep(NaN, b)
  te.all <- rep(NaN, b)
  for(i in 1:b){
    bootstrap.index <- sample(1:n, size=n, replace=TRUE)
    boot <- dataset[bootstrap.index,]
    med.lr <- lm(cholesterol ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi,
                  data = boot)
    med.lg <- glm(heart_disease ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi + cholesterol + cholesterol*smoking_status, 
                     family = binomial, data = boot)
    eff <- med.est(med.lr, med.lg)
    nde.all[i] <- eff[1]
    nie.all[i] <- eff[2]
    te.all[i] <- eff[3]
  }
  nde.conf <- quantile(nde.all, c(0.025, 0.975))
  nie.conf <- quantile(nie.all, c(0.025, 0.975))
  te.conf <- quantile(te.all, c(0.025, 0.975))
  return(data.frame(nde.conf=nde.conf, nie.conf=nie.conf, te.conf=te.conf))
}

set.seed(2023)
med.est(med4.lr, med.model)
med.conf <- med.conf.int(smht) # This takes about 30s to run.

# No interaction between mediator & exposure.
med.model_ <- glm(heart_disease ~ smoking_status + gender + age + ever_married + as.factor(work_type) + Residence_type + bmi + cholesterol, 
                 family = binomial, data = smht)
summary(med.model_)
med.est(med4.lr, med.model_, interaction=FALSE)
med.conf_ <- med.conf.int(smht) # This takes about 30s to run.

# Confoundings of gender.
rrau.gender <- exp(ps.model$coefficients['gender'])
rruy.gernder <- exp(adj1.model$coefficients['gender'])

