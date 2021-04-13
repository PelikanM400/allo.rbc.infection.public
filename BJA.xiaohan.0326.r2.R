########################################################
# preparations before analysis 
# set working folder
setwd("")

# load necessary libraries 
if(!require('MatchIt')) install.packages('MatchIt')
if(!require('stddiff')) install.packages('stddiff')
if(!require('Rcpp')) install.packages('Rcpp')
if(!require('car')) install.packages('car')
if(!require('MASS')) install.packages('MASS')
if(!require('brant')) install.packages('brant')
if(!require('metafor')) install.packages('metafor')
library('MatchIt')
library('stddiff')
library('MASS')
library('brant')
library('car')
library('metafor')

# load data 
df <- read.csv('df.analysis.final.1025.csv')
df <- read.csv('df.analysis.final.202101.csv', encoding = 'UTF-8')
head(df)

#suppl. baseline
table(df$infect, df$group)
round(prop.table(table(df$infect, df$group), 2), 3)
sum(df$sum)

mean(df$age)
sd(df$age)
quantile(df$age)
aggregate(df$age, by = list(df$group), FUN = quantile)

table(df$sex)
round(prop.table(table(df$sex)), 3)

mean(df$BMI)
sd(df$BMI)

table(df$asa.clas)
round(prop.table(table(df$asa.clas)), 3)

table(df$HTN)
round(prop.table(table(df$HTN)), 3)

table(df$DM)
round(prop.table(table(df$DM)), 3)

table(df$emer.)
round(prop.table(table(df$emer.)), 3)

table(df$anes.gener)
round(prop.table(table(df$anes.gener)), 3)

hist(df$surg.dura)
mean(df$surg.dura, na.rm = TRUE)
sd(df$surg.dura, na.rm = TRUE)

table(df$smoke)
round(prop.table(table(df$smoke)), 3)

table(df$Af)
round(prop.table(table(df$Af)), 3)

table(df$COPD)
round(prop.table(table(df$COPD)), 3)

table(df$CI)
round(prop.table(table(df$CI)), 3)

table(df$SAH)
round(prop.table(table(df$SAH)), 3)

table(df$radio)
round(prop.table(table(df$radio)), 3)

table(df$immu)
round(prop.table(table(df$immu)), 3)

table(df$chemo)
round(prop.table(table(df$chemo)), 3)

table(df$CAD)
round(prop.table(table(df$CAD)), 3)

table(df$anem.)
round(prop.table(table(df$anem.)), 3)

table(df$open)
round(prop.table(table(df$open)), 3)

table(df$lower)
round(prop.table(table(df$lower)), 3)

mean(df$ebl.ratio)
quantile(df$ebl.ratio)

table(df$cancer)
round(prop.table(table(df$cancer)), 3)
table(df$cancer, df$group)
round(prop.table(table(df$cancer, df$group), 2), 3)
chisq.test(df$cancer, df$group)
stddiff.binary(data = df, 
               gcol = which(colnames(df) == 'group'), 
               vcol = which(colnames(df) == 'cancer'))

table(df$coagu)
round(prop.table(table(df$coagu)), 3)
table(df$coagu, df$group)
round(prop.table(table(df$coagu, df$group), 2), 3)
chisq.test(df$coagu, df$group)
stddiff.binary(data = df, 
               gcol = which(colnames(df) == 'group'), 
               vcol = which(colnames(df) == 'coagu'))

table(df$liver)
round(prop.table(table(df$liver)), 3)
table(df$liver, df$group)
round(prop.table(table(df$liver, df$group), 2), 3)
chisq.test(df$liver, df$group)
stddiff.binary(data = df, 
               gcol = which(colnames(df) == 'group'), 
               vcol = which(colnames(df) == 'liver'))

table(df$renal)
round(prop.table(table(df$renal)), 3)
table(df$renal, df$group)
round(prop.table(table(df$renal, df$group), 2), 3)
chisq.test(df$renal, df$group)
stddiff.binary(data = df, 
               gcol = which(colnames(df) == 'group'), 
               vcol = which(colnames(df) == 'renal'))

table(df$alb)
round(prop.table(table(df$alb)), 3)
table(df$alb, df$group)
round(prop.table(table(df$alb, df$group), 2), 3)
chisq.test(df$alb, df$group)
stddiff.binary(data = df, 
               gcol = which(colnames(df) == 'group'), 
               vcol = which(colnames(df) == 'alb'))

table(df$anem.)
round(prop.table(table(df$anem.)), 3)
table(df$alb, df$group)
round(prop.table(table(df$anem., df$group), 2), 3)
chisq.test(df$anem., df$group)
stddiff.binary(data = df, 
               gcol = which(colnames(df) == 'group'), 
               vcol = which(colnames(df) == 'anem.'))


# primary outcome, logistic regression
model1 <- glm(infect ~ group, 
              family = binomial('logit'), 
              data = df)
summary(model1) # P values 
exp(summary(model1)[['coefficients']][, 1]) # ORs
exp(confint(model1)) # 95% CIs of ORs
p <- summary(model1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model1)))
ci
or <- data.frame(exp(summary(model1)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results1.infect.csv')

model2 <- glm(infect ~ group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu, 
              family = binomial('logit'), 
              data = df)
summary(model2) # P values 
exp(summary(model2)[['coefficients']][, 1]) # ORs
exp(confint(model2)) # 95% CIs of ORs
p <- summary(model2)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model2)))
ci
or <- data.frame(exp(summary(model2)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results2.infect.csv')

model3 <- glm(infect ~ group + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,              
              family = binomial('logit'), 
              data = df)
summary(model3) # P values 
exp(summary(model3)[['coefficients']][, 1]) # ORs
exp(confint(model3)) # 95% CIs of ORs
p <- summary(model3)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model3)))
ci
or <- data.frame(exp(summary(model3)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results3.infect.csv')

model4 <- glm(infect ~group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.infect.csv')

# primary outcome, ps
# section 1: ps match
df.full <- subset(df, select = - c(trans.gap,trans.other,trans.leuko,trans.age,trans.timing))
match1.out <- matchit(group ~ age + sex + BMI + asa.clas + emer. + anes.gener + smoke + surg.dura + 
                        + Af + COPD + CI + SAH + radio + immu + chemo + HTN + DM + CAD + anem. + coagu 
                      + liver + renal + open + lower + ebl.ratio + cancer + alb, 
                      method = "nearest", 
                      distance = "logit", 
                      ratio = 1, 
                      caliper = 0.05, 
                      data = df.full) 
summary(match1.out)
plot(match1.out, type = "jitter")
plot(match1.out, type = "hist")
df.match1 <- match.data(match1.out)
table (df.match1$infect,df.match1$group)
write.csv(df.match1, "df.match1.csv")

aggregate(df.match1$age, by = list(df.match1$group), FUN = mean)
aggregate(df.match1$age, by = list(df.match1$group), FUN = sd)
aggregate(df.match1$age, by = list(df.match1$group), FUN = quantile)
t.test(age ~ group, data = df.match1)
stddiff.numeric(data = df.match1, 
                gcol = which(colnames(df.match1) == 'group'), 
                vcol = which(colnames(df.match1) == 'age'))

table(df.match1$sex)
round(prop.table(table(df$sex)), 3)
table(df.match1$sex, df.match1$group)
round(prop.table(table(df.match1$sex, df.match1$group), 2), 3)
chisq.test(df.match1$sex, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'sex'))

aggregate(df.match1$BMI, by = list(df.match1$group), FUN = mean)
aggregate(df.match1$BMI, by = list(df.match1$group), FUN = sd)
aggregate(df.match1$BMI, by = list(df.match1$group), FUN = quantile)
t.test(BMI ~ group, data = df.match1)
stddiff.numeric(data = df.match1, 
                gcol = which(colnames(df.match1) == 'group'), 
                vcol = which(colnames(df.match1) == 'BMI'))

table(df.match1$asa.clas)
round(prop.table(table(df$asa.clas)), 3)
table(df.match1$asa.clas, df.match1$group)
round(prop.table(table(df.match1$asa.clas, df.match1$group), 2), 3)
chisq.test(df.match1$asa.clas, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'asa.clas'))

table(df.match1$HTN)
round(prop.table(table(df.match1$HTN)), 3)
table(df.match1$HTN, df.match1$group)
round(prop.table(table(df.match1$HTN, df.match1$group), 2), 3)
chisq.test(df.match1$HTN, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'HTN'))

table(df.match1$DM)
round(prop.table(table(df.match1$DM)), 3)
table(df.match1$DM, df.match1$group)
round(prop.table(table(df.match1$DM, df.match1$group), 2), 3)
chisq.test(df.match1$DM, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'DM'))

table(df.match1$CAD)
round(prop.table(table(df.match1$CAD)), 3)
table(df.match1$CAD, df.match1$group)
round(prop.table(table(df.match1$CAD, df.match1$group), 2), 3)
chisq.test(df.match1$CAD, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'CAD'))

table(df.match1$Af)
round(prop.table(table(df.match1$Af)), 3)
table(df.match1$Af, df.match1$group)
round(prop.table(table(df.match1$Af, df.match1$group), 2), 3)
chisq.test(df.match1$Af, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'Af'))

table(df.match1$CI)
round(prop.table(table(df.match1$CI)), 3)
table(df.match1$CI, df.match1$group)
round(prop.table(table(df.match1$CI, df.match1$group), 2), 3)
chisq.test(df.match1$CI, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'CI'))

table(df.match1$SAH)
round(prop.table(table(df.match1$SAH)), 3)
table(df.match1$SAH, df.match1$group)
round(prop.table(table(df.match1$SAH, df.match1$group), 2), 3)
chisq.test(df.match1$SAH, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'SAH'))

table(df.match1$COPD)
round(prop.table(table(df.match1$COPD)), 3)
table(df.match1$COPD, df.match1$group)
round(prop.table(table(df.match1$COPD, df.match1$group), 2), 3)
chisq.test(df.match1$COPD, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'COPD'))

table(df.match1$smoke)
round(prop.table(table(df.match1$smoke)), 3)
table(df.match1$smoke, df.match1$group)
round(prop.table(table(df.match1$smoke, df.match1$group), 2), 3)
chisq.test(df.match1$smoke, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'smoke'))

table(df.match1$anem.)
round(prop.table(table(df.match1$anem.)), 3)
table(df.match1$anem., df.match1$group)
round(prop.table(table(df.match1$anem., df.match1$group), 2), 3)
chisq.test(df.match1$anem., df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'anem.'))

table(df.match1$coagu)
round(prop.table(table(df.match1$coagu)), 3)
table(df.match1$coagu, df.match1$group)
round(prop.table(table(df.match1$coagu, df.match1$group), 2), 3)
chisq.test(df.match1$coagu, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'coagu'))

table(df.match1$liver)
round(prop.table(table(df.match1$liver)), 3)
table(df.match1$liver, df.match1$group)
round(prop.table(table(df.match1$liver, df.match1$group), 2), 3)
chisq.test(df.match1$liver, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'liver'))

table(df.match1$renal)
round(prop.table(table(df.match1$renal)), 3)
table(df.match1$renal, df.match1$group)
round(prop.table(table(df.match1$renal, df.match1$group), 2), 3)
chisq.test(df.match1$renal, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'renal'))

table(df.match1$alb)
round(prop.table(table(df.match1$alb)), 3)
table(df.match1$alb, df.match1$group)
round(prop.table(table(df.match1$alb, df.match1$group), 2), 3)
chisq.test(df.match1$alb, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'alb'))

table(df.match1$chemo)
round(prop.table(table(df.match1$chemo)), 3)
table(df.match1$chemo, df.match1$group)
round(prop.table(table(df.match1$chemo, df.match1$group), 2), 3)
chisq.test(df.match1$chemo, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'chemo'))

table(df.match1$cancer)
round(prop.table(table(df.match1$cancer)), 3)
table(df.match1$cancer, df.match1$group)
round(prop.table(table(df.match1$cancer, df.match1$group), 2), 3)
chisq.test(df.match1$cancer, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'cancer'))

table(df.match1$radio)
round(prop.table(table(df.match1$radio)), 3)
table(df.match1$radio, df.match1$group)
round(prop.table(table(df.match1$radio, df.match1$group), 2), 3)
chisq.test(df.match1$radio, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'radio'))

table(df.match1$immu)
round(prop.table(table(df.match1$immu)), 3)
table(df.match1$immu, df.match1$group)
round(prop.table(table(df.match1$immu, df.match1$group), 2), 3)
chisq.test(df.match1$immu, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'immu'))

table(df.match1$anes.gener)
round(prop.table(table(df.match1$anes.gener)), 3)
table(df.match1$anes.gener, df.match1$group)
round(prop.table(table(df.match1$anes.gener, df.match1$group), 2), 3)
chisq.test(df.match1$anes.gener, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'anes.gener'))

aggregate(df.match1$surg.dura, by = list(df.match1$group), FUN = mean)
aggregate(df.match1$surg.dura, by = list(df.match1$group), FUN = sd)
aggregate(df.match1$surg.dura, by = list(df.match1$group), FUN = quantile)
t.test(surg.dura ~ group, data = df.match1)
stddiff.numeric(data = df.match1, 
                gcol = which(colnames(df.match1) == 'group'), 
                vcol = which(colnames(df.match1) == 'surg.dura'))

table(df.match1$open)
round(prop.table(table(df.match1$open)), 3)
table(df.match1$open, df.match1$group)
round(prop.table(table(df.match1$open, df.match1$group), 2), 3)
chisq.test(df.match1$open, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'open'))

table(df.match1$lower)
round(prop.table(table(df.match1$lower)), 3)
table(df.match1$lower, df.match1$group)
round(prop.table(table(df.match1$lower, df.match1$group), 2), 3)
chisq.test(df.match1$lower, df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'lower'))

table(df.match1$emer.)
round(prop.table(table(df.match1$emer.)), 3)
table(df.match1$emer., df.match1$group)
round(prop.table(table(df.match1$emer., df.match1$group), 2), 3)
chisq.test(df.match1$emer., df.match1$group)
stddiff.binary(data = df.match1, 
               gcol = which(colnames(df.match1) == 'group'), 
               vcol = which(colnames(df.match1) == 'emer.'))

aggregate(df.match1$ebl.ratio, by = list(df.match1$group), FUN = mean)
aggregate(df.match1$ebl.ratio, by = list(df.match1$group), FUN = sd)
aggregate(df.match1$ebl.ratio, by = list(df.match1$group), FUN = quantile)
t.test(ebl.ratio ~ group, data = df.match1)
stddiff.numeric(data = df.match1, 
                gcol = which(colnames(df.match1) == 'group'), 
                vcol = which(colnames(df.match1) == 'ebl.ratio'))

model.match1 <- glm(infect ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1)) # 95% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1)))
ci
or <- data.frame(exp(summary(model.match1)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.csv')

# section 2: IPTW
df.full <- subset(df, select = - c(trans.gap,trans.other,trans.leuko,trans.age,trans.timing))
head (df.full)
ps <- glm(group ~ age + sex + BMI + asa.clas + emer. + anes.gener + surg.dura + smoke 
          + Af + COPD + CI + SAH + radio + immu + chemo + HTN + DM + CAD + anem. + coagu 
          + liver + renal + open + lower + ebl.ratio + cancer + alb, 
          family = binomial('logit'), 
          data = df.full)
summary(ps)
df.full$score.ps <- predict(ps)
head(df.full)
hist(df.full$score.ps)

df.full$score.ps.IPTW <- predict(ps, type = "response") # calculate score.ps.IPTW
hist(df.full$score.ps.IPTW)
table(df.full$group)
head(df.full$score.ps.IPTW)
quantile(df.full$score.ps.IPTW)
quantile(df.full$score.ps.IPTW, prob = seq(0, 1, length = 11))
df.full$score.ps.IPTW.cate <- cut(df.full$score.ps.IPTW, c(quantile(df.full$score.ps.IPTW, prob = seq(0, 1, length = 11))[2], 
                                                           quantile(df.full$score.ps.IPTW, prob = seq(0, 1, length = 11))[10], 
                                                           quantile(df.full$score.ps.IPTW, prob = seq(0, 1, length = 11))[11]))
table(df.full$score.ps.IPTW.cate)
levels(df.full$score.ps.IPTW.cate)
levels(df.full$score.ps.IPTW.cate) <- c('1', '2')
df.full$score.ps.IPTW.cate <- as.character(df.full$score.ps.IPTW.cate)
df.full$score.ps.IPTW.cate <- as.numeric(df.full$score.ps.IPTW.cate)
df.full.quantile10 <- subset(df.full, df.full$score.ps.IPTW.cate == 1)
table(df.full.quantile10$infect,df.full.quantile10$group)
table(df.full.quantile10$group)
##  percentage (group=1)=0.4932
df.full.quantile10$iptw[df.full.quantile10$group == 1] <- 0.4932 / df.full.quantile10$score.ps.IPTW[df.full.quantile10$group == 1]
df.full.quantile10$iptw[df.full.quantile10$group == 0] <- (1 - 0.4932) / (1 - df.full.quantile10$score.ps.IPTW[df.full.quantile10$group == 0])
head(df.full.quantile10$iptw)
hist(df.full.quantile10$iptw)
quantile(df.full.quantile10$iptw,prob=seq(0,1,length=11))
max(df.full.quantile10$iptw)
head(df.full.quantile10)

model.ps.iptw <- glm(infect ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw)) # 95% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.ps.iptw.csv')

########################################################################################
# secondary outcomes

#SSI
table(df$SSI,df$group)
round(prop.table(table(df$SSI,df$group), 2), 3)

model4 <- glm(SSI ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.995)) # 99.5% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.995)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.SSI.csv')

model.ps.iptw <- glm(SSI ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.SSI.ps.iptw.csv')

table(df.match1$SSI,df.match1$group)
round(prop.table(table(df.match1$SSI,df.match1$group), 2), 3)
model.match1 <- glm(SSI ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1,level=0.995)))
ci
or <- data.frame(exp(summary(model.match1)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

#PNEU
table(df$PNEU,df$group)
round(prop.table(table(df$PNEU,df$group), 2), 3)

model4 <- glm(PNEU ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.995)) # 99.5% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.995)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.PNEU.csv')

model.ps.iptw <- glm(PNEU ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.PNEU.ps.iptw.csv')

table(df.match1$PNEU,df.match1$group)
round(prop.table(table(df.match1$PNEU,df.match1$group), 2), 3)

model.match1 <- glm(PNEU ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1)))
ci
or <- data.frame(exp(summary(model.match1,level=0.995)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

#BSI
table(df$BSI,df$group)
round(prop.table(table(df$BSI,df$group), 2), 3)

model4 <- glm(BSI ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.995)) # 99.5% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.995)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.BSI.csv')

model.ps.iptw <- glm(BSI ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.BSI.ps.iptw.csv')

table(df.match1$BSI,df.match1$group)
round(prop.table(table(df.match1$BSI,df.match1$group), 2), 3)

model.match1 <- glm(BSI ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1)))
ci
or <- data.frame(exp(summary(model.match1,level=0.995)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

#UTI
table(df$UTI,df$group)
round(prop.table(table(df$UTI,df$group), 2), 3)

model4 <- glm(UTI ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.995)) # 99.5% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.995)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.UTI.csv')

model.ps.iptw <- glm(UTI ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.UTI.ps.iptw.csv')

table(df.match1$UTI,df.match1$group)
round(prop.table(table(df.match1$UTI,df.match1$group), 2), 3)

model.match1 <- glm(UTI ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level = 0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1)))
ci
or <- data.frame(exp(summary(model.match1 = 0.995)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

#GI
table(df$GI,df$group)
round(prop.table(table(df$GI,df$group), 2), 3)

model4 <- glm(GI ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI +  COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.995)) # 99.5% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.995)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.GI.csv')

model.ps.iptw <- glm(GI ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.GI.ps.iptw.csv')

table(df.match1$GI,df.match1$group)
round(prop.table(table(df.match1$GI,df.match1$group), 2), 3)

model.match1 <- glm(GI ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1)))
ci
or <- data.frame(exp(summary(model.match1,level=0.995)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

#overlap
table(df$overlap,df$group)
round(prop.table(table(df$overlap,df$group), 2), 3)

model4 <- glm(overlap ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.995)) # 99.5% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.995)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.overlap.csv')

model.ps.iptw <- glm(overlap ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.overlap.ps.iptw.csv')

table(df.match1$overlap,df.match1$group)
round(prop.table(table(df.match1$overlap,df.match1$group), 2), 3)
model.match1 <- glm(overlap ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1)))
ci
or <- data.frame(exp(summary(model.match1,level=0.995)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

#ICU
table(df$icu,df$group)
round(prop.table(table(df$icu,df$group), 2), 3)

model4 <- glm(icu ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.995)) # 99.5% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.995)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.icu.csv')

model.ps.iptw <- glm(icu ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.icu.ps.iptw.csv')

table(df.match1$icu,df.match1$group)
round(prop.table(table(df.match1$BSI,df.match1$group), 2), 3)

model.match1 <- glm(icu ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1,level=0.995)))
ci
or <- data.frame(exp(summary(model.match1)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

#reopen
table(df$reopen,df$group)
round(prop.table(table(df$reopen,df$group), 2), 3)

model4 <- glm(reopen ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.995)) # 99.5% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.995)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.reopen.csv')

model.ps.iptw <- glm(reopen ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.reopen.ps.iptw.csv')

table(df.match1$reopen,df.match1$group)
round(prop.table(table(df.match1$reopen,df.match1$group), 2), 3)

model.match1 <- glm(reopen ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1)))
ci
or <- data.frame(exp(summary(model.match1,level=0.995)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

#post.length
# median(df.secondary.exposures$post.length)
# quantile(df.secondary.exposures$post.length)
df.control<-subset(df, df$group==0)
head(df.control)
median(df.control$post.length)
quantile(df.control$post.length)

quantile(df$post.length)
df$post.length.cate <- cut(df$post.length, c(quantile(df$post.length)[1], quantile(df$post.length)[2], quantile(df$post.length)[3], quantile(df$post.length)[4], quantile(df$post.length)[5]))
table(df$post.length)
levels(df$post.length.cate)
levels(df$post.length.cate) <- c('1', '2', '3', '4')
df$post.length.cate<- as.character(df$post.length.cate)
df$post.length.cate<- as.numeric(df$post.length.cate)
f.post.length.cate<-factor(df$post.length.cate)

model4<-polr(f.post.length.cate~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
             + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
             + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
             data = df, 
             Hess = T)
summary(model4)
coef<-coef(summary(model4))
coef
p <- pnorm(abs(coef[, "t value"]), 
           lower.tail = FALSE) * 2
p
ci<-exp(confint(model4,level=0.995))
ci
or<-exp(coef(model4))
or
brant(model4) #test for parallel assumption
p <-data.frame(p) [-27,]
p <-data.frame(p) [-27,]
p <-data.frame(p) [-27,] 
p <- data.frame(p)
p
ci <- data.frame(ci)
ci
or <- data.frame(or)
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.post.length.csv')

quantile(df.full.quantile10$post.length)
df.full.quantile10$post.length.cate <- cut(df.full.quantile10$post.length, c(quantile(df.full.quantile10$post.length)[1], quantile(df.full.quantile10$post.length)[2], quantile(df.full.quantile10$post.length)[3], quantile(df.full.quantile10$post.length)[4], quantile(df.full.quantile10$post.length)[5]))
levels(df.full.quantile10$post.length.cate)
levels(df.full.quantile10$post.length.cate) <- c('1', '2', '3', '4')
df.full.quantile10$post.length.cate<- as.character(df.full.quantile10$post.length.cate)
df.full.quantile10$post.length.cate<- as.numeric(df.full.quantile10$post.length.cate)
f.post.length.cate<-factor(df.full.quantile10$post.length.cate)
model4<-polr(f.post.length.cate~ group, 
             data = df.full.quantile10,
             Hess = T,
             weights = iptw)
summary(model4)
coef<-coef(summary(model4))
coef
p <- pnorm(abs(coef[, "t value"]), 
           lower.tail = FALSE) * 2
p
ci<-exp(confint(model4,level=0.995))
ci
or<-exp(coef(model4))
or
brant(model4) #test for parallel assumption
p <-data.frame(p) [-27,]
p <-data.frame(p) [-27,]
p <-data.frame(p) [-27,] 
p <- data.frame(p)
p
ci <- data.frame(ci)
ci
or <- data.frame(or)
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.post.length.iptw.csv')

quantile(df.match1$post.length)
df.match1$post.length.cate <- cut(df.match1$post.length, c(quantile(df.match1$post.length)[1], quantile(df.match1$post.length)[2], quantile(df.match1$post.length)[3], quantile(df.match1$post.length)[4], quantile(df.match1$post.length)[5]))
levels(df.match1$post.length.cate)
levels(df.match1$post.length.cate) <- c('1', '2', '3', '4')
df.match1$post.length.cate <- as.character(df.match1$post.length.cate)
df.match1$post.length.cate <- as.numeric(df.match1$post.length.cate)
f.post.length.cate<-factor(df.match1$post.length.cate)

model4<-polr(f.post.length.cate~ group, 
             data = df.match1,
             Hess = T)
summary(model4)
coef<-coef(summary(model4))
coef
p <- pnorm(abs(coef[, "t value"]), 
           lower.tail = FALSE) * 2
p
ci<-exp(confint(model4,level=0.995))
ci
or<-exp(coef(model4))
or
brant(model4) #test for parallel assumption

#death
table(df.full$death,df.full$group)
round(prop.table(table(df.full$death,df.full$group), 2), 3)
table(df.full.quantile10$death,df.full.quantile10$group)
round(prop.table(table(df.full.quantile10$death,df.full.quantile10$group), 2), 3)
table (df.match1$death, df.match1$group)
round(prop.table(table(df.match1$death,df.match1$group), 2), 3)

model.ps.iptw <- glm(death ~ group, 
                     family = quasibinomial('logit'), 
                     data = df.full.quantile10,
                     weights = iptw)
summary(model.ps.iptw) # P values 
exp(summary(model.ps.iptw)[['coefficients']][, 1]) # ORs
exp(confint(model.ps.iptw,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.ps.iptw)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.ps.iptw,level=0.995)))
ci
or <- data.frame(exp(summary(model.ps.iptw)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.death.ps.iptw.csv')

model.match1 <- glm(death ~ group, 
                    family = binomial('logit'), 
                    data = df.match1)
summary(model.match1) # P values 
exp(summary(model.match1)[['coefficients']][, 1]) # ORs
exp(confint(model.match1,level=0.995)) # 99.5% CIs of ORs
p <- summary(model.match1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model.match1,level=0.995)))
ci
or <- data.frame(exp(summary(model.match1)[['coefficients']][, 1]))
or
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results.match1.1.csv')

table(df.match1$death, df.match1$group)
rma.peto(ai = as.matrix(table(df.match1$death, df.match1$group))[4], 
         n1i = length(df.match1$id[df.match1$group == 1]), 
         ci = as.matrix(table(df.match1$death, df.match1$group))[2], 
         n2i = length(df.match1$id[df.match1$group == 0]),
         level = 99.5)

###########################################################################################
# dose-response
df.secondary.exposures<-subset(df, df$group == 1)
quantile(df.secondary.exposures$trans.unit)
df$trans.unit.cate <- cut(df$trans.unit, c(0, 0.1, 2, 4, 6, 8, 142),include.lowest =T)                          
table(df$trans.unit.cate)
head(df$trans.unit.cate)
levels(df$trans.unit.cate)
levels(df$trans.unit.cate) <- c('0','1', '2', '3', '4','5')
df$trans.unit.cate <- as.character(df$trans.unit.cate)
df$trans.unit.cate <- as.numeric(df$trans.unit.cate)
table(df$infect,df$trans.unit.cate)
round(prop.table(table(df$infect,df$trans.unit.cate), 2), 3)

model1 <- glm(infect ~ factor(trans.unit.cate), 
              family = binomial('logit'), 
              data = df)
summary(model1) # P values 
exp(summary(model1)[['coefficients']][, 1]) # ORs
exp(confint(model1)) # 95% CIs of ORs
p <- summary(model1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model1)))
ci
or <- data.frame(exp(summary(model1)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results1.dose.csv')

model2 <- glm(infect ~ factor(trans.unit.cate)+ age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu, 
              family = binomial('logit'), 
              data = df)
summary(model2) # P values 
exp(summary(model2)[['coefficients']][, 1]) # ORs
exp(confint(model2)) # 95% CIs of ORs
p <- summary(model2)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model2)))
ci
or <- data.frame(exp(summary(model2)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results2.dose.csv')

model3 <- glm(infect ~factor(trans.unit.cate) + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,              
              family = binomial('logit'), 
              data = df)
summary(model3) # P values 
exp(summary(model3)[['coefficients']][, 1]) # ORs
exp(confint(model3)) # 95% CIs of ORs
p <- summary(model3)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model3)))
ci
or <- data.frame(exp(summary(model3)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results3.dose.csv')

model4 <- glm(infect ~ factor(trans.unit.cate) +age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.dose.csv')

# p for trend
df$trans.unit.median[df$trans.unit.cat==0] <- median(df$trans.unit[df$trans.unit.cat==0])
df$trans.unit.median[df$trans.unit.cat==1] <- median(df$trans.unit[df$trans.unit.cat==1])
df$trans.unit.median[df$trans.unit.cat==2] <- median(df$trans.unit[df$trans.unit.cat==2])
df$trans.unit.median[df$trans.unit.cat==3] <- median(df$trans.unit[df$trans.unit.cat==3])
df$trans.unit.median[df$trans.unit.cat==4] <- median(df$trans.unit[df$trans.unit.cat==4])
df$trans.unit.median[df$trans.unit.cat==5] <- median(df$trans.unit[df$trans.unit.cat==5])
table(df$trans.unit.median)
df$trans.unit.median <- as.character(df$trans.unit.median)
df$trans.unit.median <- as.numeric(df$trans.unit.median)

model1 <- glm(infect ~ trans.unit.median, 
              family = binomial('logit'), 
              data = df)
summary(model1) # P values 
exp(summary(model1)[['coefficients']][, 1]) # ORs
exp(confint(model1)) # 95% CIs of ORs
p <- summary(model1)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model1)))
ci
or <- data.frame(exp(summary(model1)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results1.median.csv')

model2 <- glm(infect ~ trans.unit.median+ age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu, 
              family = binomial('logit'), 
              data = df)
summary(model2) # P values 
exp(summary(model2)[['coefficients']][, 1]) # ORs
exp(confint(model2)) # 95% CIs of ORs
p <- summary(model2)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model2)))
ci
or <- data.frame(exp(summary(model2)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results2.median.csv')

model3 <- glm(infect ~trans.unit.median + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,              
              family = binomial('logit'), 
              data = df)
summary(model3) # P values 
exp(summary(model3)[['coefficients']][, 1]) # ORs
exp(confint(model3)) # 95% CIs of ORs
p <- summary(model3)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model3)))
ci
or <- data.frame(exp(summary(model3)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results3.median.csv')

model4 <- glm(infect ~ trans.unit.median + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.trans.unit.cate.median.csv')

##############################################################################################
#subgroup analysis
df.depart1 <- subset(df, df$depart == 1)
head(df.depart1)
length(df.depart1$id)
length(unique(df.depart1$id))
sum(df.depart1$infect)

table(df.depart1$group)
round(prop.table(table(df.depart1$group)), 3)
table(df.depart1$group, df.depart1$infect)
round(prop.table(table(df.depart1$group, df.depart1$infect), 2), 3)

model4 <- glm(infect ~ group +age + sex + BMI + asa.clas + HTN + DM + CAD + CI + COPD
              + smoke + anem. + coagu + liver + renal + cancer + radio
              + surg.dura + open, 
              family = binomial('logit'), 
              data = df.depart1)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.983)) # 98.3% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.983)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.depart1.csv')

#depart2
df.depart2 <- subset(df, df$depart == 2)
head(df.depart2)
length(df.depart2$id)
length(unique(df.depart2$id))
sum(df.depart2$infect)

table(df.depart2$group)
round(prop.table(table(df.depart2$group)), 3)
table(df.depart2$infect, df.depart2$group)
round(prop.table(table(df.depart2$infect, df.depart2$group), 2), 3)

table(df.depart2$infect, df.depart2$anes.gener)
table(df.depart2$infect, df.depart2$cancer)
table(df.depart2$infect, df.depart2$SAH)
table(df.depart2$infect, df.depart2$Af)
table(df.depart2$infect, df.depart2$radio)
table(df.depart2$infect, df.depart2$emer.)
table(df.depart2$infect, df.depart2$lower)

model4 <- glm(infect ~ group + age + sex + BMI + asa.clas + HTN + DM + CAD + CI + COPD + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + surg.dura + open + ebl.ratio, 
              family = binomial('logit'), 
              data = df.depart2)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 2]) # ORs
exp(confint(model4,level=0.983)) # 98.3% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.983)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 2]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.depart2.csv')

#depart34
df.depart3 <- subset(df, df$depart == 3)
df.depart4 <- subset(df, df$depart == 4)
df.depart34 <- rbind.data.frame(df.depart3,df.depart4)
head(df.depart34)
length(df.depart34$id)
length(unique(df.depart34$id))
sum(df.depart34$infect)

table(df.depart34$group)
round(prop.table(table(df.depart34$group)), 3)
table(df.depart34$infect, df.depart34$group)
round(prop.table(table(df.depart34$infect, df.depart34$group), 2), 3)

model4 <- glm(infect ~ group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH + COPD + smoke + immu + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df.depart34)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.983)) # 98.3% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.983)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.depart34.csv')

#depart56
df.depart5 <- subset(df, df$depart == 5)
df.depart6 <- subset(df, df$depart == 6)
df.depart56 <- rbind.data.frame(df.depart5,df.depart6)
head(df.depart56)
length(df.depart56$id)
length(unique(df.depart56$id))
sum(df.depart56$infect)

table(df.depart56$group)
round(prop.table(table(df.depart56$group)), 3)
table(df.depart56$infect, df.depart56$group)
round(prop.table(table(df.depart56$infect, df.depart56$group), 2), 3)

model4 <- glm(infect ~ group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH + COPD + smoke + immu + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio, 
              family = binomial('logit'), 
              data = df.depart56)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4,level=0.983)) # 98.3% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4,level=0.983)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.depart56.csv')

# Sensitivity analysis
# pod 30
df.post30 <- subset(df, df$gap.post.30 == 1)
table(df.post30$group)
round(prop.table(table(df.post30$group)), 3)
table(df.post30$infect, df.post30 $group)
round(prop.table(table(df.post30$infect, df.post30$group), 2), 3)

model4 <- glm(infect ~group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,
              family = binomial('logit'), 
              data = df.post30)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.post30.csv')

# pod 14 days
df.post14<- subset(df, df$gap.post.14 == 1)
table(df.post14$group)
round(prop.table(table(df.post14$group)), 3)
table(df.post14$infect, df.post14$group)
round(prop.table(table(df.post14$infect, df.post30$group), 2), 3)

model4 <- glm(infect ~group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,
              family = binomial('logit'), 
              data = df.post14)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.post14.csv')

# pod 7 days
df.post7 <- subset(df, df$gap.post.7 == 1)
table(df.post7$group)
round(prop.table(table(df.post7$group)), 3)
table(df.post7$infect, df.post7$group)
round(prop.table(table(df.post7$infect, df.post7$group), 2), 3)

model4 <- glm(infect ~group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,
              family = binomial('logit'), 
              data = df.post7)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.post7.csv')

# pod 3 days
df.post3 <- subset(df, df$gap.post.3 == 1)
table(df.post3$group)
round(prop.table(table(df.post3$group)), 3)
table(df.post3$infect, df.post3$group)
round(prop.table(table(df.post3$infect, df.post3$group), 2), 3)

model4 <- glm(infect ~group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,
              family = binomial('logit'), 
              data = df.post3)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.post3.csv')

# intra
df.intra<- subset(df, df$trans.timing == 0)
table(df.intra$group)
round(prop.table(table(df.intra$group)), 3)
table(df.intra$infect, df.intra$group)
round(prop.table(table(df.intra$infect, df.intra$group), 2), 3)

model4 <- glm(infect ~group + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,
              family = binomial('logit'), 
              data = df.intra)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.intra.csv')

# pre vs. intro vs. post
df.pre <- subset(df, df$trans.timing == 1)
table(df.pre$anem., df.pre$group)
round(prop.table(table(df.pre$anem., df.pre$group), 2), 3)
df$trans.timing [df$group == 0] <- 4
df$trans.timing [df$trans.timing == 0] <- 3 #intraoperative=3, pre-=1, post-=2
df$trans.timing [df$group == 0] <- 0 # control=0
table (df$infect, df$trans.timing)
round(prop.table(table(df$infect, df$trans.timing), 2), 3)

model4 <- glm(infect ~factor(trans.timing) + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
presults <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.trans.timing.csv')

# leuco-depleted
df$trans.leuko [is.na(df$trans.leuko) == TRUE] <- 3 # control 3, leuko 1, non-leuko 2
df$trans.leuko [df$trans.leuko == 0] <- 4
df$trans.leuko [df$trans.leuko == 3] <- 0
df$trans.leuko [df$trans.leuko == 4] <- 2 # control 0, leuko 1, non-leuko 2
table (df$trans.leuko)
table (df$infect, df$trans.leuko)
round(prop.table(table(df$infect, df$trans.leuko), 2), 3)

model4 <- glm (infect ~factor(trans.leuko) + age + sex + BMI + asa.clas + HTN + DM + CAD + Af + CI + SAH+ COPD
              + smoke + anem. + coagu + liver + renal + alb+ cancer + chemo + radio + immu
              + anes.gener + surg.dura + open + lower + emer.+ ebl.ratio,
              family = binomial('logit'), 
              data = df)
summary(model4) # P values 
exp(summary(model4)[['coefficients']][, 1]) # ORs
exp(confint(model4)) # 95% CIs of ORs
p <- summary(model4)[['coefficients']][, 4]
p <- data.frame(p)
p
ci <- data.frame(exp(confint(model4)))
ci
or <- data.frame(exp(summary(model4)[['coefficients']][, 1]))
or 
results <- data.frame(c(data.frame(rownames(or)), round(or, 3), round(ci, 3), round(p, 3)))
write.csv(results, 'results4.trans.leuko.csv')

