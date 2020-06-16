library(tableone)
library(MatchIt)
library(ggplot2)


#load dataset
load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))


# create new dataset: only variables that will be used, convert character to numeric

ARF <- as.numeric(rhc$cat1 == 'ARF')
CHF<- as.numeric(rhc$cat1=='CHF')
Cirr<- as.numeric(rhc$cat1 == 'Cirrhosis')
colcan <- as.numeric(rhc$cat1 == 'Colon Cancer')
Coma<- as.numeric(rhc$cat1 == 'Coma')
COPD <- as.numeric(rhc$cat1 == 'COPD')
lungcan <- as.numeric(rhc$cat1 == 'Lung Cancer')
MOSF <- as.numeric(rhc$cat1 == 'MOSF w/Malignancy')
sepsis<- as.numeric(rhc$cat1 == 'MOSF w/Sepsis')
female <- as.numeric(rhc$sex == 'Female')
died <- as.numeric(rhc$death == 'Yes')
age <- rhc$age
treatment <- as.numeric(rhc$swang1 == 'RHC')
meanbp1 <- rhc$meanbp1
aps <- rhc$aps1


mydata<- cbind(ARF, CHF, Cirr, colcan, Coma, lungcan, MOSF, sepsis, age, female, meanbp1, aps, treatment, died)
mydata <- data.frame(mydata)


#compute propensity score using logistic models
psmodel <- glm(treatment~ ARF+ CHF+ Cirr+ colcan + Coma+ lungcan+ MOSF + sepsis + age + female + meanbp1 + aps,
               family = binomial(), data = mydata)
#show coefficients
summary(psmodel)
# that's propensity score
pscore <- psmodel$fitted.values

#now look at the plot of propensity score, prematching
hist_control<- hist(pscore[mydata$treatment==0])
hist_treat<- hist(pscore[mydata$treatment==1])
plot(hist_control, col= 'red',, main = 'Plot of Propensity Socre, Pre-matching', xlab= 'score')
plot(hist_treat, add= TRUE, col='blue')
plot(hist_control, col=rgb(1,0,0,0.5), add=T)

#use matchit for propensity score, nearest neigbor matching

m.out <- matchit(treatment~ ARF+CHF+Cirr+ colcan + Coma + lungcan + MOSF + sepsis +
                   age+ female + meanbp1+ aps, data= mydata, method = 'nearest')

summary(m.out)

#propensity score plots
plot(m.out, type =  'jitter') #make sure the distribution looks similar for the matched units
plot(m.out, type = 'hist')

#use differnt package

#No caliper used here:
library(Matching)
psmatch<- Match(Tr= mydata$treatment, M=1, X= log(pscore), replace= FALSE)
matched <- mydata[unlist(psmatch[c('index.treated', 'index.control')]),]
xvars <- c('ARF', 'CHF', 'Cirr', 'colcan', 'Coma', 'lungcan', 'MOSF', 'sepsis',
           'age', 'female', 'meanbp1', 'aps')
matchedtab1<- CreateTableOne(vars= xvars, strata = 'treatment', data = matched, test= FALSE)
print(matchedtab1,smd= TRUE)

#use caliper
psmatch<- Match(Tr= mydata$treatment, M=1, X= log(pscore), replace= FALSE, caliper = 0.2)
matched <- mydata[unlist(psmatch[c('index.treated', 'index.control')]),]
xvars <- c('ARF', 'CHF', 'Cirr', 'colcan', 'Coma', 'lungcan', 'MOSF', 'sepsis',
           'age', 'female', 'meanbp1', 'aps')

matchedtab1<- CreateTableOne(vars= xvars, strata = 'treatment', data = matched, test= FALSE)
print(matchedtab1,smd= TRUE)







