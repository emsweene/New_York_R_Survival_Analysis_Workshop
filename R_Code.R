###############################################################
## Code chunk 1: investigating the ovarian dataset  
###############################################################

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(xtable) 
library(JM)
library(cmprsk)


##load the ovarian dataset from the survival package 
data(ovarian)

##check the help for the ovarian dataset 
?ovarian

##look at the data 
View(ovarian)

##arrange by rx 
##add an index variable
##rename the rx and fustat variables 
ovarian <- ovarian %>% 
  arrange(rx) %>% 
  mutate(index = 1:n())


##plot the data 
ggplot(ovarian, 
       aes(xend = 0, 
           y = index, 
           x = futime, 
           yend = index, 
           colour = factor(rx),
           shape = factor(fustat))) + 
  geom_segment() + 
  geom_point() 


###############################################################
## Code chunk 2: survival function for exponential data 
###############################################################

##set the values of x 
exponential.data <- data.frame( x = 1:100) 

## set the expected lifetime (1/\lambda)
lambda <- 1/10

##calculate the pdf, cdf, and survival function  
exponential.data <- 
  exponential.data %>% 
  mutate(f = lambda * exp(- lambda * x),
         F = 1 - exp(- lambda * x), 
         S = exp(- lambda * x)) 

##plot the pdf 
ggplot(exponential.data, aes(x = x, y = f)) + 
  geom_line() 

##plot the cdf 
ggplot(exponential.data, aes(x = x, y = F)) + 
  geom_line() 


##plot the survival function  
ggplot(exponential.data, aes(x = x, y = S)) + 
  geom_line() 



###############################################################
## In-class exercise 1 solution 
###############################################################

##set the values of x 
weibull.data <- data.frame( x = 1:100) 

## set the parameters 
lambda <- 1/20
alpha <- 1.5


##calculate the pdf, cdf, and survival function  
weibull.data <- 
  weibull.data %>% 
  mutate(f = alpha * lambda * x ^(alpha - 1) * exp(- lambda * x ^ alpha),
         F = 1 - exp(- lambda * x ^ alpha), 
         S = exp(- lambda * x^alpha))  

##plot the pdf 
ggplot(weibull.data, aes(x = x, y = f)) + 
  geom_line() 

##plot the cdf 
ggplot(weibull.data, aes(x = x, y = F)) + 
  geom_line() 

##plot the survival function  
ggplot(weibull.data, aes(x = x, y = S)) + 
  geom_line() 


###############################################################
## Code chunk 3: hazard function for exponential data 
###############################################################

##calculate the hazard 
exponential.data <- exponential.data %>% 
  mutate(h = f/S)

##plot the hazard 
ggplot(exponential.data, aes(x = x, y = h)) + 
  geom_line() 


###############################################################
## In-class exercise 2 solution 
###############################################################

##calculate the hazard 
weibull.data <- weibull.data %>% 
  mutate(h = f/S)

##plot the hazard 
ggplot(weibull.data, aes(x = x, y = h)) + 
  geom_line() 


###############################################################
## Code chunk 4: short lubridate tutorial 
## Reference: https://data.library.virginia.edu/working-with-dates-and-time-in-r-using-the-lubridate-package/
###############################################################

library(lubridate)

##type in dates or beginning an end in two different formats 
begin <- c("May 11, 1996", "September 12, 2001", "July 1, 1988")
end <- c("7/8/97","10/23/02","1/4/91")

class(begin)
class(end)

?mdy 

##parse the dates with year, month, and day components 
begin <- mdy(begin)
end <- mdy(end)

class(begin)
class(end)

## The dates now have class “Date” and are printed in year-month-day format. 
## They may appear to still be character data when printed, but they are in 
## fact numbers. The “Date” class means dates are stored as the number of 
## days since January 1, 1970, with negative values for earlier dates. We can 
## use the as.numeric function to view the raw values.

as.numeric(begin)

as.numeric(end)

## We can now subtract the dates! 

diff = end - begin 




###############################################################
## Code chunk 5: the jasa dataset 
###############################################################


jasa <- jasa %>% 
  arrange(transplant) %>% 
  mutate(index = 1:n())

ggplot(jasa, 
       aes(xend = 0, 
           y = index, 
           x = futime, 
           yend = index, 
           colour = factor(transplant),
           shape = factor(fustat))) + 
  geom_segment() + 
  geom_point() 


###############################################################
## In class exercise 3 -- the veterans dataset 
###############################################################


data(veteran)
?veteran

veteran <- veteran %>% 
  arrange(trt) %>% 
  mutate(index = 1:n())

##plot the data 
ggplot(veteran, 
       aes(xend = 0, 
           y = index, 
           x = time, 
           yend = index, 
           colour = factor(trt),
           shape = factor(status))) + 
  geom_segment() + 
  geom_point() 



###############################################################
## Code chunk 6: survival curve for the psychiatric dataset  
###############################################################

psych.survival <- data.frame(t = c(0, 1,2,11,14,22,24,25,26,28,32,35,39, 40),
                             S = c(1, .923, .885, .847, .809, .732, .693, .654, .616, .578, .525, .42, .21, 0))


ggplot(psych.survival, aes(x = t, y = S)) + 
  geom_step() + 
  geom_hline(yintercept = .5, col = 'red')


###############################################################
## Code chunk 6: logrank test for the ovarian dataset   
###############################################################

surv_object_ovarian <- Surv(time = ovarian$futime, event = ovarian$fustat)

##note that this is a chi-square distribution with 1 degree of freedom 
##which is a standard normal; for more groups this will be a higher 
##degree of freedom 

survdiff(surv_object_ovarian ~ rx, data = ovarian)


###############################################################
## In class exercise 4 -- logrank test for the veterans dataset 
## and the jasa dataset 
###############################################################

surv_object_veteran <- Surv(time = veteran$time, 
                            event = veteran$status)

survdiff(surv_object_veteran ~ trt, data = veteran)

##we can also look at the logrank test for the different celltypes 
survdiff(surv_object_veteran ~ celltype, data = veteran)



surv_object_jasa <- Surv(time = jasa$futime,
                         event = jasa$fustat)

survdiff(surv_object_jasa ~ transplant, data = jasa)


###############################################################
## Code chunk 7: Kaplan Meier curves    
###############################################################

fit_ovarian_1 <- survfit(surv_object_ovarian ~ 1, data = ovarian)
 
ggsurvplot(fit_ovarian_1, 
           data = ovarian, 
           pval = FALSE,
           conf.int = TRUE)

fit_ovarian_rx <-  survfit(surv_object_ovarian ~ rx, data = ovarian)

ggsurvplot(fit_ovarian_rx , 
           data = ovarian, 
           pval = TRUE,
           conf.int = TRUE)



fit_veteran_1 <- survfit(surv_object_veteran ~ 1, data = veteran)

ggsurvplot(fit_veteran_1, 
           data = veteran, 
           pval = TRUE,
           conf.int = TRUE)


fit_veteran_trt <- survfit(surv_object_veteran ~ trt, data = veteran)

ggsurvplot(fit_veteran_trt, 
           data = veteran, 
           pval = TRUE,
           conf.int = TRUE)




fit_veteran_celltype <- survfit(surv_object_veteran ~ celltype, data = veteran)

ggsurvplot(fit_veteran_celltype, 
           data = veteran, 
           pval = TRUE,
           conf.int = TRUE)

###############################################################
## In class exercise 5 -- KM ciurves for the jasa dataset  
###############################################################


fit_jasa_1 <- survfit(surv_object_jasa ~ 1, data = jasa)

ggsurvplot(fit_jasa_1, 
           data = jasa, 
           pval = FALSE,
           conf.int = TRUE)


fit_jasa_transplant <- survfit(surv_object_jasa ~ transplant, data = jasa)

ggsurvplot(fit_jasa_transplant, 
           data = jasa, 
           pval = TRUE,
           conf.int = TRUE)


jasa <- jasa %>% 
  mutate(age.50 = age >= 50)

fit_jasa_age_50 <- survfit(surv_object_jasa ~ age.50, data = jasa)

ggsurvplot(fit_jasa_age_50, 
           data = jasa, 
           pval = TRUE,
           conf.int = TRUE)




###############################################################
## Code chunk 8: Cox proportional hazards model     
###############################################################

?ovarian 


ovarian$rx <- factor(ovarian$rx, 
                     levels = c("1", "2"), 
                     labels = c("A", "B"))
ovarian$resid.ds <- factor(ovarian$resid.ds, 
                           levels = c("1", "2"), 
                           labels = c("no", "yes"))
ovarian$ecog.ps <- factor(ovarian$ecog.ps, 
                          levels = c("1", "2"), 
                          labels = c("good", "bad"))


ovarian <- ovarian %>% mutate(age_group = ifelse(age >=50, "old", "young"))
ovarian$age_group <- factor(ovarian$age_group)


fit.ovarian <- coxph(surv_object_ovarian ~ rx + resid.ds + age_group + ecog.ps, 
                   data = ovarian)


summary(fit.ovarian )

xtable(fit.ovarian)

###############################################################
## Code chunk 9: ggforest    
###############################################################


ggforest(fit.ovarian, data = ovarian)


###############################################################
## Code chunk 10: Cox proportional hazards model for veteran     
###############################################################

veteran$trt <- factor(veteran$trt,
                      levels = c("1", "2"),
                      labels = c("standard", "test"))

veteran$prior <- factor(veteran$prior,
                        levels = c("0", "10"),
                        labels = c("no", "yes"))

coxph_veteran <- coxph(surv_object_veteran ~ trt + 
        celltype + 
        karno + 
        diagtime + 
        age + 
        prior, data = veteran) 


ggforest(coxph_veteran, data = veteran)


###############################################################
## In-Class Exercise 6 -- fit a cox proportional hazards model 
## to the jasa dataset 
###############################################################

##model comparing treatment 
jasa$transplant <- factor(jasa$transplant,
                          levels = c("0", "1"),
                          labels = c("no", "yes")) 
  
jasa$surgery <- factor(jasa$surgery,
                       levels = c("0", "1"),
                       labels = c("no", "yes")) 

fit_jasa_treatment <- coxph(surv_object_jasa ~ 
                              transplant +
                              age.50 +
                              surgery, data = jasa)


ggforest(fit_jasa_treatment, data = jasa)


##model in the group that got transplants 

jasa_transplant <- jasa %>% 
  filter(transplant == 'yes')

surv_object_jasa_transplant  <- Surv(time = jasa_transplant$futime,
                                     event = jasa_transplant$fustat)
  
jasa_transplant$reject <- factor(jasa_transplant$reject,
                                 levels = c("0", "1"),
                                 labels = c("no", "yes")) 


fit_jasa_transplant_1  <- coxph(surv_object_jasa_transplant ~ 
                                  age.50 +
                                  wait.time + 
                                  mscore + 
                                  reject + 
                                  surgery, 
                                data = jasa_transplant)


ggforest(fit_jasa_transplant_1 , data = jasa_transplant)


###############################################################
## Code chunk 11: Testing the proportional hazards model using 
## the scaled Schoenfeld residuals 
###############################################################

##ovarian 
zph.test.ovarian <- cox.zph(fit.ovarian)

ggcoxzph(zph.test.ovarian )

##veteran 
zph.test.veteran <- cox.zph(coxph_veteran)

ggcoxzph(zph.test.veteran)

plot(cox.zph(coxph_veteran), var = 5)


##jasa models 
zph.test.jasa.1 <-  cox.zph(fit_jasa_treatment)

ggcoxzph(zph.test.jasa.1)




###############################################################
## Code chunk 12: Adding in a time-varying coefficient to the model 
###############################################################

veteran2 <- survSplit(surv_object_veteran ~ trt + 
                        celltype + 
                        karno + 
                        diagtime + 
                        age + 
                        prior, 
                      data = veteran,
                      cut = c(70), 
                      episode = "time_group")

coxph_veteran.2 <- coxph(surv_object_veteran ~ trt + 
                         celltype + 
                         karno:strata(time_group) + 
                         diagtime + 
                         age + 
                         prior, data = veteran2) 

zph.test.veteran.2 <- cox.zph(coxph_veteran.2)

###############################################################
## Code chunk 13: time varying covariates  
###############################################################

data("aids", package = "JM")

surv_object_aids <- Surv(aids$start, aids$stop, aids$event)

coxph_aids <- coxph(surv_object_aids ~ CD4, data = aids)

summary(coxph_aids)



###############################################################
## Code chunk 14: cause-specific survival 
###############################################################


?BMT

## Label diseases
BMT$dis <- factor(BMT$dis, 
                  levels = c(0,1), 
                  labels = c("ALL", "AML"))
## Label status
BMT$status <- factor(BMT$status, 
                     levels = c(0,1,2), 
                     labels = c("Censored","Mortality","Relapse"))



surv_object_BMT_relapse <- Surv(BMT$ftime, BMT$status == "Relapse")


fit_jsurv_object_BMT_relapse <- survfit(surv_object_BMT_relapse ~ dis, data = BMT)

ggsurvplot(fit_jsurv_object_BMT_relapse, 
           data = BMT, 
           pval = TRUE,
           conf.int = TRUE)


###############################################################
## Code chunk 15: CIF 
###############################################################


## Calculate the grouped cumulative incidence functions (CIF)
bmt_cum_incidence <- cuminc(ftime   = BMT$ftime,  # failure time variable
                            fstatus = BMT$status,  # variable with distinct codes for different causes of failure
                            group   = BMT$dis,  # estimates will calculated within groups
                            cencode = "Censored", # value of fstatus variable which indicates the failure time is censored
                            )

plot(bmt_cum_incidence)



###############################################################
## Code chunk 17: Competing risk regression 
###############################################################


bmt_dis_mat <- matrix(as.numeric(BMT$dis == "AML"))
colnames(bmt_dis_mat) <- "dis"

##use relapse as the outcome of interest 
bmt_crr_relapse <- crr(ftime = BMT$ftime, # vector of failure/censoring times
                      fstatus  = BMT$status, # vector with a unique code for each failure type and censoring
                      cov1     = bmt_dis_mat, #  matrix (nobs x ncovs) of fixed covariates
                      failcode = "Relapse", # code of fstatus that denotes the failure type of interest
                      cencode  = "Censored" # code of fstatus that denotes censored observations
)

summary(bmt_crr_relapse)


##use mortality as the outcome of interest 
bmt_crr_mortality <- crr(ftime = BMT$ftime, # vector of failure/censoring times
                       fstatus  = BMT$status, # vector with a unique code for each failure type and censoring
                       cov1     = bmt_dis_mat, #  matrix (nobs x ncovs) of fixed covariates
                       failcode = "Mortality", # code of fstatus that denotes the failure type of interest
                       cencode  = "Censored" # code of fstatus that denotes censored observations
)

summary(bmt_crr_mortality) 








