# Introduction to survival analysis workshop
# practical session with R
# CBB, 2024/03/19

# BEFORE THE START: install packages survival and surminer
# install.packages(c("survival", "survminer"))
install.packages("survival")
install.packages("survminer")

library("survival")
library("survminer")

# IMPORT OR LOAD DATA
# One can use the drop-down menu or the R command
# the same dataset can be loaded from a folder, local or from the web. A csv file
# can be found at
# lung <- read.csv("https://raw.githubusercontent.com/davval-ki/IntroSurv/main/lung_cancer.csv")

# We’ll use the lung cancer data available in the survival package.
data("cancer")
# how is the data structured? Variables, variable types, observations
? help
head(lung) 

#First, it is necessary to create a "Survival" object
?Surv
# It can do much more than this (i.e. time dependent variables, etc)
Surv(lung$time, lung$status)
# the output contains a surv object, with times and censoring status. Takes 0/1, T/F or 1/2
# one can also specify the event of interest...
Surv(lung$time, lung$status==1)

# We want to compute the survival curves: survfit()
? survfit
fit1 <- survfit(Surv(time, status==1) ~ 1, data = lung)
 print(fit1)

#if you want to display a more complete summary of the survival curves, type this:
# Summary of survival curves
summary(fit1)
# Access to the sort summary table
summary(fit1)$table

# Optional: One can access and extract different values from this object.
# to see which ones, ?survfit.object
# fit1$n.event

#what we can simply do now is plotting the KM curve
plot(fit1)

#But let's try a nicer system, using survminer
ggsurvplot(fit1)
? ggsurvplot


# now let's compute the survival separately for two groups, for example by sex at birth and and 
# plot the corresponding curves,
fit <- survfit(Surv(time, status) ~ sex, data = lung)
print(fit)
summary(fit)
summary(fit)$table

# and draw a even more informative plot 
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

# many possibilities, among the different parameters
# 1) xlim, to shorten your plot

  ggsurvplot(fit,
             conf.int = TRUE,
             risk.table.col = "strata", # Change risk table color by groups
             ggtheme = theme_bw(), # Change ggplot2 theme
             palette = c("#E7B800", "#2E9FDF"),
             xlim = c(0, 600))
  
# fun="event" , to plot the event function (relapses, ect)
  
  ggsurvplot(fit,
             conf.int = TRUE,
             risk.table.col = "strata", # Change risk table color by groups
             ggtheme = theme_bw(), # Change ggplot2 theme
             palette = c("#E7B800", "#2E9FDF"),
             fun = "event")
  

  # Log-Rank test comparing survival curves: survdiff()
  surv_diff <- survdiff(Surv(time, status) ~ sex, data = lung)
  surv_diff
  
# What about age? we can compare survival by different age groups..
  lung$agegrp <- NA
  lung$agegrp[lung$age < 55] <- "less than 55"
  lung$agegrp[lung$age >= 55 & lung$age < 70] <- "55-69"
  lung$agegrp[lung$age >= 70] <- "70 and above"
  
  
  # do it yourself. Try to compute survival and draw a KM curve by age groups
  fitAge <- survfit(Surv(time, status) ~ agegrp, data = lung)
  ggsurvplot(fitAge,
             conf.int = T,
             risk.table.col = "strata", # Change risk table color by groups
             ggtheme = theme_bw())
  # Now, try to "turn off" the C.I
  
    
 ##### R function to compute the Cox model: coxph()
  
  # 1) Compute the Cox model
  # We’ll fit the Cox regression using the following covariates: 
  # age, sex, ph.ecog and wt.loss.
  
  # We start by computing univariate Cox analyses for all these variables; 
  # then we’ll fit multivariate cox analyses using two variables 
  # to describe how the factors jointly impact on survival.
  
  # Univariate Cox regression
  res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
  res.cox
  
  #again for a more complete report
  summary(res.cox)
  
  # Exercise: do the same (univariate) for variables 
  #age, ph.ecog, wt.loss
  # the result will be
  #          beta HR (95% CI for HR) wald.test p.value
  #age      0.019            1 (1-1)       4.1   0.042
  #sex      -0.53   0.59 (0.42-0.82)        10  0.0015
  #ph.ecog   0.48        1.6 (1.3-2)        18 2.7e-05
  #wt.loss 0.0013         1 (0.99-1)      0.05    0.83
  
  # Multivariate Cox regression analysis
  
  res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)
  summary(res.cox)
  
  #Visualizing the estimated distribution of survival times
  # Plot the baseline survival function
  ggsurvplot(survfit(res.cox), data=lung, palette = "#2E9FDF",
             ggtheme = theme_minimal())
  
  # Plot the hazard ratios
  ggforest(res.cox, data=lung)
  
  
#######  C)	Testing Cox assumptions
#Testing the proportional hazards assumption: function cox.zph

  test.ph <- cox.zph(res.cox)
  test.ph

  #From the output above, the test is not statistically significant for each of the covariates, and the global test is also not statistically significant. 
  #Therefore, we can assume the proportional hazards.
  
  
#  Plot the scaled Schoenfeld residuals against the transformed time 
  dev.new()
  plot(test.ph[1], lwd=1.5, col="red")
  abline(h=0, col="blue")
  plot(test.ph[2], lwd=1.5, col="red")
  abline(h=0, col="blue")
  plot(test.ph[3], lwd=1.5, col="red")
  abline(h=0, col="blue")
  # Schoenfeld residuals "can essentially be thought of as the observed minus 
  #the expected values of the covariates at each failure time" 
  #The plot of Schoenfeld residuals against time for any covariate should 
  #not show a pattern of changing residuals for that covariate.
  
  #Schoenfeld residuals are the differences between observed
  #and expected values of a covariate, adjusted for the rest of 
  #the model. The plot gives The smoothed plot is thus an estimate 
  #of the time dependence of the coefficient for the covariate
  #If the Schoenfeld residuals show a pattern over time 
  #(e.g., if they exhibit a non-linear trend), it suggests that the proportional 
  #hazards assumption might be violated for that covariate.
  #if there's a noticeable trend, it indicates that the hazard ratio
  #for that covariate is changing over time. 
  #This could mean that the effect of the covariate on survival is not constant over time, 
  
  # what to do if assumptions are violated? It depends.
  # -if for one variable only, one can try to stratify for that variable
  # -maybe there is a time dependent covariate?
  # -different model. or even, different methods (splines?)
  
  
  ggcoxdiagnostics(res.cox, type = "dfbeta",
                   linear.predictions = FALSE, ggtheme = theme_bw())
  
  
  