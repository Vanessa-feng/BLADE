library(readxl)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(reshape2)

clinic_data <- read.csv("data/epoc_pathology_image_data/epoc_clinic.csv", header=T)

clinic_data$PFI_surv <- as.factor(clinic_data$PFI_surv)
surv_obj4 <- Surv(time = clinic_data$PFI_time, 
                  event = as.numeric(clinic_data$PFI_surv))
cox4 <- coxph(surv_obj4 ~ median_numlayer_output + cluster(case), data = clinic_data)
summary(cox4)


clinic_data <- clinic_data %>%
  mutate(
    gender = factor(gender),
    First.group = factor(First.group),
    HistBaselineR3 = factor(HistBaselineR3),
    Prior.oral.cancer = factor(Prior.oral.cancer.v2..1.yes..0.no.),
    Concurrent.oral.cancer = factor(Concurrent.oral.cancer.v2..1.yes..0.no.),
    Leukoplakia.group = factor(Leukoplakia.group.v2..1.Leukoplakia..0.Other.),
    LOH_Eligible = factor(LOH_Eligible..1.LOH..2.no.LOH.),
    TP53.mut_status = factor(TP53.mut_status)
  )

cox4_covar <- coxph(
  formula = surv_obj4 ~ median_numlayer_output + age + gender +
    Prior.oral.cancer + Concurrent.oral.cancer + Leukoplakia.group + LOH_Eligible + cluster(case), 
  data = clinic_data
)
summary(cox4_covar)
