


library(data.table)
library(stringr)
library(dplyr)
library(survival)
library(survminer)
library(tableone)
library(forestplot)
library(formattable)
library(ggplot2)

options(datatable.fread.datatable=FALSE)

setwd("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/")
PRS_int = fread("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/PRSCs_interaction.txt")
PRS_noint = fread("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/PRSCS_noInteraction.csv")
# PRS_int_BP = fread("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/gwPRS_PRSCS_SBP_DBP.csv")

PRS_int_sig = PRS_int %>% filter(p_interaction_cox<0.05/13)
PRS_int_sig = PRS_int_sig$Risk_factor
paste(PRS_int_sig, collapse="','")

# PRS_int = fread("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/PRS_indep_raw.csv")
# PRS_int_BP = fread("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/PRS_indep_raw_SBP_DBP.csv")

# PRS_int = bind_rows(PRS_int, PRS_int_BP)
rr = function(x, d) return(round(x,d))

head(PRS_int)

colnames(PRS_int)

PRS_int_1 = PRS_int %>% 
  filter(Risk_factor != "CKD_has_disease") %>%
  mutate(c_interaction_cox_out = paste(round(c_interaction_cox, 3), " (", round(c_interaction_cox_lowerCI, 3), "-", round(c_interaction_cox_upperCI,3), ")", sep="")) %>%
  mutate(c_riskfactorOnly_cox_out = paste(round(c_riskfactorOnly_cox, 3), " (", round(c_riskfactorOnly_cox_lowerCI, 3), "-", round(c_riskfactorOnly_cox_upperCI,3), ")", sep="")) %>%
  mutate(hr_rf = paste0(round(hr_riskfactorOnly_cox, 2), " [", round(hr_riskfactorOnly_cox_lowerCI, 2), "-", round(hr_riskfactorOnly_cox_upperCI, 2), "]")) %>%
  mutate(hr_int_rf = paste0(round(hr_riskfactor_interaction_cox, 2), " [", round(hr_riskfactor_interaction_cox_lowerCI, 2), "-", round(hr_riskfactor_interaction_cox_upperCI, 2), "]")) %>%
  mutate(hr_noint_PRS = paste0(round(hr_PRS_noninteraction_cox, 2), " [", round(hr_PRS_noninteraction_cox_lowerCI, 2), "-", round(hr_PRS_noninteraction_cox_upperCI, 2), "]")) %>%
  mutate(hr_int_PRS = paste0(round(hr_PRS_interaction_cox, 2), " [", round(hr_PRS_interaction_cox_lowerCI, 2), "-", round(hr_PRS_interaction_cox_upperCI, 2), "]")) %>%
  mutate(hr_int = paste0(round(hr_interaction_cox, 2), " [", round(hr_interaction_cox_lowerCI, 2), "-", round(hr_interaction_cox_upperCI, 2), "]")) %>%
  mutate(p_riskfactor_riskfactorOnly_cox_out = format(p_riskfactor_riskfactorOnly_cox, digit=3)) %>%
  mutate(p_interaction_cox_out = format(p_interaction_cox, digit=3)) %>%
  mutate(p_PRSint = format(p_PRS_interaction_cox, digit=3)) %>%
  mutate(p_PRSnonint = format(p_PRS_noninteraction_cox, digit=3)) %>%
  mutate(Risk_factor = recode(Risk_factor,
                              "lnlpa"="Lp(a)",
                              "chol" = "Total cholesterol",
                              "sex" = "Sex",
                              "ldl" = "LDL-C",
                              "CRP" = "CRP",
                              "eversmoker" = "Ever smoker",
                              "T2D_has_disease" = "Type 2 diabetes",
                              "BMI" = "BMI",
                              "FHheart" = "Family history of heart diseases",
                              "trig" = "Triglycerides",
                              "currentsmoker" = "Current smoker",
                              "hdl" = "HDL-C",
                              "enroll_age" = "Age"))

pp = PRS_int_1 %>% filter(p_interaction_cox < 0.05/13)
pp = pp$Risk_factor

PRS_int_22 = PRS_int_1 %>%
  # select(Risk_factor, c_riskfactorOnly_cox, c_interaction_cox) %>%
  mutate(deltac = c_interaction_cox - c_riskfactorOnly_cox) %>%
  mutate(deltac_lowerCI = c_interaction_cox_lowerCI - c_riskfactorOnly_cox_lowerCI) %>%
  mutate(deltac_upperCI = c_interaction_cox_upperCI - c_riskfactorOnly_cox_upperCI) %>%
  mutate(deltac_out = paste(round(deltac, 3), " (", round(deltac_lowerCI, 3), "-", round(deltac_upperCI,3), ")", sep=""))

mm = mean(PRS_int_22$c_riskfactorOnly_cox)

res_sum = c(mean(PRS_int_22$c_riskfactorOnly_cox), sd(PRS_int_22$c_riskfactorOnly_cox))
res_sum_out = res_sum[1]
res_sum_out_lowerCI = res_sum[1] - 1.96*res_sum[2]
res_sum_out_upperCI = res_sum[1] + 1.96*res_sum[2]

print(paste0(rr(res_sum_out, 4), " (", rr(res_sum_out_lowerCI, 4), "-", rr(res_sum_out_upperCI, 4), ")"))

PRS_int_22 %>%
  select(Risk_factor, c_riskfactorOnly_cox_out, c_interaction_cox_out, deltac_out)

PRS_int_2 = PRS_int_1 %>% 
  select(Risk_factor, hr_rf, c_riskfactorOnly_cox_out, c_interaction_cox_out)
PRS_int_2


PRS_int_1 = PRS_int_1 %>% 
  select(Risk_factor, hr_rf, p_riskfactor_riskfactorOnly_cox_out, hr_int_rf, hr_int_PRS, hr_int, p_interaction_cox_out)
PRS_int_1 = PRS_int_1[order(as.numeric(PRS_int_1$p_interaction_cox_out)),]
head(PRS_int_1)
changeSciNot <- function(n) {
  # n = 5e-9
  if (as.numeric(n) > 0.001) 
    return(round(as.numeric(n),3))
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", " x 10", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

PRS_int_IntnoInt = PRS_int_1 %>%
  select(Risk_factor, hr_noint_PRS, p_PRSnonint, hr_int_PRS, p_PRSint)
PRS_int_IntnoInt$p_PRSint = sapply(PRS_int_IntnoInt$p_PRSint, changeSciNot)
PRS_int_IntnoInt$p_PRSnonint = sapply(PRS_int_IntnoInt$p_PRSnonint, changeSciNot)
colnames(PRS_int_IntnoInt) = c("Risk factor", 
                               "HR of Risk factor (95% CI)", "P-value", 
                               "HR of Risk factor (95% CI)", "P-value"
                               )
write.csv(PRS_int_IntnoInt, file="PRSrisk_withandwithoutRiskadj.csv", row.names=F)



colnames(PRS_int_1) = c("Risk factor", "HR of Risk factor (95% CI)", "P-value", "HR of carrier (95% CI)", "HR of PRS", "HR of interaction", "Pinteraction")
PRS_int_1$`P-value` = sapply(PRS_int_1$`P-value`, changeSciNot)
PRS_int_1$Pinteraction = sapply(PRS_int_1$Pinteraction, changeSciNot)
head(PRS_int_1)

# write.csv(PRS_int_1, file="gwPRS_riskInt_summary.csv", row.names=F)
write.csv(PRS_int_1, file="PRSindep_riskInt_summary.csv", row.names=F)




#####################################
multivar_risk = read.csv("multivar_risk.csv") %>%
  mutate(hr_rf_multi = paste0(round(HR, 2), " [", round(lowerCI, 2), "-", round(upperCI, 2), "]")) %>% 
  mutate(p_multi = format(P, digit=3)) %>%
  mutate(Risk_factor = recode(Risk_factor,
                              "lnlpa"="Lp(a)",
                              "chol" = "Total cholesterol",
                              "sex" = "Sex",
                              "ldl" = "LDL-C",
                              "CRP" = "CRP",
                              "eversmoker" = "Ever smoker",
                              "T2D_has_disease" = "Type 2 diabetes",
                              "BMI" = "BMI",
                              "FHheart" = "Family history of heart diseases",
                              "trig" = "Triglycerides",
                              "currentsmoker" = "Current smoker",
                              "hdl" = "HDL-C",
                              "enroll_age" = "Age"))

multivar_risk
multivar_risk = multivar_risk %>%
  select(Risk_factor, hr_rf_multi, p_multi)


PRS_int_uni_multi = merge(PRS_int_1, multivar_risk, by.x='Risk factor', by.y='Risk_factor')
PRS_int_uni_multi = PRS_int_uni_multi %>%
  select("Risk factor", "HR of Risk factor (95% CI)", "P-value", "hr_rf_multi", "p_multi")

colnames(PRS_int_uni_multi)[c(2:5)] = c("HR (95% CI, univariate)", "P-value (univariate)", "HR (95% CI, multivariate)", "P-value (multivariate)")
PRS_int_uni_multi
PRS_int_uni_multi = PRS_int_uni_multi[order(as.numeric(PRS_int_uni_multi$`P-value (multivariate)`)),]
PRS_int_uni_multi$`P-value (multivariate)` = sapply(PRS_int_uni_multi$`P-value (multivariate)`, changeSciNot)
PRS_int_uni_multi
write.csv(PRS_int_uni_multi, "multivar_cor_out.csv", row.names=F)

