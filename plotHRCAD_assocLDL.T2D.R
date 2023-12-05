

library(data.table)
library(stringr)
library(dplyr)
library(survival)
library(survminer)
library(tableone)
library(forestplot)
library(formattable)
library(ggplot2)
library(plyr)

options(datatable.fread.datatable=FALSE)

setwd("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/")
data = fread("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/res_HR.CAD_logor.T2D_beta.LDL_noint.csv")

data = data %>%
  mutate(color = ifelse(Risk_factor_p_int_cox_ldl <= 0.05, "P<=0.05", "P>0.05")) %>%
  mutate(label = ifelse(Risk_factor_p_int_cox_ldl <= 0.05 & hr_PRS_cox > 1.2, RefSeq.Gene.name, ""))
table(data$color)
dim(data)
head(data)

data1 = data %>% filter(color == "P<=0.05")
data2 = data %>% filter(color != "P<=0.05")


cor.test(data$coef_LDL_log, data$hr_PRS_cox)
cor.test(data1$coef_LDL_log, data1$hr_PRS_cox)
cor.test(data2$coef_LDL_log, data2$hr_PRS_cox)


pdf("LDL_CAD.pdf", width = 10, height = 8)
ggplot(data=data, aes(x=coef_LDL_log, y=hr_PRS_cox, color=color)) +
  geom_point(size=3) +
  # geom_label() +
  # geom_text(check_overlap = TRUE) +
  geom_smooth(method='lm', formula=y~x, se=F) + 
  ggrepel::geom_label_repel(aes(label=label), size=3, show.legend = FALSE) +
  labs(x="Association with LDL (Beta)", y="Association with CAD (Hazard ratio)", color="Interaction P-value") +
  theme_bw(base_size = 20)
dev.off()

data = data %>%
  mutate(color = ifelse(Risk_factor_p_int_cox_T2D_has_disease <= 0.05, "P<=0.05", "P>0.05")) %>%
  mutate(label = ifelse(Risk_factor_p_int_cox_T2D_has_disease <= 0.05 & hr_PRS_cox > 1.2, RefSeq.Gene.name, ""))
table(data$color)
dim(data)
head(data)


data1 = data %>% filter(color == "P<=0.05")
data2 = data %>% filter(color != "P<=0.05")


cor.test(data$coef_LDL_log, data$hr_PRS_cox)
cor.test(data1$coef_LDL_log, data1$hr_PRS_cox)
cor.test(data2$coef_LDL_log, data2$hr_PRS_cox)



pdf("T2D_CAD.pdf", width = 10, height = 8)
ggplot(data=data, aes(x=exp(coef_T2D_log), y=hr_PRS_cox, color=color)) +
  geom_point(size=3) +
  # geom_label() +
  # geom_text(check_overlap = TRUE) +
  geom_smooth(method='lm', formula=y~x, se=F) + 
  ggrepel::geom_label_repel(aes(label=label), size=3, show.legend = FALSE) +
  labs(x="Association with T2D (Odds ratio)", y="Association with CAD (Hazard ratio)", color="Interaction P-value") +
  theme_bw(base_size = 20)
dev.off()



