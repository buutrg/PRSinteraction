
library(reshape2)
library(dplyr)
library(tidyr)
library(tableone)
library(data.table)
library(forestplot)
library(plotly)

# setwd("/Users/btruong/Library/CloudStorage/OneDrive-TheBroadInstitute/Works/Broad/Projects/snp_interaction/draft/Tables/")
setwd("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/")

PRS_int = fread("/Volumes/GoogleDrive/My Drive/Works/Broad/Projects/snp_interaction/PRSCs_interaction.txt")
PRS_int = PRS_int %>% filter(Risk_factor != "CKD_has_disease")
PRS_int = PRS_int[order(PRS_int$p_interaction_cox),]

data = read.table("RiskHR_prsStrat.tsv", header=T, sep="\t")
colnames(data)

data$qt = as.factor(data$qt)
data = data %>% filter(Risk_factor != "CKD_has_disease")
data$color = ifelse(PRS_int$p_interaction_cox[match(data$Risk_factor, PRS_int$Risk_factor)] < 0.05/15, 1, 0)
data$color = as.factor(data$color)

intersect(PRS_int$Risk_factor, data$Risk_factor)

colnames(data)
pd = position_dodge(0.5)


risk_list = unique(PRS_int$Risk_factor)
plotlist = list()
for (risk_i in risk_list) {
  p = ggplot(data, aes(x=qt, y=hr_riskfactorOnly_cox, group=Risk_factor)) +
    theme_bw() +
    geom_line(position = pd) +
    geom_errorbar(aes(ymin = hr_riskfactorOnly_cox_lowerCI, ymax = hr_riskfactorOnly_cox_upperCI), position = pd) +
    geom_point(aes(color=color), position = pd) +
    scale_colour_manual(
      values = c("1" = "red", "0" = "black")
    )
  plotlist = append(plotlist, list(p))
}  
ggarrange(plotlist = plotlist)
ggplotly(p)
