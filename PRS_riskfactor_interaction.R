


library("optparse")
# ==========================================================================

# ==========================================================================
option_list = list(
  make_option(c("--snp"), type = "character"),
  make_option(c("--out"), type = "character"),
  make_option(c("--isPRS"), type = "logical")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);





library(data.table)
library(stringr)
library(dplyr)
library(survival)
library(survminer)


options(datatable.fread.datatable=FALSE)


# opt = data.frame(snp="rs9349379", isPRS=F, out="rs9349379_interaction.txt", stringsAsFactors=F)
# opt = data.frame(snp="rs10455872", isPRS=F, out="rs10455872_interaction.txt", stringsAsFactors=F)
# opt = data.frame(snp="ukb_cardiogramSNP", isPRS=T, out="cardiogramSNP_interaction.txt", stringsAsFactors=F)
opt = data.frame(snp="PRSCS", isPRS=T, out="PRSCs_interaction.txt", stringsAsFactors=F)


setwd("/medpop/esp2/btruong/Projects/snp_risk_interaction/data")


source('/medpop/esp2/btruong/scripts/master_functions.R')


allcleanpheno = as.data.frame(fread("/medpop/esp2/aniruddh/SouthAsians/allcleanpheno.txt",header=T))


pca = as.data.frame(fread("/medpop/esp2/btruong/Projects/snp_risk_interaction/data/ukb_sqc_v2_7089.tsv", header=T))


# if (opt$isPRS) {
# 	if (opt$snp == 'ukb_cardiogramSNP') {
# 			prs = read.table(paste0('/medpop/esp2/btruong/Projects/snp_risk_interaction/data/ukb_cad_snps/ukb_cardiogramSNP.profile'), header=T, stringsAsFactors=F)
# 			prs = prs %>% rename(PRS = "SCORESUM")
# 		} else {
# 			prs = read.table(paste0("/medpop/esp2/yruan/projects/tgilliland.prs/data/sumprs/PRS-CS.CAD_phi1e-04.sumprs"), header=T)
# 			prs = prs %>%
# 				rename("PRS"='sum.prs')}
# } else {
# 	prs = read.table(paste0('/medpop/esp2/btruong/Projects/snp_risk_interaction/data/ukb_cad_snps/', opt$snp, '.profile'), header=T, stringsAsFactors=F)
# 	prs = prs %>% rename(PRS = 'CNT2')
# }

prs = fread(paste0("/medpop/esp2/wallace/projects/Amit/PRS/MultiAncestryPRS/Meta2021CAD_noUKB/UKB_LDref/Scores_LDpred2_UKB.data/results/CAD-80000-ldpred2.scores.txt.gz"))
colnames(prs)[1] = "IID"
colnames(prs)[2:ncol(prs)] = paste0("prs", 1:(ncol(prs)-1))


# CAD_data = as.data.frame(fread("/medpop/esp2/aniruddh/Lpa/UKBB-Coronary_artery_disease_HARD_202006.tsv"))
CAD_data = as.data.frame(fread("/medpop/esp2/projects/UK_Biobank/Phenotype_Library/phenoV3_r202006/disease/Coronary_Artery_Disease_SOFT.tab.tsv.gz"))


CAD_data_df = CAD_data[, c("sample_id", "has_disease", "incident_disease", "prevalent_disease", "censor_age")]
CAD_data_df = CAD_data_df %>%
	rowwise() %>%
	mutate(incident_disease = ifelse(is.na(incident_disease), 0, incident_disease)) %>%
	mutate(prevalent_disease = ifelse(is.na(prevalent_disease), 0, prevalent_disease))
colnames(CAD_data_df)[2:5] = paste0("CAD_", colnames(CAD_data_df)[2:5])



############################ 

treatment_field = fread("/medpop/esp2/btruong/Projects/tentative/colchicine/data/treatment_field_6153.txt")
treatment_field_untreat = apply(treatment_field[,2:5], 1, function(x) all(x %in% c(NA, "-7")))
notx_iid = treatment_field$eid[which(treatment_field_untreat)]



tot_race = as.data.frame(fread('/medpop/esp2/aniruddh/Lpa/LpaPrimSecDatasetnew.txt', check.names=F))
tot_race = tot_race[which(tot_race$race ==1), ]

tot_race = tot_race %>%
	rename(any_of(c(sex="X31.0.0")))

tot_race = merge(tot_race, CAD_data_df, by.x="sample_id", by.y='sample_id', all.x=T)
tot_race = merge(tot_race, pca, by.x="sample_id", by.y='eid', all.x=T)

tot_df = merge(tot_race, prs, by.x='sample_id', by.y='IID', all.x=T)
tot_df$notreat = tot_df$sample_id %in% notx_iid 

allcleanpheno1 = allcleanpheno[, -which(colnames(allcleanpheno) %in% intersect(colnames(allcleanpheno), colnames(tot_df)))]
tot_df = merge(tot_df, allcleanpheno1, by.x='sample_id', by.y='eid', all.x=T)


tot_df$difference = tot_df$CAD_censor_age - tot_df$enroll_age
tot_df$prevalent_disease = tot_df$CAD_prevalent_disease
tot_df$incident_disease = tot_df$CAD_incident_disease

tot_df = tot_df %>%
	filter(lpa>0)

# r2 = NULL

# for (score_i in 2:ncol(prs)) {
# 	# score_i = 2
# 	mod = lm(as.formula(paste0("CAD_incident_disease~", colnames(prs)[score_i], " + enroll_age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=tot_df)
# 	a = summary(mod)
# 	r2=c(r2, a$r.squared)
# }

# nn = colnames(prs)[which.max(r2)+1]
nn = "prs27"
colnames(tot_df)[which(colnames(tot_df)==nn)] = "PRS"

cat_var = c("eversmoker", "currentsmoker", "DM", "FHheart")
cont_var = c("lnlpa", "hdl", "ldl", "trig", "chol", "BMI", "CRP", "SBP", "DBP", "enroll_age", "sex")
variable_list = c(cat_var, cont_var)
# variable_list = c("Overall", "eversmoker", "currentsmoker", "lpanew", "hdladj", "ldl", "trigadj", "choladj", "DM", "BMI", "CRP", 'FHheart')

covariates_cont = c("enroll_age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
covariates_cat = c("sex", "genotyping_array")
covariates = c(covariates_cont, covariates_cat)



other_covar = c("prevalent_disease", "incident_disease", 'difference', "PRS", 'anylipidmed0', 'statin0')
# other_covar = c("prevalent_disease", "incident_disease", 'difference', 'anylipidmed0', 'statin0')

tot_df = tot_df[, unique(c(variable_list, covariates, other_covar, "notreat"))]
# tot_df = tot_df[complete.cases(tot_df), ]


tot_df = tot_df %>% 
	filter(anylipidmed0 == 0)
tot_df$txstt = 1-as.numeric(tot_df$notreat)
dim(tot_df)

e1 = list()


g_cox = tot_df %>% 
	filter(prevalent_disease==0)
g_cox[, c(covariates_cont)] = scale(g_cox[, c(covariates_cont)])
g_cox_save = g_cox
if (opt$isPRS) 
	g_cox$PRS = scale(g_cox$PRS)

g_logistic = tot_df %>% 
	filter(incident_disease==0)
g_logistic[, c(covariates_cont)] = scale(g_logistic[, c(covariates_cont)])
g_logistic_save = g_logistic
if (opt$isPRS) 
	g_logistic$PRS = scale(g_logistic$PRS)
	

summary(g_cox[, covariates_cont])


e2 = list()

formula = as.formula(paste0("Surv(difference, incident_disease) ~ PRS + enroll_age + sex + genotyping_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))



mod = coxph(formula, data=g_cox)
f_cox = summary(mod)
f_cox
e2$snp = opt$snp
e2$Risk_factor = "Overall"
e2$incident_rate = paste0(sum(g_cox$incident_disease),' / ',nrow(g_cox),' (',round(sum(g_cox$incident_disease)/nrow(g_cox)*100,1),' %)')
tmp = f_cox$coefficients
e2$p_interaction_cox = tmp[1, 5]
e2$coef_interaction_cox = tmp[1, 1]





formula = as.formula(paste0("prevalent_disease ~ PRS + enroll_age + sex + genotyping_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
mod_log = glm(formula, data=g_logistic, family="binomial")
f_log = summary(mod_log)
f_log
e2$prevalent_rate = paste0(sum(g_logistic$prevalent_disease),' / ',nrow(g_logistic),' (',round(sum(g_logistic$prevalent_disease)/nrow(g_logistic)*100,1),' %)')
tmp = f_log$coefficients
e2$p_interaction_log = tmp[2, 4]
e2$coef_interaction_log = tmp[2, 1]
e2


rr = function(a,b=2) round(a,b)

hrest = function(coef, se) {
	mm = exp(coef)
	lowerCI = mm - 1.97*se
	upperCI = mm + 1.97*se
	return(paste0(rr(mm), " [", rr(lowerCI), "-", rr(upperCI), "]"))
}



e1 = list(e2)
e1_save = e1
not_transform = c("lnlpa")
e1 = e1_save
cox_res = NULL

for (variable_i in 1:length(variable_list)) {

	# variable_i = 8
	variable = variable_list[variable_i]
	print(variable)
	
	
	formula = as.formula(paste0("Surv(difference, incident_disease) ~ PRS*", variable, " + txstt + enroll_age + sex + genotyping_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
	
	g_cox_new = g_cox
	if (variable %in% c("ldl")) {
		g_cox_new = g_cox %>%
			filter(anylipidmed0==0)
	} else if (variable %in% cont_var && (!variable %in% not_transform)) {
		g_cox_new[, variable] = irnt(g_cox_new[, variable])	
	}
		
	g_logistic_new = g_logistic
	if (variable %in% c("ldl")) {
		g_logistic_new = g_logistic %>%
			filter(anylipidmed0==0)
		} else if (variable %in% cont_var && (!variable %in% not_transform)) {
			g_logistic_new[, variable] = irnt(g_logistic_new[, variable])	
		}
	
	mod = coxph(formula, data=g_cox_new)
	tmp = summary(mod)$coefficients
	f = summary(mod)
	f
	e2 = data.frame(
			snp = opt$snp,
			Risk_factor = variable_list[variable_i],
			incident_rate = paste0(sum(g_cox_new$incident_disease),' / ',nrow(g_cox_new),' (',round(sum(g_cox_new$incident_disease)/nrow(g_cox_new)*100,1),' %)'),
			coef_interaction_cox = tmp[nrow(tmp), 1], p_interaction_cox = tmp[nrow(tmp), 5], hr_cox = hrest(tmp[nrow(tmp), 1], tmp[nrow(tmp), 3]),
			coef_prs_cox = tmp[1, 1], p_prs_cox = tmp[1, 5], hr_prs_cox = hrest(tmp[1, 1], tmp[1, 3]),
			coef_riskfactor_cox = tmp[2, 1], p_riskfactor_cox = tmp[2, 5], hr_riskfactor_cox = hrest(tmp[2, 1], tmp[2, 3])
			
	)
	cox_res = rbind(cox_res, e2)
	
	formula = as.formula(paste0("prevalent_disease ~ PRS*", variable, " + enroll_age + sex + genotyping_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
	
		
	mod_log = glm(formula, data=g_logistic_new, family="binomial")
	f_log = summary(mod_log)
	f_log
	e2$prevalent_rate = paste0(sum(g_logistic_new$prevalent_disease),' / ',nrow(g_logistic_new),' (',round(sum(g_logistic_new$prevalent_disease)/nrow(g_logistic_new)*100,1),' %)')
	tmp = f_log$coefficients
	e2$p_interaction_log = tmp[nrow(tmp), 4]
	e2$coef_interaction_log = tmp[nrow(tmp), 1]
	e2
	
	e1 = append(e1, list(e2))
	

}

e1 = do.call(bind_rows, e1)
e1$HR = exp(e1$coef_interaction_cox)
e1_out = e1[order(e1$p_interaction_cox),]
e1_out = as.data.frame(e1_out)
e1_out = e1_out[!is.na(e1_out$hr_prs_cox),]
e1_out

cox_res = cox_res[order(cox_res$p_interaction_cox),]

cox_res_out = cox_res
cox_res_out$p_interaction_cox = format.pval(cox_res_out$p_interaction_cox, digits=1, eps=1e-16)
cox_res_out = cox_res_out %>%
	mutate(risk = recode(Risk_factor,
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
cox_res_out = cox_res_out %>% select(risk, hr_riskfactor_cox, hr_prs_cox, hr_cox, p_interaction_cox)
cox_res_out

fwrite(cox_res_out, "cox_res_out_noANYtx_adjCovar.txt", row.names=F, quote=F, sep="\t")

if (opt$isPRS) {
	e1_out$gene = NA
	e1_out$beta_cardiogram = NA
	e1_out$AF = NA
	e1_out$SE_cardiogram = NA
	e1_out$P_cardiogram = NA
} else {
	gwas_indep = as.data.frame(fread("gwas_indep.txt"))
	gwas_indep_sel = gwas_indep %>% filter(snptestid == opt$snp)
	e1_out$gene = as.character(gwas_indep_sel['RefSeq.Gene.name'])
	e1_out$beta_cardiogram = as.numeric(gwas_indep_sel['logOR'])
	e1_out$AF = as.numeric(gwas_indep_sel['effect_allele_freq'])
	e1_out$SE_cardiogram = as.numeric(gwas_indep_sel['se_gc'])
	e1_out$P_cardiogram = as.numeric(gwas_indep_sel['p-value_gc'])	
}

e1_out = e1_out[c("snp", 'gene', 'Risk_factor', 'incident_rate', 'prevalent_rate', 'AF', 'beta_cardiogram', 'SE_cardiogram', 'P_cardiogram', 
	"coef_interaction_log", "p_interaction_log", 
	"coef_interaction_cox", "p_interaction_cox")]
e1_out

# e1_out_df = as.data.frame(e1_out)
# e1_out_df1 = e1_out_df[order(unlist(e1_out_df$p_interaction_cox)), ]

fwrite(e1_out, file=opt$out, row.names=F, quote=F, sep="\t")










