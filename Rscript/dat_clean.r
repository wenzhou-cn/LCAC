
library(data.table)

dat_org0 <- fread("SEER_2022Sub_17Res_CRC.csv")
vari_tab <- fread("vari_table.csv")
dat_org <- dat_org0[, .SD, .SDcols = vari_tab[use == 1]$name_org]
setnames(dat_org, vari_tab[use == 1]$name_org, vari_tab[use == 1]$name_set)

############ Data filtering ############
dat_filter1 <- dat_org[behavior == "Malignant"]
dat_filter2 <- dat_filter1[seq.num %in% c("One primary only", "1st of 2 or more primaries") & first.malig == "Yes"]
dat_filter3 <- dat_filter2[diag.confirm == "Positive histology"]
dat_filter3[, age_num:= as.numeric(ifelse(age_org == "Unknown", NA, ifelse(age_org == "90+ years", 90, gsub(" years", "", age_org, fixed = TRUE))))]
dat_filter4 <- dat_filter3[age_num >= 18]
dat_filter5 <- dat_filter4[race != "Non-Hispanic Unknown Race"]

### Re-define Stage
dat_filter5[year.diag >= 2018, dstage:= dStage.2018]
dat_filter5[year.diag %in% c(2016, 2017), dstage:= dStage.2016]
dat_filter5[year.diag >= 2010 & year.diag <= 2015, dstage:= dStage.2010]
dat_filter5[year.diag >= 2004 & year.diag <= 2009, dstage:= dStage.2004]
dat_filter5[year.diag <= 2003, dstage:= ajcc3.Stage.1988]

### Re-define T
dat_filter5[year.diag >= 2018, dT:= dT.2018]
dat_filter5[year.diag %in% c(2016, 2017), dT:= dT.2016]
dat_filter5[year.diag >= 2010 & year.diag <= 2015, dT:= dT.2010]
dat_filter5[year.diag >= 2004 & year.diag <= 2009, dT:= dT.2004]
dat_filter5[year.diag <= 2003, dT:= dT.1988]

dat_filter6 <- dat_filter5[!(dstage %in% c(88, 90, 98, 99, "NA", "Not applicable", "UNK Stage", "0", "Stage 0") | dT %in% c("c0", "p0", "pIS", "T0", "Tis", "Tis(LAMN)"))]
dat_filter7 <- dat_filter6[!(site %in% c("Large Intestine, NOS", "Appendix"))]
dat_filter8 <- dat_filter7[surv.month != "Unknown"]
dat_final <- dat_filter8[!(cod.site == "State DC not available or state DC available but no COD")]

dat_final[year.diag <= 2017, grade:= grade.2017]
dat_final[year.diag >= 2018, grade:= grade.path.2018]
save(dat_final, file = paste0(out_path, "dat_finaluse.RData"))
