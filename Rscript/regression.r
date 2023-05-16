pkgs_overall <- c("data.table", "reshape2", "ggplot2", "ggrepel", "scales", "paletteer")
pkgs_regress <- c("survival", "survminer", "RISCA", "cmprsk", "forestplot") #"fastcmprsk", 
pkgs_treat <- c("survival", "survminer", "forestplot")
sapply(pkgs_overall, require, character.only = TRUE, quietly = TRUE)
sapply(pkgs_regress, require, character.only = TRUE, quietly = TRUE)
sapply(pkgs_treat, require, character.only = TRUE, quietly = TRUE)
source(paste0(run_path, "codes/Revise/FuncsCols_Revise.r"))

load(paste0(out_path, "dat_surv.RData"))

### Cox regression
idat <- dat_treat[vclass == iclass & dsurg_prim == "Surgery" & dstage == istage & site == isite]
ifit <- coxph(Surv(surv.newmonth, surv.new) ~ radiation, data = idat)
res_coef <- summary(ifit)$coef[1, c(1, 3, 4, 5)]
res_conf <- summary(ifit)$conf[1, c(1, 3, 4)]
ires <- c(iclass, istage, isite, table(idat$radiation), res_coef, res_conf)
names(ires) <- c("Class", "Stage", "Site", "nonRadio", "Radio", "beta", "se", "zval", "pval", "hr", "l95", "u95")
ires <- as.data.table(ires)
ires[, hr:= as.numeric(hr)]
ires[, l95:= as.numeric(l95)]
ires[, u95:= as.numeric(u95)]
ires[, pval:= as.numeric(pval)]
ires[, hr_ci:= paste0(sprintf("%.2f", hr), " (", sprintf("%.2f", l95), "-", sprintf("%.2f", u95), ")")]
ires[, pval_use:= sapply(pval, function(x) trans_pval(x))]

### Multi-state model
dat_regress <- copy(dat_surv)
fit_regress_uni <- coxph(Surv(surv.year, caus.death) ~ vclass, dat_regress, id = patient.id)
fit_regress_multi <- coxph(Surv(surv.year, caus.death) ~ vclass + age + sex + race + site + dstage, dat_regress, id = patient.id)

res_crr <- list(dat = dat_regress, fit_uni = fit_regress_uni, fit_multi = fit_regress_multi)
save(res_crr, file = paste0(out_path, "res_crr.RData"))
