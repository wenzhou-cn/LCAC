pkgs_overall <- c("data.table", "reshape2", "ggplot2", "ggrepel", "scales", "paletteer")
sapply(pkgs_overall, require, character.only = TRUE, quietly = TRUE)
source(paste0(run_path, "codes/Revise/FuncsCols_Revise.r"))

pkgs_surv <- c("poLCA", "survival", "survminer", "cmprsk")
sapply(pkgs_surv, require, character.only = TRUE, quietly = TRUE)

load("dat_surv_TCGA.RData")

### Original survival
ifit <- eval(parse(text = paste0("survfit(as.formula(Surv(", surv.type, ".new.year, ", surv.type, ".new) ~ group), data = dat_score)")))
idiff <- eval(parse(text = paste0("survdiff(as.formula(Surv(", surv.type, ".new.year, ", surv.type, ".new) ~ group), data = dat_score)")))
iplgr <- get_psurv(idiff)
plg_sig <- ifelse(iplgr < 0.001, "***", ifelse(iplgr < 0.01, "**", ifelse(iplgr < 0.05, "*", "")))

icox <- eval(parse(text = paste0("coxph(as.formula(Surv(", surv.type, ".new.year, ", surv.type, ".new) ~ group), data = dat_score)")))
icoef <- as.data.table(summary(icox)$coef, keep.rownames = TRUE)
setnames(icoef, "Pr(>|z|)", "pval")
icoef[, psig:= ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))]
iconf <- as.data.table(summary(icox)$conf, keep.rownames = TRUE)

### LCAC-adjusted survival
idat_adj <- copy(dat_score)
index <- with(idat_adj, cbind(as.numeric(vclass), as.numeric(group)))
tab1 <- with(idat_adj, table(vclass))/ nrow(idat_adj)
tab2 <- with(idat_adj, table(vclass, group))/nrow(idat_adj)
tab3 <- with(idat_adj, table(group)) / nrow(idat_adj)
rwt <- rep(tab1, length(tab3))/tab2
idat_adj$rwt <- rwt[index] # add per subject weights to the data set

wtscale <- table(idat_adj$group)/tapply(idat_adj$rwt, idat_adj$group, sum)
idat_adj$wt_use <- c(idat_adj$rwt * wtscale[idat_adj$group])
ifit_adj <- eval(parse(text = paste0("survfit(as.formula(Surv(", surv.type, ".new.year, ", surv.type, ".new) ~ group), data = idat_adj, weight=wt_use)")))
iid <- 1:nrow(idat_adj)
icox_adj <- eval(parse(text = paste0("coxph(as.formula(Surv(", surv.type, ".new.year, ", surv.type, ".new) ~ group), data = idat_adj, cluster=iid, weight=wt_use)")))
iplgr_adj <- trans_pval(summary(icox_adj)$robscore["pvalue"])
plg_sig_adj <- ifelse(iplgr_adj < 0.001, "***", ifelse(iplgr_adj < 0.01, "**", ifelse(iplgr_adj < 0.05, "*", "")))
icoef_adj <- as.data.table(summary(icox_adj)$coef, keep.rownames = TRUE)
setnames(icoef_adj, "Pr(>|z|)", "pval")
icoef_adj[, psig:= ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))]
iconf_adj <- as.data.table(summary(icox_adj)$conf, keep.rownames = TRUE)
