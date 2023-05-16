pkgs_overall <- c("data.table", "reshape2", "ggplot2", "ggrepel", "scales", "paletteer")
pkgs_surv <- c("survival", "survminer", "cmprsk", "adjustedCurves")
sapply(pkgs_overall, require, character.only = TRUE, quietly = TRUE)
sapply(pkgs_surv, require, character.only = TRUE, quietly = TRUE)
source(paste0(run_path, "codes/Revise/FuncsCols_Revise.r"))
load(paste0(out_path, "dat_base.RData"))

dat_surv <- copy(dat_base)
dat_surv[, surv.year:= surv.month/12]
dat_surv[caus.death == "Alive", status:= 0]
dat_surv[caus.death == "CRC", status:= 1]
dat_surv[caus.death == "nonCRC", status:= 2]
dat_surv[, caus.death:= factor(caus.death, levels = c("Alive", "CRC", "nonCRC"))]

### K-M CIF
fit_overall <- eval(parse(text=paste0("survfit(as.formula(Surv(surv.year, surv.death) ~ ", igroup, "), data = dat_surv)")))
if (length(fit_overall$strata) > 0) {
  fit_strata <- c()
  for (i in 1:length(fit_overall$strata)) {
    fit_strata <- c(fit_strata, rep(gsub(paste0(igroup, "="), "", names(fit_overall$strata)[i]), fit_overall$strata[i]))
  }
} else {
  fit_strata <- NA
}
dat_fit_overall <- data.table(time = fit_overall$time, event = fit_overall$n.event, surv = fit_overall$surv, strata = fit_strata)
dat_fit_overall[, cif:= 1 - surv]

### Aalen-Johansen CIF (multi-state)
fit_aj <- eval(parse(text=paste0("survfit(as.formula(Surv(surv.year, caus.death) ~ ", igroup, "), data = dat_surv)")))
#fit_aj <- eval(parse(text=paste0("survfit(as.formula(Surv(surv.year, status, type = 'mstate') ~ ", igroup, "), data = dat_surv)")))
if (length(fit_aj$strata) > 0) {
  fit_strata_aj <- c()
  for (j in 1:length(fit_aj$strata)) {
    fit_strata_aj <- c(fit_strata_aj, rep(gsub(paste0(igroup, "="), "", names(fit_aj$strata)[j]), fit_aj$strata[j]))
  }
} else {
  fit_strata_aj <- NA
}
dat_fit_aj_raw <- data.table(time = fit_aj$time)
dat_fit_aj_raw <- cbind(dat_fit_aj_raw, fit_aj$n.event[, -1])
setnames(dat_fit_aj_raw, c(paste0("V", 1:(ncol(fit_aj$n.event) - 1))), c(paste0("event_", 1:(ncol(fit_aj$n.event) - 1))))
dat_fit_aj_raw <- cbind(dat_fit_aj_raw, fit_aj$pstate[, -1])
setnames(dat_fit_aj_raw, c(paste0("V", 1:(ncol(fit_aj$pstate) - 1))), c(paste0("cif_", 1:(ncol(fit_aj$pstate) - 1))))
dat_fit_aj_raw <- cbind(dat_fit_aj_raw, data.table(strata = fit_strata_aj))

if (ncol(fit_aj$n.event) > 1) {
  dat_fit_aj <- data.table()
  for (j in 1:(ncol(fit_aj$n.event) - 1)) {
    dat_j <- dat_fit_aj_raw[, .SD, .SDcols = c("time", paste0(c("event_", "cif_"), j), "strata")]
    setnames(dat_j, paste0(c("event_", "cif_"), j), c("event", "cif"))
    dat_j[, group:= j]
    dat_fit_aj <- rbind(dat_fit_aj, dat_j)
  }
} else {
  dat_fit_aj <- copy(dat_fit_aj_raw)
}
#setnames(dat_fit_aj, "strata", igroup)
dat_fit_aj[, group:= factor(group, levels = c(1, 2), labels = c("CRC", "nonCRC"))]

### Adjust LCAC
if (length(summary(dat_surv[, .SD, .SDcols = igroup])) == 2) {
  fit_psm <- eval(parse(text=paste0("glm(as.formula(", igroup, " ~ vclass), data=dat_surv, family='binomial'(link='logit'))")))
} else {
  fit_psm <- eval(parse(text=paste0("nnet::multinom(as.formula(", igroup, " ~ vclass), data=dat_surv)")))
}
fit_adj <- adjustedcif(data=dat_surv, variable=igroup, ev_time="surv.year", event="surv.death", method="iptw", cause=1, treatment_model=fit_psm, conf_int=FALSE)
dat_fit_iptw <- as.data.table(fit_adj$adjcif)
setnames(dat_fit_iptw, "group", "strata")

### Adjust variables using adjustedCurves
vari_org <- c("age", "sex", "race", "site", "dstage")
if (length(grep(igroup, vari_org)) > 0) {
  vari_use <- vari_org[-which(vari_org == igroup)]
} else {
  vari_use <- vari_org
}
formu_vari <- paste(vari_use, collapse = " + ")

if (length(summary(dat_surv[, .SD, .SDcols = igroup])) == 2) {
  vwidth <- 12
  fit_psm2 <- eval(parse(text=paste0("glm(as.formula(", igroup, " ~ ", formu_vari, "), data=dat_surv, family='binomial'(link='logit'))")))
} else {
  vwidth <- 18
  fit_psm2 <- eval(parse(text=paste0("nnet::multinom(as.formula(", igroup, " ~ ", formu_vari, "), data=dat_surv)")))
}
fit_adj2 <- adjustedcif(data=dat_surv, variable=igroup, ev_time="surv.year", event="surv.death", method="iptw", cause=1, treatment_model=fit_psm2, conf_int=FALSE)
dat_fit_iptw2 <- as.data.table(fit_adj2$adjcif)
setnames(dat_fit_iptw2, "group", "strata")

dat_fit_overall[, type:= "K-M"]
dat_fit_aj[, type:= "A-J"]
dat_fit_iptw[, type:= "LCA"]
dat_fit_iptw2[, type:= "Adjusted"]

dat_comp_group <- rbind(dat_fit_overall[, .(time, cif, strata, type)], dat_fit_aj[group == "CRC", .(time, cif, strata, type)], dat_fit_iptw[, .(time, cif, strata, type)], dat_fit_iptw2[, .(time, cif, strata, type)])
dat_comp_group[, type:= factor(type, levels = c("A-J", "K-M", "Adjusted", "LCA"), labels = c("Aalen-Johansen", "Kaplan-Meier", "Variable-adjusted", "LCAC-adjusted"))]
cols_comp <- cols_sankey[c(7, 6, 8, 9)]
plot_cmp_adjust2 <- ggplot(data = dat_comp_group, aes(x = time, y = cif, color = type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols_comp) +
  theme_bw() +
  facet_grid(.~strata) +
  labs(x = "Years after diagnosis", y = "Cumulative incidence", color = "") +
  theme(axis.title=element_text(size=18), axis.text=element_text(size=13), legend.text=element_text(size=13), strip.text=element_text(size=13), panel.grid=element_blank())
ggsave(plot_cmp_adjust2, file = paste0(out_path, "plots/survival/CompFacet_cmprsk_", igroup, "_nclass", n_use, "_adjust.jpeg"), dpi=600, width=vwidth, height=6, device="jpeg")
