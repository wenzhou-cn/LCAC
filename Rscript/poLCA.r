pkgs_lca <- c("poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA")
sapply(pkgs_overall, require, character.only = TRUE, quietly = TRUE)
source(paste0(run_path, "codes/Revise/FuncsCols_Revise.r"))

sapply(pkgs_lca, require, character.only = TRUE, quietly = TRUE)

load(paste0(out_path, "dat_cleaned.RData"))

vari_lca <- c("age", "sex", "race", "site", "dstage")
dat_lca <- dat_clean[, .SD, .SDcols = c("patient.id", vari_lca)]
dat_lca <- dat_lca[complete.cases(dat_lca)]

res_corr <- data.frame()
for (ivari in vari_lca) {
  icor <- c()
  for (jvari in vari_lca) {
    if (ivari == jvari) {
      jcor <- 1
    } else {
      jdat <- copy(dat_lca)
      setnames(jdat, c(ivari, jvari), c("xvari", "yvari"))
      jcor <- DescTools::UncertCoef(jdat$xvari, jdat$yvari, direction = "row")
    }
    icor <- cbind(icor, jcor)
  }
  res_corr <- rbind(res_corr, icor)
}
colnames(res_corr) <- c("Age", "Sex", "Race", "Site", "Stage")
rownames(res_corr) <- c("Age", "Sex", "Race", "Site", "Stage")

jpeg(file = paste0(out_path, "lca/corrplot.jpeg"), res = 600, width = 10, height = 10, units = "in")
corrplot(as.matrix(res_corr), col.lim = c(0, 1), method = "color", col = COL1('YlGn'), tl.col = 'black', tl.srt = 45, addgrid.col = 'white', addCoef.col = 'grey50', number.digits = 3)
dev.off()

f_lca <- eval(parse(text = paste0("as.formula(cbind(", paste(vari_lca, collapse = ", "), ") ~ 1)"))) #caus.death * year.diag

for (i in n_start:n_end) {
  mi_start <- poLCAParallel::poLCA(formula = f_lca, data = dat_lca, nclass = i, maxiter = 1000, nrep = 30, graph = FALSE, calc.se = TRUE, verbose = FALSE)
  vprobs.start <- poLCA.reorder(mi_start$probs.start, order(mi_start$P, decreasing = TRUE))
  mi <- poLCAParallel::poLCA(formula = f_lca, data = dat_lca, nclass = i, maxiter = 1000, probs.start = vprobs.start, graph = FALSE, calc.se = TRUE, verbose = FALSE)
  save(mi_start, mi, file = paste0(out_path, "lca/lca_m", i, ".RData"))
}

save(dat_lca, file = paste0(out_path, "dat_lca.RData"))
