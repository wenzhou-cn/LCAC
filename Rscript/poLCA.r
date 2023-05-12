pkgs_lca <- c("poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA")
sapply(pkgs_overall, require, character.only = TRUE, quietly = TRUE)
source(paste0(run_path, "codes/Revise/FuncsCols_Revise.r"))

sapply(pkgs_lca, require, character.only = TRUE, quietly = TRUE)

load(paste0(out_path, "dat_cleaned.RData"))

dat_lca <- dat_clean[, .SD, .SDcols = c("patient.id", vari_lca)]
dat_lca <- dat_lca[complete.cases(dat_lca)]

f_lca <- eval(parse(text = paste0("as.formula(cbind(", paste(vari_lca, collapse = ", "), ") ~ 1)"))) #caus.death * year.diag

for (i in n_start:n_end) {
  mi_start <- poLCAParallel::poLCA(formula = f_lca, data = dat_lca, nclass = i, maxiter = 1000, nrep = 30, graph = FALSE, calc.se = TRUE, verbose = FALSE)
  vprobs.start <- poLCA.reorder(mi_start$probs.start, order(mi_start$P, decreasing = TRUE))
  mi <- poLCAParallel::poLCA(formula = f_lca, data = dat_lca, nclass = i, maxiter = 1000, probs.start = vprobs.start, graph = FALSE, calc.se = TRUE, verbose = FALSE)
  save(mi_start, mi, file = paste0(out_path, "lca/lca_m", i, ".RData"))
}

save(dat_lca, file = paste0(out_path, "dat_lca.RData"))
