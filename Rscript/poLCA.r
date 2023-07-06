run_blrt <- FALSE

pkgs_overall <- c("data.table", "reshape2", "ggplot2", "ggrepel", "scales", "paletteer")
pkgs_lca <- c("poLCA", "networkD3", "scatterpie", "corrplot", "tidyLPA")
sapply(pkgs_overall, require, character.only = TRUE, quietly = TRUE)
source(paste0(run_path, "codes/Revise/FuncsCols_Revise.r"))
sapply(pkgs_lca, require, character.only = TRUE, quietly = TRUE)

load(paste0(out_path, "dat_cleaned.RData"))

vari_lca <- c("age", "sex", "race", "site", "dstage")
dat_lca <- dat_clean[, .SD, .SDcols = c("patient.id", vari_lca)]
dat_lca <- dat_lca[complete.cases(dat_lca)]

### Correlation between indicator variables
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

### Run poLCA
f_lca <- eval(parse(text = paste0("as.formula(cbind(", paste(vari_lca, collapse = ", "), ") ~ 1)"))) #caus.death * year.diag

for (i in n_start:n_end) {
  mi_start <- poLCAParallel::poLCA(formula = f_lca, data = dat_lca, nclass = i, maxiter = 1000, nrep = 30, graph = FALSE, calc.se = TRUE, verbose = FALSE)
  vprobs.start <- poLCA.reorder(mi_start$probs.start, order(mi_start$P, decreasing = TRUE))
  mi <- poLCAParallel::poLCA(formula = f_lca, data = dat_lca, nclass = i, maxiter = 1000, probs.start = vprobs.start, graph = FALSE, calc.se = TRUE, verbose = FALSE)
  save(mi_start, mi, file = paste0(out_path, "lca/lca_m", i, ".RData"))
}

save(dat_lca, file = paste0(out_path, "dat_lca.RData"))

p_use <- "p_lrt"
vy_thres <- 0.5

### LCA performance
res_lca <- data.table()
for (j in n_start:n_end) {
  load(paste0(out_path, "lca/lca_m", j, ".RData"))
  prop <- table(mi$predclass)/mi$Nobs*100
  pre_prop <- mi$P*100

  val_post <- apply(mi$posterior, 1, max)
  sum_post <- summary(val_post)
  if (j == 1) {
    mean_post = 1
  } else {
    dat_post <- data.table(posterior = val_post, class = mi$predclass)
    mean_post <- dat_post[, mean(posterior), by = "class"]$V1
  }

  nadj=(mi$Nobs+2)/24
  sabic = -2 * mi$llik + mi$npar * log(nadj)

  if (file.exists(paste0(out_path, "lca/lca_m", j, ".txt"))) {
    p_file <- readLines(paste0(out_path, "lca/lca_m", j, ".txt"))
    p_blrt <- as.numeric(strsplit(gsub("\"", "", p_file, fixed = FALSE), ": ", fixed = TRUE)[[1]][2])
  } else {
    if (run_blrt & j > 1) {
      p_blrt <- blrt_lca(vnum = j)
      sink(paste0(out_path, "lca/lca_m", j, ".txt"))
      print(paste0("P value for BLRT: ", p_blrt))
      sink()
    } else {
      p_blrt <- NA
    }
  }
  
  res_j <- data.table(n.class = ncol(mi$posterior), sample.size = mi$Nobs, n.param = mi$npar, loglike = mi$llik, AIC = mi$aic, BIC = mi$bic, SABIC = sabic, entropy = poLCA.entropy(mi), p.blrt = p_blrt,
    min.post = as.numeric(sum_post[1]), median.post = as.numeric(sum_post[3]), mean.post = as.numeric(sum_post[4]), max.post = as.numeric(sum_post[6]), mean.per.post = paste(round(mean_post, 2), collapse = "|"), min.per.post = min(mean_post),
    pre.prop = paste(round(pre_prop, 2), collapse = "|"), min.preprop = min(pre_prop), class.prop = paste(round(prop, 2), collapse = "|"), min.prop = min(prop)
  )
  res_lca <- rbind(res_lca, res_j)
}

res_lca[, p_lrt:= NA]
for (n in 2:n_end) {
  plrt <- calc_lrt(n = res_lca$sample.size[n], null_ll = res_lca$loglike[n-1], null_param = res_lca$n.param[n-1], null_classes = res_lca$n.class[n-1], alt_ll = res_lca$loglike[n], alt_param = res_lca$n.param[n], alt_classes = res_lca$n.class[n])[4]
  res_lca$p_lrt[n] <- plrt
}
fwrite(res_lca, file = paste0(out_path, "res_lca.csv"))
