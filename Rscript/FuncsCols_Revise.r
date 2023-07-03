############ Functions ############
trans_pval <- function(pval) {
	if (is.na(pval)) {
		vpval <- NA
	} else if (pval == 0) {
		vpval <- "< 0.001"
	} else if (pval > 0 & pval < 0.01) {
		vpval <- gsub("e", "E", format(pval, scientific=T, digits=3))
	#} else if (pval<0.05) {
	#	vpval<-paste0(sprintf("%.3f",pval),"*")
	} else {
		vpval <- sprintf("%.3f", pval)
	}
	return(vpval)
}

get_psurv <- function(vsurv, formated = TRUE) {
  df <- length(vsurv$n) - 1
  p_surv <- pchisq(q = vsurv$chisq, df = df, lower.tail = FALSE)
  if (formated) {
    if (p_surv == 0) {
      p_plot <- "< 2E-16"
    } else if (p_surv < 0.001) {
      p_plot <- gsub("e", "E", format(p_surv, scientific=T, digits=3))
    } else {
      p_plot <- sprintf("%.3f", p_surv)
    }
  } else {
    p_plot <- p_surv
  }
  return(p_plot)
}

cols_sankey <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#EFC000", "#B09C85", "#BC80BD", "#FDB462", "#B3DE69", "#FCCDE5")
col_cms <- c("#DD871BCC", "#356EA1CC", "#D080B1CC", "#2FA88ACC")
col_mol <- ggsci::pal_locuszoom()(4)
col_msi <- ggsci::pal_cosmic("hallmarks_light")(4)[-1]
