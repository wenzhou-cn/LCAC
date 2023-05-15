
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

vari_use <- c("patient.id", "age_num", "sex", "race", "marital", "rural", "income", 
	"year.diag", "site", "behavior", "histology", "dstage", "dT", "dN", "dM", "grade",
	"radiation", "chemotherapy", "dsurg_prim", "surg_rad", "surg_sys", "reason_nosurg",
	"monDiag2Treat", "cea.pre", "per.inva", "tumor.depo", "cod.site",
	"surv.month", "vital"
)
dat_use <- dat_final[, .SD, .SDcols = vari_use]

############ Data re-code ############
dat_clean <- copy(dat_use)
dat_clean[age_num >= 18 & age_num <= 44, age:= "18-44 years"]
dat_clean[age_num >= 45 & age_num <= 69, age:= "45-69 years"]
dat_clean[age_num >= 70, age:= "70+ years"]
dat_clean[, age:= factor(age, levels = c("18-44 years", "45-69 years", "70+ years"))]

dat_clean[, sex:= factor(sex, levels = c("Female", "Male"))]
#dat_clean[, race:= factor(race, levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic (All Races)", "Non-Hispanic Asian or Pacific Islander", "Non-Hispanic American Indian/Alaska Native"), labels = c("White", "Black", "Hispanic", "Asian", "Native"))]
dat_clean[, marital:= factor(marital, levels = c("Married (including common law)", "Unmarried or Domestic Partner", "Separated", "Single (never married)", "Divorced", "Widowed"), labels = c("Married", "Partner", "Separated", "Single", "Divorced", "Widowed"))]
dat_clean[, rural:= factor(rural, levels = c("Counties in metropolitan areas ge 1 million pop", "Counties in metropolitan areas of 250,000 to 1 million pop", "Counties in metropolitan areas of lt 250 thousand pop", "Nonmetropolitan counties adjacent to a metropolitan area", "Nonmetropolitan counties not adjacent to a metropolitan area", "Unknown/missing/no match (Alaska or Hawaii - Entire State)"), labels = c(">1Mmetro", "<1Mmetro", "<250Kmetro", "adjMetro", "nonMetro", "Alaska"))]
dat_clean[, chemotherapy:= factor(chemotherapy, levels = c("No/Unknown", "Yes"), labels = c("No", "Yes"))]
dat_clean[, cea.pre:= factor(cea.pre, levels = c("CEA negative/normal; within normal limits", "Borderline", "CEA positive/elevated"), labels = c("Negative", "Borderline", "Positive"))]
dat_clean[, per.inva:= factor(per.inva, levels = c("Perineural invasion not identified/not present", "Perineural invasion identified/present"), labels = c("nopresent", "present"))]
dat_clean[, vital:= factor(vital, levels = c("Alive", "Dead"))]

dat_clean[race %in% c("Non-Hispanic Asian or Pacific Islander", "Non-Hispanic American Indian/Alaska Native"), race:= "Other"]
dat_clean[, race:= factor(race, levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic (All Races)", "Other"), labels = c("White", "Black", "Hispanic", "Other"))]

dat_clean[income %in% c("$35,000 - $39,999", "$40,000 - $44,999", "$45,000 - $49,999"), income:= "$35,000 - $49,999"]
dat_clean[income %in% c("$50,000 - $54,999", "$55,000 - $59,999", "$60,000 - $64,999", "$65,000 - $69,999", "$70,000 - $74,999"), income:= "$50,000 - $74,999"]
dat_clean[, income:= factor(income, levels = c("< $35,000", "$35,000 - $49,999", "$50,000 - $74,999", "$75,000+"), labels = c("< 35K", "35K-49K", "50K-74K", "75K+"))]

dat_clean[site %in% c("Cecum", "Ascending Colon", "Hepatic Flexure", "Transverse Colon"), site:= "Right.sided.Colon"]
dat_clean[site %in% c("Splenic Flexure", "Descending Colon", "Sigmoid Colon"), site:= "Left.sided.Colon"]
dat_clean[site %in% c("Rectosigmoid Junction", "Rectum"), site:= "Rectum"]
dat_clean[, site:= factor(site, levels = c("Right.sided.Colon", "Left.sided.Colon", "Rectum"))]

dat_clean[, histology:= as.character(histology)]
dat_clean[histology %in% as.character(8140:8147), histology:= "adenocarcinoma"]
dat_clean[histology %in% as.character(8260:8265), histology:= "papillary"]
dat_clean[histology %in% as.character(8210:8213), histology:= "adenomatous.polyp"]
dat_clean[histology %in% as.character(8480:8482), histology:= "mucinous"]
dat_clean[histology %in% as.character(8240:8249), histology:= "carcinoid"]
dat_clean[histology == as.character(8490), histology:= "signet.ring"]
dat_clean[histology %in% as.character(8070:8078), histology:= "squamous"]
dat_clean[histology == as.character(8510), histology:= "medullary"]
dat_clean[histology %in% as.character(8040:8043), histology:= "small.cell"]
dat_clean[!(histology %in% c("adenocarcinoma", "papillary", "adenomatous.polyp", "mucinous", "carcinoid", "signet.ring", "squamous", "medullary", "small.cell")), histology:= "other"]
dat_clean[, histology:= factor(histology, levels = c("adenocarcinoma", "papillary", "adenomatous.polyp", "mucinous", "carcinoid", "signet.ring", "squamous", "medullary", "small.cell", "other"))]

#dat_clean[dstage %in% c("0", "Stage 0"), stageb:= "Insitu"]
dat_clean[dstage %in% c("1A", "1A1", "1A2", "1A3", "IA"), stageb:= "Stage.IA"]
dat_clean[dstage %in% c("1B", "IB"), stageb:= "Stage.IB"]
dat_clean[dstage %in% c("21", "2A", "IIA"), stageb:= "Stage.IIA"]
dat_clean[dstage %in% c("22", "2B", "IIB"), stageb:= "Stage.IIB"]
dat_clean[dstage %in% c("2C", "IIC"), stageb:= "Stage.IIC"]
dat_clean[dstage %in% c("31", "3A", "IIIA"), stageb:= "Stage.IIIA"]
dat_clean[dstage %in% c("32", "3B", "IIIB"), stageb:= "Stage.IIIB"]
dat_clean[dstage %in% c("3C", "IIIC"), stageb:= "Stage.IIIC"]
dat_clean[dstage %in% c("4A", "IVA"), stageb:= "Stage.IVA"]
dat_clean[dstage %in% c("4B", "IVB"), stageb:= "Stage.IVB"]
dat_clean[dstage %in% c("4C", "IVC"), stageb:= "Stage.IVC"]
dat_clean[, stageb:= factor(stageb, levels = c("Stage.IA", "Stage.IB", "Stage.IIA", "Stage.IIB", "Stage.IIC", "Stage.IIIA", "Stage.IIIB", "Stage.IIIC", "Stage.IVA", "Stage.IVB", "Stage.IVC"))]

#dat_clean[dstage %in% c("0", "Stage 0"), dstage:= "Insitu"]
dat_clean[dstage %in% c("1", "10", "1A", "1A1", "1A2", "1A3", "1B", "I", "IA", "IB"), dstage:= "Stage.I"]
dat_clean[dstage %in% c("2", "20", "21", "22", "2A", "2B", "2C", "II", "IIA", "IIB", "IIC", "IINOS"), dstage:= "Stage.II"]
dat_clean[dstage %in% c("3", "30", "31", "32", "3A", "3B", "3C", "III", "IIIA", "IIIB", "IIIC", "IIINOS"), dstage:= "Stage.III"]
dat_clean[dstage %in% c("4", "40", "4A", "4B", "4C", "IV", "IVA", "IVB", "IVC", "IVNOS"), dstage:= "Stage.IV"]
dat_clean[, dstage:= factor(dstage, levels = c("Stage.I", "Stage.II", "Stage.III", "Stage.IV"))] #"Insitu", 

#dat_clean[dT %in% c("c0", "p0", "pIS", "T0", "Tis", "Tis(LAMN)"), dT:= "T0"]
dat_clean[dT %in% c("c1", "c1A", "c1B", "c1C", "c1MI", "p1", "p1A", "p1B", "p1C", "p1MI", "T1", "T1a", "T1b", "T1c", "T1mi", "T1mic", "T1NOS"), dnT:= "T1"]
dat_clean[dT %in% c("c2", "c2A", "c2B", "p2", "p2A", "p2B", "T2", "T2a", "T2b", "T2NOS"), dnT:= "T2"]
dat_clean[dT %in% c("c3", "p3", "T3"), dnT:= "T3"]
dat_clean[dT %in% c("cX", "pX", "TX", "TXa", "TXc"), dnT:= "TX"]
dat_clean[dT %in% c("c4", "p4", "T4", "T4NOS"), dnT:= "T4NOS"]
dat_clean[dT %in% c("c4A", "p4A", "T4a"), dnT:= "T4a"]
dat_clean[dT %in% c("c4B", "p4B", "T4b"), dnT:= "T4b"]
dat_clean[dT %in% c("c4C", "p4C", "T4c"), dnT:= "T4c"]
dat_clean[dT %in% c("c4D", "p4D", "T4d"), dnT:= "T4d"]
fwrite(dat_clean[, .(patient.id, dnT)], file = paste0(out_path, "dat_T4b.csv"))

dat_clean[dT %in% c("c1", "c1A", "c1B", "c1C", "c1MI", "p1", "p1A", "p1B", "p1C", "p1MI", "T1", "T1a", "T1b", "T1c", "T1mi", "T1mic", "T1NOS"), dT:= "T1"]
dat_clean[dT %in% c("c2", "c2A", "c2B", "p2", "p2A", "p2B", "T2", "T2a", "T2b", "T2NOS"), dT:= "T2"]
dat_clean[dT %in% c("c3", "p3", "T3"), dT:= "T3"]
dat_clean[dT %in% c("c4", "c4A", "c4B", "c4C", "c4D", "p4", "p4A", "p4B", "p4C", "p4D", "T4", "T4a", "T4b", "T4c", "T4d", "T4NOS"), dT:= "T4"]
dat_clean[dT %in% c("cX", "pX", "TX", "TXa", "TXc"), dT:= "TX"]
dat_clean[, dT:= factor(dT, levels = c("T1", "T2", "T3", "T4", "TX"))] #"T0", 

dat_clean[dN %in% c("c0", "N0", "N0(i-)", "N0(i+)", "N0(mol-)", "N0(mol+)", "p0", "p0I-", "p0I+", "p0M-", "p0M+"), dN:= "N0"]
dat_clean[dN %in% c("c1", "c1A", "c1B", "c1C", "N1", "N1a", "N1b", "N1c", "N1mi", "N1NOS", "N1x", "p1", "p1A", "p1B", "p1C", "p1MI"), dN:= "N1"]
dat_clean[dN %in% c("c2", "c2A", "c2B", "N2", "N2a", "N2b", "N2NOS", "p2", "p2A", "p2B"), dN:= "N2"]
dat_clean[dN %in% c("c3", "c3A", "c3B", "c3C", "N3", "N3a", "N3b", "N3c", "N3NOS", "p3", "p3A", "p3B", "p3C"), dN:= "N3"]
dat_clean[dN %in% c("cX", "NX", "NXr", "NXu", "pX"), dN:= "NX"]
dat_clean[, dN:= factor(dN, levels = c("N0", "N1", "N2", "N3", "NX"))]

dat_clean[dM %in% c("c0", "c0I+", "M0", "M0(i+)"), dM:= "M0"]
dat_clean[dM %in% c("c1", "c1A", "c1B", "M1", "M1a", "M1b", "M1c", "M1NOS", "p1", "p1A", "p1B"), dM:= "M1"]
dat_clean[dM %in% c("MX"), dM:= "MX"]
dat_clean[, dM:= factor(dM, levels = c("M0", "M1", "MX"))]

dat_clean[grade %in% c("1", "Well differentiated; Grade I"), grade:= "G1"]
dat_clean[grade %in% c("2", "Moderately differentiated; Grade II"), grade:= "G2"]
dat_clean[grade %in% c("3", "Poorly differentiated; Grade III"), grade:= "G3"]
dat_clean[grade %in% c("4", "Undifferentiated; anaplastic; Grade IV"), grade:= "G4"]
dat_clean[, grade:= factor(grade, levels = c("G1", "G2", "G3", "G4"))]

dat_clean[radiation %in% c("Beam radiation", "Radiation, NOS  method or source not specified", "Radioactive implants (includes brachytherapy) (1988+)", "Radioisotopes (1988+)", "Combination of beam with implants or isotopes"), radiation:= "Yes"]
dat_clean[radiation %in% c("None/Unknown", "Refused (1988+)", "Recommended, unknown if administered"), radiation:= "No"]
dat_clean[, radiation:= factor(radiation, levels = c("No", "Yes"))]

levels_surg_rad <- c("No radiation and/or cancer-directed surgery", "Radiation prior to surgery", "Intraoperative radiation", "Radiation after surgery", "Radiation before and after surgery", "Surgery both before and after radiation", "Sequence unknown, but both were given")
labels_surg_rad <- c("None", "RadPri", "RadIntra", "RadAfter", "RadBoth", "SurgBoth", "Unknown")
dat_clean[surg_rad %in% c("Intraoperative rad with other rad before/after surgery", "Intraoperative radiation"), surg_rad:= "Intraoperative radiation"]
dat_clean[, surg_rad:= factor(surg_rad, levels = levels_surg_rad, label = labels_surg_rad)]

levels_surg_sys <- c("No systemic therapy and/or surgical procedures", "Systemic therapy before surgery", "Intraoperative systemic therapy", "Systemic therapy after surgery", "Systemic therapy both before and after surgery", "Surgery both before and after systemic therapy", "Sequence unknown")
labels_surg_sys <- c("None", "SysPri", "SysIntra", "SysAfter", "SysBoth", "SurgBoth", "Unknown")
dat_clean[surg_sys %in% c("Intraop systemic rx & oth systemic rx before/after surg", "Intraoperative systemic therapy"), surg_sys:= "Intraoperative systemic therapy"]
dat_clean[, surg_sys:= factor(surg_sys, levels = levels_surg_sys, label = labels_surg_sys)]

levels_reason <- c("Surgery performed", "Recommended but not performed, patient refused", "Not performed, patient died prior to recommended surgery", "Recommended but not performed, unknown reason", "Not recommended", "UnknownSurg")
labels_reason <- c("Performed", "Refused", "Died", "UnknownReas", "NotRecomm", "UnknownSurg")
dat_clean[reason_nosurg %in% c("Not recommended, contraindicated due to other cond; autopsy only (1973-2002)", "Not recommended"), reason_nosurg:= "Not recommended"]
dat_clean[reason_nosurg %in% c("Recommended, unknown if performed", "Unknown; death certificate; or autopsy only (2003+)"), reason_nosurg:= "UnknownSurg"]
dat_clean[, reason_nosurg:= factor(reason_nosurg, levels = levels_reason, labels_reason)]

dat_clean[monDiag2Treat == "Blank(s)", monDiag2Treat:= NA]
dat_clean[, monDiag2Treat:= as.numeric(monDiag2Treat)]

dat_clean[tumor.depo == "81 or more Tumor Deposits", tumor.depo:= 81]
dat_clean[tumor.depo %in% c("Blank(s)", "Not documented/assessed; Indeterminate; No mention in path report; No resection", "Tumor Deposits identified, number unknown"), tumor.depo:= NA]
dat_clean[tumor.depo == "No tumor deposits", tumor.depo:= 0]
dat_clean[, tumor.depo:= as.numeric(tumor.depo)]


cod_tab <- fread("/home/zhouwen/project/population/codes/Revise/cod_type.csv")
for (itype in unique(cod_tab$Type)) {
  dat_clean[cod.site %in% cod_tab[Type == itype]$cod, cod:= itype]
}
dat_clean[, cod:= factor(cod, levels = c("Alive", "Benign.or.unknown", "Other.causes", "External.causes", "Infectious", "Repiratory", "Cardiovascular", "Diabetes.Mellitus", "Alzheimers", "Renal", "GI.and.liver", "Other.cancers", "CRC"))]

dat_clean[, caus.death:= factor(ifelse(cod == "Alive", "Alive", ifelse(cod == "CRC", "CRC", "nonCRC")), levels = c("Alive", "CRC", "nonCRC"))]

dat_clean[surv.month == "Unknown", surv.month:= NA]
dat_clean[, surv.month:= as.numeric(surv.month)]
dat_clean[, surv.death:= ifelse(vital == "Dead", 1, 0)]

save(dat_clean, file = paste0(out_path, "dat_cleaned.RData"))
