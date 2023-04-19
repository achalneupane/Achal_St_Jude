## Extract echo measures
rm(list=ls())
# setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//')
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/")
## Phenotype data for GWAS
pheno_gwas = read.table('pheno/sjlife_ttn_bag3.pheno', header = TRUE)

## Echo measures with data from most recent visit
echo_data = read.table('../echo_data/SJLIFE_echo_measures_most_recent_visit.txt', header = TRUE)

## Merge both data
pheno_final = merge(pheno_gwas, echo_data, by.x = 'FID', by.y = 'sjlid')

## Genotypes of TTN and BAG3 top SNPs
bag3 = read.table('sjlife_chr10_119670121_T_C.raw', header = TRUE)
bag3[c(2:6)] = NULL
ttn = read.table('sjlife_chr2_178562809_T_C.raw', header = TRUE)
ttn[c(2:6)] = NULL
geno = merge(bag3, ttn, by="FID")

## Merge phenotype and genotype data
dat_final = merge(pheno_final, geno, by = 'FID')

## Analyze both variants wrt each echo measures
# BAG3 SNP
res_bag3 = NULL
echo_phenotypes = colnames(echo_data)[-c(1:4)]
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~chr10.119670121.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~chr10.119670121.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  res_bag3 = rbind(res_bag3, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
res_bag3$p_BH = p.adjust(res_bag3$pval, method = "BH")

# TTN SNP
res_ttn = NULL
for (p in 1:length(echo_phenotypes)){
  dat_eff = subset(dat_final, !is.na(dat_final[echo_phenotypes[p]]))
  # Remove values with >3 sd away from the absolute value
  dat_eff = subset(dat_eff, abs(scale(dat_eff[,echo_phenotypes[p]])) <= 6)
  # Adjust for BSA; two (LVMass2D_Index and LA_Volume_Index) are already done
  if (!(echo_phenotypes[p] == "LVMass2D_Index" | echo_phenotypes[p] == "LA_Volume_Index"| 
        echo_phenotypes[p] == "LV_Ejection_Fraction_3D" | echo_phenotypes[p] == "LV_Cardiac_Output_3D" | echo_phenotypes[p] == "LV_GLPS_AVG")){
    fit = lm(scale(dat_eff[,echo_phenotypes[p]]/BodySurfaceArea)~chr2.178562809.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  } else {
    fit = lm(scale(dat_eff[,echo_phenotypes[p]])~chr2.178562809.T.C_C+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dat_eff)
  }
  beta = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  pval = summary(fit)$coef[2,4]
  res_ttn = rbind(res_ttn, data.frame(phenotype=echo_phenotypes[p], beta, se, pval, n=nrow(dat_eff)))
}
res_ttn$p_BH = p.adjust(res_ttn$pval, method = "BH")

print(res_bag3); print(res_ttn)

#################
## Forest plot ##
#################
res_bag3$phenotype <- gsub("_", " ", res_bag3$phenotype)
res_ttn$phenotype <- gsub("_", " ", res_ttn$phenotype)

## BAG3
summary_table <- res_bag3[c("phenotype", "beta", "se", "p_BH", "n")]
summary_table$P <- summary_table[,"p_BH"]
summary_table$OR <- exp(summary_table[,"beta"])


# # For Beta
summary_table$lower <- round(summary_table[,"beta"] -1.96*summary_table[,"se"],2)
summary_table$upper <- round(summary_table[,"beta"] +1.96*summary_table[,"se"], 2)

# # For OR
# summary_table$lower <- round(exp(summary_table[,"beta"]) -1.96*summary_table[,"se"],2)
# summary_table$upper <- round(exp(summary_table[,"beta"]) +1.96*summary_table[,"se"], 2)



library(ggplot2)
library(meta)

# # Create a new data frame with non-missing n values
summary_table_n_bag3 <- subset(summary_table, !is.na(n))
# 
# # Plot the data with labels for n
# ggplot(data = summary_table_n_bag3, aes(x = beta, y = phenotype)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_text(aes(x = -0.4, y = phenotype, label = paste0("n = ", n)), hjust = 0, size = 3) +
#   scale_x_continuous(limits = c(-0.4, 0.25), breaks = seq(-0.4, 0.25, 0.1)) +
#   xlab("Estimates (95% CI)") +
#   ylab("Echo measures") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())


## TTN
summary_table <- res_ttn[c("phenotype", "beta", "se", "p_BH", "n")]
summary_table$P <- summary_table[,"p_BH"]
summary_table$OR <- exp(summary_table[,"beta"])

# # For Beta
summary_table$lower <- round(summary_table[,"beta"] -1.96*summary_table[,"se"],2)
summary_table$upper <- round(summary_table[,"beta"] +1.96*summary_table[,"se"], 2)

# # For OR
# summary_table$lower <- round(exp(summary_table[,"beta"]) -1.96*summary_table[,"se"],2)
# summary_table$upper <- round(exp(summary_table[,"beta"]) +1.96*summary_table[,"se"], 2)


library(ggplot2)
library(meta)

# Create a new data frame with non-missing n values
summary_table_n_ttn <- subset(summary_table, !is.na(n))

# # Plot the data with labels for n
# p2 <- ggplot(data = summary_table_n_ttn, aes(x = beta, y = phenotype)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   geom_text(aes(x = -0.4, y = phenotype, label = paste0("n = ", n)), hjust = 0, size = 3) +
#   scale_x_continuous(limits = c(-0.4, 0.25), breaks = seq(-0.4, 0.25, 0.2)) +
#   xlab("Estimates (95% CI)") +
#   ylab("Echo measures") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())


## Combined plot
# library(gridExtra)
# grid.arrange(p1, p2, ncol = 2)

# Add a 'variable' column to each data frame
summary_table_n_bag3$variable <- "BAG3"
summary_table_n_ttn$variable <- "TTN"

# Combine the data frames
summary_table <- rbind(summary_table_n_bag3, summary_table_n_ttn)

# # Create the combined plot with facet_wrap
# combined_plot <- ggplot(summary_table, aes(x = beta, y = phenotype)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_text(aes(x = -0.4, y = phenotype, label = paste0("n = ", n)), hjust = 0, size = 3) +
#   scale_x_continuous(limits = c(-0.4, 0.25), breaks = seq(-0.4, 0.25, 0.2)) +
#   xlab("Estimates (95% CI)") +
#   ylab("Echo measures") +
#   facet_wrap(~ variable, ncol = 2) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank())

# # N in phenotype
# summary_table$phenotype <- as.factor(summary_table$phenotype)
# combined_plot <- ggplot(summary_table, aes(x = beta, y = paste0(phenotype, " (n = ", n, ")"))) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_text(aes(x = -0.4, y = paste0(phenotype, " (n = ", n, ")"), label = ""), hjust = 0, size = 3) +
#   scale_x_continuous(limits = c(-0.4, 0.25), breaks = seq(-0.4, 0.25, 0.2)) +
#   xlab("Estimates (95% CI)") +
#   ylab("Echo measures") +
#   facet_wrap(~ variable, ncol = 2) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank())

summary_table$phenotype <- as.factor(summary_table$phenotype)
combined_plot <- ggplot(summary_table, aes(x = beta, y = phenotype)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-0.4, 0.25), breaks = seq(-0.4, 0.25, 0.2)) +
  xlab("Estimates (95% CI)") +
  ylab("Echo measures") +
  facet_wrap(~ variable, ncol = 2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

combined_plot
ggsave("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Figures/echo_measures.tiff", combined_plot, dpi = 600, width = 5, height = 3, units = "in")




##################################################################################################################
# ############################################## 
# #################
# ## Forest plot ##
# #################
# res_bag3$phenotype <- gsub("_", " ", res_bag3$phenotype)
# res_ttn$phenotype <- gsub("_", " ", res_ttn$phenotype)
# 
# ## BAG3
# summary_table <- res_bag3[c("phenotype", "beta", "se", "p_BH", "n")]
# summary_table$P <- summary_table[,"p_BH"]
# summary_table$OR <- exp(summary_table[,"beta"])
# 
# # For OR
# summary_table$lower <- round(exp(summary_table[,"beta"]) -1.96*summary_table[,"se"],2)
# summary_table$upper <- round(exp(summary_table[,"beta"]) +1.96*summary_table[,"se"], 2)
# 
# 
# 
# library(ggplot2)
# library(meta)
# 
# # Create a new data frame with non-missing n values
# summary_table_n_bag3 <- subset(summary_table, !is.na(n))
# 
# # Plot the data with labels for n
# p1 <- ggplot(data = summary_table_n_bag3, aes(x = OR, y = phenotype)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   geom_text(aes(x = 0.5, y = phenotype, label = paste0("n = ", n)), hjust = 0, size = 3) +
#   scale_x_continuous(limits = c(0.5, 1.5), breaks = seq(0.5, 1.5, 0.3)) +
#   xlab("Odds ratio (95% CI)") +
#   ylab("Echo measures") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())
# 
# 
# ## TTN
# summary_table <- res_ttn[c("phenotype", "beta", "se", "p_BH", "n")]
# summary_table$P <- summary_table[,"p_BH"]
# summary_table$OR <- exp(summary_table[,"beta"])
# 
# summary_table$lower <- round(exp(summary_table[,"beta"]) -1.96*summary_table[,"se"],2)
# summary_table$upper <- round(exp(summary_table[,"beta"]) +1.96*summary_table[,"se"], 2)
# 
# 
# library(ggplot2)
# library(meta)
# 
# # Create a new data frame with non-missing n values
# summary_table_n_ttn <- subset(summary_table, !is.na(n))
# 
# # Plot the data with labels for n
# p2 <- ggplot(data = summary_table_n_ttn, aes(x = OR, y = phenotype)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   geom_text(aes(x = 0.5, y = phenotype, label = paste0("n = ", n)), hjust = 0, size = 3) +
#   scale_x_continuous(limits = c(0.5, 1.5), breaks = seq(0.5, 1.5, 0.3)) +
#   xlab("Odds ratio (95% CI)") +
#   ylab("Echo measures") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())
# 
# 
# ## Combined plot
# # library(gridExtra)
# # grid.arrange(p1, p2, ncol = 2)
# 
# # Add a 'variable' column to each data frame
# summary_table_n_bag3$variable <- "BAG3"
# summary_table_n_ttn$variable <- "TTN"
# 
# # Combine the data frames
# summary_table <- rbind(summary_table_n_bag3, summary_table_n_ttn)
# 
# # Create the combined plot with facet_wrap
# combined_plot <- ggplot(summary_table, aes(x = OR, y = phenotype)) +
#   geom_point(size = 3) +
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   geom_text(aes(x = 0.5, y = phenotype, label = paste0("n = ", n)), hjust = 0, size = 3) +
#   scale_x_continuous(limits = c(0.5, 1.5), breaks = seq(0.5, 1.5, 0.3)) +
#   xlab("Odds ratio (95% CI)") +
#   ylab("Echo measures") +
#   facet_wrap(~ variable, ncol = 2) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank())
# 
# combined_plot
# 
# library(sjPlot)
# plot_model(fit)

