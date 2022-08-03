
library(data.table)
library(susieR)
set.seed(1)
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/")

dat1 <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/data/small_data_11.rds")
dat3 <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main//data/small_data_11_sim_gaussian_pve_n_8_get_sumstats_n_1.rds")
maf  <- dat1$maf$in_sample
bhat <- dat3$sumstats$bhat
shat <- dat3$sumstats$shat
z    <- bhat/shat

dat2 <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/data/small_data_11_sim_gaussian_pve_n_8.rds")
b <- drop(dat2$meta$true_coef)
vars <- which(b != 0)
vars

ldinfile <- "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/data/small_data_11_sim_gaussian_pve_n_8_get_sumstats_n_1.ld_sample_n_file.in_n.ld"
Rin <- as.matrix(fread(ldinfile))
fit1 <- susie_rss(z,Rin,n = 800,min_abs_corr = 0.1,refine = FALSE,
                  verbose = TRUE)

print(fit1$sets[c("cs","purity")])


par(mar = c(4,4,1,1))
cs1 <- fit1$sets$cs$L1
plot(1:1001,fit1$pip,pch = 20,cex = 0.8,ylim = c(0,0.1),
     xlab = "SNP",ylab = "susie PIP")
points(cs1,fit1$pip[cs1],pch = 1,cex = 1,col = "cyan")
points(vars,fit1$pip[vars],pch = 2,cex = 0.8,col = "tomato")


run_finemap <- function (bhat, shat, maf, prefix) {
  p <- length(b)
  dat <- data.frame(rsid       = 1:p,
                    chromosome = rep(1,p),
                    position   = rep(1,p),
                    allele1    = rep("A",p),
                    allele2    = rep("C",p),
                    maf        = round(maf,digits = 6),
                    beta       = round(bhat,digits = 6),
                    se         = round(shat,digits = 6))
  outfile <- paste(prefix,"z",sep = ".")
  masterfile <- paste(prefix,"master",sep = ".")
  write.table(dat,outfile,quote = FALSE,col.names = TRUE,row.names = FALSE)
  out <- system(paste("./finemap_v1.4.1_x86_64 --sss --log --n-causal-snps 5", "--in-files",masterfile),intern = TRUE)
  out <- out[which(!grepl("Reading",out,fixed = TRUE))]
  out <- out[which(!grepl("Computing",out,fixed = TRUE))]
  out <- out[which(!grepl("evaluated",out,fixed = TRUE))]
  cat(out,sep = "\n")
  return(invisible(out))
}


setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/output")
run_finemap(bhat,shat,maf,prefix = "sim1")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/analysis/")


par(mar = c(4,4,1,1))
finemap <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/output/sim1.cred2",header = TRUE)
pip   <- rep(0,1001)
cs1   <- finemap$cred1
rows1 <- which(!is.na(cs1))
cs1   <- cs1[rows1]
pip[cs1] <- pip[cs1] + finemap$prob1[rows1]
plot(1:1001,pip,pch = 20,cex = 0.8,ylim = c(0,0.1),
     xlab = "SNP",ylab = "finemap PIP")
points(cs1,pip[cs1],pch = 1,cex = 1,col = "cyan")
points(vars,pip[vars],pch = 2,cex = 0.8,col = "tomato")


ldoutfile <- "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/data/small_data_11.ld_refout_file.refout.ld"
Rout <- as.matrix(fread(ldoutfile))
fit2 <- susie_rss(z,Rout,n = 800,min_abs_corr = 0.1,refine = FALSE,
                 verbose = TRUE)

print(fit2$sets[c("cs","purity")])


all(fit1$sets$cs$L1 == fit2$sets$cs$L1)

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/output")
run_finemap(bhat,shat,maf,prefix = "sim2")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/analysis")


cat(readLines("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/helsinki/finemap-main/output/sim2.cred5"),sep = "\n")


[dsc]: https://github.com/stephenslab/dsc_susierss
[finemap]: http://www.christianbenner.com
[susie]: https://github.com/stephenslab/susieR