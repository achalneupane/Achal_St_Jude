rm(list=ls())

library(survival)
library(splines)
library(ggplot2)
library(grid)
library(ggpubr)

# for cumulative incidence plots
library(survminer)
library(cmprsk)
library(dplyr)
library(utile.visuals) # geom_stepconfint

options(digits=10)

maindir = 'R:/Biostatistics/Biostatistics/cim1/GxTPRS_snpnet/'

# Sex, age at dx; dx year (every 5 years)
# max RT dose to any body region (None, cut every 10Gy with right closed intervals; >=50 highest)
# alkylating agents: CED dose (None; 1-3999; 4000-7999; >=8000)
# Anthracycline dose (None; 0-100; 101-300; >300)
# Epipodophyllotoxin dose (None; 1-1000; 1001-4000; >4000)
# Add ancestry PCs: note that PCs in the datasets are based on cohort-specific PCA
# Note that "PC1...PCX" = PCs from cohort- and batch-specific PCA (e.g., CCSS Original only)

##############################################################################
# Get your PCs (load your own)

load(paste0(maindir,'BCC/EURpcs.RData'))
str(newpcs) 
# first 10 EUR PCs computed for CCSS Ori and Exp combined and then for SJLIFE B1 and B2 combined (example)

##############################################################################
# PRS load

load(paste0(maindir,'BCC/samplePRS.RData'))
str(gpprs.comb) # example PRSs computed with my harmonized WGS/array dataset; labeled with PGS Catalog numbers

# Standard normal PRS

prscomb = merge(newpcs,gpprs.comb,by='studyid') # merge PCs and PRSs
prscomb$sn_anycancer = as.vector(scale(prscomb[,'SCORESUM_PGS000356_combined'],center=TRUE,scale=TRUE))
prscomb$sn_bcc = as.vector(scale(prscomb[,'SCORESUM_PGS000459_combined'],center=TRUE,scale=TRUE))
prscomb$sn_bcc2 = as.vector(scale(prscomb[,'SCORESUM_PGS000119_combined'],center=TRUE,scale=TRUE))
prscomb$sn_scc = as.vector(scale(prscomb[,'SCORESUM_PGS000120_combined'],center=TRUE,scale=TRUE))

# Categorical quintiled PRS

prscomb$sn_anycancer_cat5 = cut(prscomb$sn_anycancer,quantile(prscomb$sn_anycancer,seq(0,1,0.2)))
prscomb$sn_bcc_cat5 = cut(prscomb$sn_bcc,quantile(prscomb$sn_bcc,seq(0,1,0.2)))
prscomb$sn_bcc2_cat5 = cut(prscomb$sn_bcc2,quantile(prscomb$sn_bcc2,seq(0,1,0.2)))
prscomb$sn_scc_cat5 = cut(prscomb$sn_scc,quantile(prscomb$sn_scc,seq(0,1,0.2)))

str(prscomb)

################# BCC PHENOTYPE DATA LOAD ######################

load(paste0('R:/Biostatistics/Biostatistics/cim1/SMN prediction/bcc_covs_Oct4_2022.RData'))

# Long-format multiple BCC datasets; note that overlaps between SJLIFE and CCSS have been removed from CCSS
# Consistent with our NMSC prediction analyses where we are keeping SJLIFE sample size as big as possible
# Long format data is used for MULTIPLE BCC analysis!
# Note that counting process stopped after experiencing a certain number of BCCs (I think 5 BCCs)

dim(ccss.df)
length(unique(ccss.df$studyid))
table(ccss.df$event) # 2835 events in CCSS
dim(sjl.df)
length(unique(sjl.df$studyid))
table(sjl.df$event) # 630 events in SJLIFE

# FLATTEN DATA
# Age as time scale, but short format data (FIRST BCC analysis!)
# Flatten so that we have start age and end age, but not divided into 1-year time segments

ccss.dfF = ccss.df[ccss.df$evt1==1,] # filter first events
table(ccss.dfF$event) # 1170 first events
sjl.dfF = sjl.df[sjl.df$evt1==1,] # filter first events
table(sjl.dfF$event) # 256 first events

# CCSS

temp1 = unique(ccss.dfF[,which(!(colnames(ccss.dfF) %in% c('start','end','event')))]) 
temp2 = do.call(rbind,lapply(split(ccss.dfF,ccss.dfF$studyid),function(x){
	start = min(x$start)
	end = max(x$end)
	event = sum(x$event)
	return(c(unique(x$studyid),start,end,event))
}))
temp2 = as.data.frame(temp2,stringsAsFactors = F)
colnames(temp2) = c('studyid','start','end','event')
str(temp2)
table(temp2$event) # Check that # events didn't change
summary(temp2$end - temp2$start) # Check that this is non-zero (end age>start age)
ccss.dfF = merge(temp2,temp1,by='studyid')
dim(ccss.dfF) # 23219 pts in CCSS data
table(ccss.dfF$event) # 1170 first events

# SJLIFE

temp1 = unique(sjl.dfF[,which(!(colnames(sjl.dfF) %in% c('start','end','event')))]) 
temp2 = do.call(rbind,lapply(split(sjl.dfF,sjl.dfF$studyid),function(x){
	start = min(x$start)
	end = max(x$end)
	event = sum(x$event)
	return(c(unique(x$studyid),start,end,event))
}))
temp2 = as.data.frame(temp2,stringsAsFactors = F)
colnames(temp2) = c('studyid','start','end','event')
str(temp2)
temp2[,2:4] = apply(temp2[,2:4],2,as.numeric)
table(temp2$event) # Check that # events didn't change
summary(temp2$end - temp2$start) # Check that this is non-zero (end age>start age)
sjl.dfF = merge(temp2,temp1,by='studyid')
dim(sjl.dfF) # 5531 pts in the SJLIFE data
table(sjl.dfF$event) # 256 events

##############################################################################
# Merge PRS and pheno data

# Multiple events

ccss.dfg = ccss.df[ccss.df$studyid %in% prscomb$studyid,]
ccss.dfg = merge(ccss.dfg,prscomb,by='studyid',all.x=T)
length(unique(ccss.dfg$studyid)) # 6634 pts
sum(ccss.dfg$event) # 1488 events

sjl.dfg = sjl.df[sjl.df$studyid %in% prscomb$studyid,]
sjl.dfg = merge(sjl.dfg,prscomb,by='studyid',all.x=T)
length(unique(sjl.dfg$studyid)) # 3261 pts
sum(sjl.dfg$event) # 505 events

# First events

ccss.dfFg = ccss.dfF[ccss.dfF$studyid %in% prscomb$studyid,]
ccss.dfFg = merge(ccss.dfFg,prscomb,by='studyid',all.x=T)
length(unique(ccss.dfFg$studyid)) # 6634 pts
sum(ccss.dfFg$event) # 614 events

sjl.dfFg = sjl.dfF[sjl.dfF$studyid %in% prscomb$studyid,]
sjl.dfFg = merge(sjl.dfFg,prscomb,by='studyid',all.x=T)
length(unique(sjl.dfFg$studyid)) # 3261 pts
sum(sjl.dfFg$event) # 207 events

##############################################################################

# Sample Cox PH models: okay to use for multiple events and first event analysis

# example; note that since age as time scale is used so you do not need to adjust for age

df = sjl.dfFg # SJLIFE data, first event
mod = coxph(Surv(start,end,event) ~ agedx + gender + 
	maxrtdose + # max RT dose to any of 7 body regions
	aadose + anthdose + epitxndose + cisplateqdose + # doses for alkylators, anthracyclines, epipodophyllotoxins, platinums
	COMBSTPC1 + COMBSTPC2 + COMBSTPC3 + COMBSTPC4 + COMBSTPC5 + # first 5 EUR PCs
	sn_bcc2_cat5, # categorical PRS
	weights = fweight,  # sampling weights, would be used in CCSS but okay in SJLIFE because all equal to 1
	data = df)

#Call:
#coxph(formula = Surv(start, end, event) ~ agedx + gender + maxrtdose + 
#    aadose + anthdose + epitxndose + cisplateqdose + COMBSTPC1 + 
#    COMBSTPC2 + COMBSTPC3 + COMBSTPC4 + COMBSTPC5 + sn_bcc2_cat5, 
#    data = df, weights = fweight)
#
#  n= 3214, number of events= 202 
#   (47 observations deleted due to missingness)
#
#                                     coef     exp(coef)      se(coef)        z   Pr(>|z|)    
#agedx                       -3.075081e-02  9.697172e-01  1.329020e-02 -2.31380 0.02067891 *  
#genderMale                  -2.371179e-01  7.888983e-01  1.421064e-01 -1.66859 0.09519790 .  
#maxrtdose                    2.679820e-04  1.000268e+00  3.975271e-05  6.74122 1.5706e-11 ***
#aadose                       9.749086e-06  1.000010e+00  8.218713e-06  1.18621 0.23554099    
#anthdose                    -6.029794e-04  9.993972e-01  5.773474e-04 -1.04440 0.29630212    
#epitxndose                  -4.819572e-05  9.999518e-01  2.702505e-05 -1.78337 0.07452572 .  
#cisplateqdose               -1.271341e-03  9.987295e-01  8.009553e-04 -1.58728 0.11244904    
#COMBSTPC1                    8.115402e-01  2.251373e+00  4.750722e+00  0.17082 0.86436170    
#COMBSTPC2                   -4.147969e+00  1.579646e-02  4.755947e+00 -0.87216 0.38311847    
#COMBSTPC3                    6.570263e+00  7.135573e+02  6.227825e+00  1.05499 0.29143216    
#COMBSTPC4                    2.082767e+00  8.026649e+00  5.407784e+00  0.38514 0.70013187    
#COMBSTPC5                   -6.239562e+00  1.950709e-03  4.186498e+00 -1.49040 0.13611874    
#sn_bcc2_cat5(-0.846,-0.272]  2.407734e-01  1.272233e+00  2.756898e-01  0.87335 0.38247291    
#sn_bcc2_cat5(-0.272,0.241]   2.996533e-01  1.349391e+00  2.836837e-01  1.05629 0.29083417    
#sn_bcc2_cat5(0.241,0.835]    8.603794e-01  2.364057e+00  2.526944e-01  3.40482 0.00066207 ***
#sn_bcc2_cat5(0.835,3.64]     9.642111e-01  2.622718e+00  2.508433e-01  3.84388 0.00012110 ***

##############################################################################
# CIF plots (Mar 2023, top quintile vs bottom quintile)

# get deaths, necessary for cumulative incidence accounting for death competing risk calculation
load(paste0(maindir,'BCC/death_data.RData'))
str(sjldeath) # deaths in SJLIFE
str(ccssdeath) # deaths in CCSS

# combine and calculate age at death
dcomb = rbind(sjldeath[,c('studyid','DOB','DEATH_LASTDT','death')],
	ccssdeath[,c('studyid','DOB','DEATH_LASTDT','death')])
dcomb$agedeath = as.numeric(difftime(dcomb$DEATH_LASTDT, dcomb$DOB, units='days')/365.25)
summary(dcomb$agedeath)

mydata = sjl.dfFg # whatever dataset you want, this is SJLIFE first events

dcomb = dcomb[dcomb$studyid %in% mydata$studyid,] 
str(dcomb)

comb.dfFg2 = merge(mydata,dcomb[,c('studyid','death','agedeath')],by='studyid')
comb.dfFg2$fstatus = comb.dfFg2$event # keep all events
comb.dfFg2$fstatus[comb.dfFg2$event==1] = 2 # code events as 2
comb.dfFg2$fstatus[comb.dfFg2$death==1 & comb.dfFg2$event==0] = 1 # pts who die can't experience event, code as 1
table(comb.dfFg2$fstatus,comb.dfFg2$event)

# all data set up to remove pt from at risk pool at death, so age at death should be same as end age
#checkid = sample(unique(comb.dfFg2$studyid[comb.dfFg2$fstatus==1]),1)
#comb.dfFg2[comb.dfFg2$studyid==checkid,c('studyid','end','event','death','agedeath')]

####################################
# Manual competing risk CIF plot: run crplot and tabplot functions

crplot = function(sfit,lowgp,higp,t30,t30.ll,t30.ul,pval,myymax,mytitle) {

	# Plots: Start at time=0, value=0; cannot end in a vertical line
	# Last CI estimate stands at time X until time is equal to nearest tens place (e.g., if f/u ends at yr 27, keep that CI until yr 30)
	# sfit = survfit object
	# lowgp = 'PRS Q1', higp = 'PRS Q5'
	# myymax = 0.5
	# mytitle = 'BCC, SJLIFE'

	df = data.frame(time=c(summary(sfit)$time), 
		etype = rep('BCC',length(summary(sfit)$time)),
		est=c(summary(sfit)$pstate[,3]),
		lower=c(summary(sfit)$lower[,3]),
		upper=c(summary(sfit)$upper[,3]),
		Strata=summary(sfit)$strata)
	levels(df$Strata)[grep('=0',levels(df$Strata))] <- lowgp
	levels(df$Strata)[grep('=1',levels(df$Strata))] <- higp
	df$etype = factor(df$etype,levels=c('Death','BCC'))	

	maxT.prs1 = 30 # cutoff f/u time
	maxT.prs0 = 30 # cutoff f/u time
	maxT.prs1.actual = max(df$time[df$Strata==higp & df$time<=30]) # to match last CI estimate
	maxT.prs0.actual = max(df$time[df$Strata==lowgp & df$time<=30]) # to match last CI estimate

	df = df[df$time<=30 & df$etype=='BCC',]
		
	# Extend values to t=30 for survivors

	addon2 = c(maxT.prs0,levels(df$etype)[2],
		df$est[df$Strata==lowgp & df$etype==levels(df$etype)[2] & df$time==maxT.prs0.actual],
		df$lower[df$Strata==lowgp & df$etype==levels(df$etype)[2] & df$time==maxT.prs0.actual],
		df$upper[df$Strata==lowgp & df$etype==levels(df$etype)[2] & df$time==maxT.prs0.actual],
		lowgp)
	addon4 = c(maxT.prs1,levels(df$etype)[2],
		df$est[df$Strata==higp & df$etype==levels(df$etype)[2] & df$time==maxT.prs1.actual],
		df$lower[df$Strata==higp & df$etype==levels(df$etype)[2] & df$time==maxT.prs1.actual],
		df$upper[df$Strata==higp & df$etype==levels(df$etype)[2] & df$time==maxT.prs1.actual],
		higp)
	df = rbind(df,addon2,addon4)

	# Add time=0
	df = rbind(df,
		c(0,levels(df$etype)[2],0,0,0,lowgp),
		c(0,levels(df$etype)[2],0,0,0,higp))	

	df[,c('time','est','lower','upper')] = apply(df[,c('time','est','lower','upper')],2,as.numeric)
	df = df[order(df$Strata,df$time),]
	
	# Manual CR plot

	p2 = ggplot(df[df$etype=='BCC',], aes(x = time, y = est, color = Strata)) +
 		geom_step(aes(color = Strata),size=1.05)  +
		geom_stepconfint(aes(ymin = lower, ymax = upper, fill = Strata), alpha = 0.2, colour = NA) +
		ylim(0,myymax) +
		scale_color_manual(values=c('black','dodgerblue'),aesthetics = c('colour', 'fill')) +
		coord_cartesian(xlim = c(0, 30)) +
		theme_classic() +
		theme(plot.title = element_text(size = 12),
			#legend.title = element_blank(),
			legend.position = "top",
			axis.text.x = element_text(size=10,color='black'),
			axis.text.y = element_text(size=10,color='black'),
			axis.title.x = element_text(size=10, color = "black", family = "sans"),
      		axis.title.y = element_text(size=10, color = "black", family = "sans")) +
			labs(x = paste0("Years since 5 years post-diagnosis"), 
				y = "BCC probability", 
				title = mytitle)

	p2 = p2 + theme(plot.margin=margin(0,5,1,5))
	p2 = p2 +
		ggplot2::annotate("text", x = 2, y = 0.47,
			label = paste0(lowgp), size = 3.5) +
		ggplot2::annotate("text", x = 2, y = 0.43,
			label = paste0('30-year (95% CI): ',t30[1],' (',t30.ll[1],' to ',t30.ul[1],')'), size = 3.5, hjust=0) +
		ggplot2::annotate("text", x = 2, y = 0.39,
			label = paste0(higp), size = 3.5) +
		ggplot2::annotate("text", x = 2, y = 0.35,
			label = paste0('30-year (95% CI): ',t30[2],' (',t30.ll[2],' to ',t30.ul[2],')'), size = 3.5, hjust=0) +
		ggplot2::annotate("text", x = 2, y = 0.15,
			label = paste0('p = ',pval), size = 3.5, hjust=0)

	return(p2)

}

tabplot = function(sfit,data,lowgp,higp) {
	
	tab = ggsurvplot(sfit, data = data, 
		xlim=c(0,30),
		pval=F,
		fun=function(x)(1-x),
		conf.int=T,
		pval.size=3.5,
		font.main=12,font.legend=10,font.x=10,font.y=10,font.tickslab=10,
		risk.table=T,
		risk.table.title='',risk.table.fontsize=3.5,
		legend='top',,
		legend.labs=c(lowgp,higp),
		palette=c('black','dodgerblue2'),
		ylab='BCC (probability)')
	tab$table = tab$table + labs(x='') + labs(y='Number at risk') 
	tab$table$theme$text$size = 10
	tab$table$theme$axis.text.y$size = 10
	tab$table = tab$table + theme(plot.margin=margin(-20,5,1,5), axis.text.x = element_text(size=10,color='black'))
	tab$table = tab$table + theme(axis.title.x = element_text(size=10, color = "black", family='sans'))
	tab$table = tab$table + theme(axis.title.y = element_text(size=10, color = "black", family='sans'))

	return(tab$table)

}

# Check against cuminc with no weights (to check the competing risk CI is computed correctly with survfit)
#data$fstatus3 = NA
#data$fstatus3[data$fstatus2=='Censored'] = 0
#data$fstatus3[data$fstatus2=='Death'] = 1
#data$fstatus3[data$fstatus2=='SMN'] = 2
#sfit2 = cuminc(data$futime, data$fstatus3, group=data$prs5, rho=0, cencode=0)
#sfit2
#sfit.nowt = survfit(Surv(data$futime, data$fstatus2) ~ data$prs5) # no weights competing risk with survfit
#print(summary(sfit.nowt,time=30))

############################################################################
# Compute cumulative incidence at 30y post 5y dx of childhood cancer dx accounting for competing risk of death
# Deals with CCSS sampling weights
# Generates a TWO-GROUP curve plot (binary risk factor, e.g., non-carrier vs. carrier) for BCC

# PRS type? 
# If rare variant, specify the variable name (MUST CODE VBLE AS 0/1 INDICATOR!)
prstype = 'sn_bcc2_cat5'

# What dataset?
data = na.omit(comb.dfFg2[,c('studyid','start','end','event','fstatus','fweight',
	'maxrtdose','aadose','anthdose','epitxndose','cisplateqdose',
	prstype)])
data$futime = data$end-data$start
summary(data$futime) # check non-zero values

# Set up data to be BINARY top quintile vs bottom (if binary vble, skip this)
data$prs5 = ifelse(data[,prstype]==levels(data[,prstype])[5],1,0)
data = data[data[,prstype]==levels(data[,prstype])[5] | data[,prstype]==levels(data[,prstype])[1],]
dim(data)

# Compute cumulative incidence
data$fstatus2 = factor(data$fstatus, 0:2, labels=c("Censored", "Death", "BCC")) # set up factor vble with death
sfit = survfit(Surv(futime, fstatus2) ~ prs5, weights = fweight, data = data) # compute cumulative incidence with competing risk of death
(t30 = formatC(summary(sfit,time=30)$pstate[,3], digits = 2, format = "f")) # cumulative inc at t=30
(t30.ll = formatC(summary(sfit,time=30)$lower[,3], digits = 2, format = "f")) # 95% CI LL
(t30.ul = formatC(summary(sfit,time=30)$upper[,3], digits = 2, format = "f")) # 95% CI UL

# Check max follow-up time is past 30 for both risk levels (edit "data$prs5==0" and "data$prs5==1" to your vble of interest)
print(max(data[data$prs5==0,'futime']))
print(max(data[data$prs5==1,'futime']))

# to get pvalue
data$fstatus3 = NA
data$fstatus3[data$fstatus2=='Censored'] = 0
data$fstatus3[data$fstatus2=='Death'] = 1
data$fstatus3[data$fstatus2=='SMN'] = 2
sfit2 = cuminc(data$futime, data$fstatus3, group=data$prs5, rho=0, cencode=0)
sfit2
## EDIT pvalue from sfit2
pval = 0.25 

# cumulative incidence plot
thisymax = 0.5 # y axis limit upper bound
thistitle = 'BCC, SJLIFE' # plot title
mylow = 'PRS Q1' # label for the lower risk group
myhigh = 'PRS Q5' # label for the higher risk group
myplot = crplot(sfit,mylow,myhigh,t30,t30.ll,t30.ul,pval,thisymax,thistitle)

# Risk tables made separately b/c sample-weighted counts are provided
sfit.tab = survfit(Surv(data$futime,data$event)~data$prs5)
mytab = tabplot(sfit.tab,data,mylow,myhigh)

cifplot = ggarrange(myplot,mytab,heights=c(1,0.2),ncol=1,nrow=2,align="v")

#save(data,sfit,sfit2,t30,t30.ll,t30.ul,myplot,mytab,file=paste0(maindir,'BCC/AN_CIF.RData'))


