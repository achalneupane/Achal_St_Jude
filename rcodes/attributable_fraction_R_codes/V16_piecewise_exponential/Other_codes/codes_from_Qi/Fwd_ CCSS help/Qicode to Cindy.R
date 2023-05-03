
#### Explain what we need to do: see file "Explain counting process with multiple events.docx"


rm(list=ls())
load("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/V16_piecewise_exponential/codes_from_Qi/Fwd_ CCSS help/STROKE_counting_process_data_to_check_for_Qi.RData")

library("survival")

stroke_multi_count = NULL
dat <- stroke_eur_crt
dat <- dat[order(dat$sjlid, dat$strokedt, decreasing = FALSE),] # Sorted by earlier stroke events
#SJL1288701 had a row with strokegrade missing, but "stroke" variable showed it is a stroke
table(dat$stroke)
table(dat$strokegrade)
#dat$event <- ifelse(!is.na(dat$strokegrade), 1, 0)
dat$event <- ifelse(dat$stroke!="No stroke event", 1, 0)

### How many people had events?
length(unique(dat$sjlid[dat$event==1]))

# only use the variables needed so easier to check the data on screen
data <- dat[,c("sjlid","strokedt","event","dob","agelstcontact","agedx")]
data$first <- ave(data$agelstcontact, data$sjlid, FUN = seq_along)
M <- max(data$first)  ### maximum number of events
data[data$first==M,]$sjlid  #the id for the person with the maximum number of rows.
data[data$sjlid==data[data$first==M,]$sjlid,]

### Take the last row so we know the maximum number of events per person
event.number <- do.call(rbind, lapply(split(data, data$sjlid), tail, 1))[,c("sjlid","first")]
colnames(event.number) <- c("sjlid","maxE")

alldata <- merge(data,event.number,by.x="sjlid",by.y="sjlid")

alldata$start <- NULL
alldata$end <- NULL
### For those without event, start is agedx and end is Fu date === please modify based on your data, such as start from 5 years post DX
alldata$start[alldata$event==0] <- alldata$agedx[alldata$event==0]
alldata$end[alldata$event==0] <- alldata$agelstcontact[alldata$event==0]

### For the first event, start is agedx and end is first event time
alldata$start[alldata$event==1 & alldata$first==1] <- alldata$agedx[alldata$event==1 & alldata$first==1]
alldata$end[alldata$event==1 & alldata$first==1] <- as.numeric(difftime(alldata$strokedt[alldata$event==1 & alldata$first==1],alldata$dob[alldata$event==1 & alldata$first==1], units = 'days')/365.25)

#### For events that are not the first, segments are from the previous event date to this event date
alldata$previous <- as.Date(c(NA,alldata$strokedt[1:length(alldata$strokedt)-1]),origin = "1970-01-01")
alldata[1:10,]
### if first>1 (i.e, 2 or more events, previous event time remained, others are missing)
alldata$previous[alldata$first==1] <- NA
alldata$start[alldata$first>1] <- as.numeric(difftime(alldata$previous[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)
alldata$end[alldata$first>1] <- as.numeric(difftime(alldata$strokedt[alldata$first>1],alldata$dob[alldata$first>1], units = 'days')/365.25)


### If one person has only 1 event, we need to have segment from agedx to event date (handled above), and then from event date to end of FU
### If one person has multiple events, we need to add segments from the last event to end of Fu.
row_add <- alldata[alldata$event==1 & alldata$first==alldata$maxE,] ## this includes (1) people with only 1 event (2) The last row for people with multiple events
row_add$start <- as.numeric(difftime(row_add$strokedt, row_add$dob, units = 'days')/365.25)
row_add$end <- row_add$agelstcontact;
#### remember to make these rows event as 0, as we are adding the time after the last event date.
row_add$event <- 0
### If we do first event analysis, we would not need these. I need to add an indicator here, in case I only need to get the first event analysis segments.
row_add$evt1 <- 0;

alldata$evt1 <- 0;
alldata$evt1[alldata$first==1] <- 1
adata <- rbind(alldata,row_add)
#### order by person and start date
adata <- adata[order(adata$sjlid, adata$start, decreasing = FALSE),]
table(adata$event)## Double check event numebr is correct

###any stop time <=start time
any <- adata[adata$end<=adata$start,]
### 2 people had the stroke date on the Fu date. In this case, either get rid of the 2 lines [i.e, no time is follow-up after the last event date], or add 1 day on end date of these 2 segments, assuming there were followed up 1 more day. Will not make much difference. 1 day out of 365 days is 0.0027.
dim(adata)
final <- adata[adata$end>adata$start,]

minimum  <-  min(final$start, na.rm = TRUE) - 1
if (minimum < 0) minimum <- 0
maximum  <-  max(final$end, na.rm = TRUE) + 1
stroke_py <-  survSplit(final, cut=seq(minimum, maximum, 1), end="end",start="start",event="event") 
table(stroke_py$event)   ## Double check event numebr is correct
#### If you need the rows for first event analysis, take evt1=1
length(unique(stroke_py$sjlid[stroke_py$event==1]))
length(unique(stroke_py$sjlid[stroke_py$event==1 & stroke_py$evt1==1 ]))

stroke_py$PY <- stroke_py$end-stroke_py$start


geeglm(formula = event ~ agedx, family = "poisson", offset = log(stroke_py$PY), id = stroke_py$sjlid, corstr = "independence",  std.err = "san.se", data = stroke_py)
