library(tidyr)
library('ggplot2')
library(dplyr)

tbl=read.table("C:\\Users\\BLenny\\Documents\\final.5.Q_header2_SJLIFE_only",header = TRUE)
tbl_long<-tbl %>% gather(Ancestry,Value,2:6)
p<-ggplot(tbl_long,aes(x=INDIVIDUAL,y=Value,fill=Ancestry))+geom_col(position=position_stack())
show(p)

