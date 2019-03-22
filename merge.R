#! /usr/bin/R



setwd("/home/analysisTemp/TinyPlus/reference/pcoa/ZSL/control_tsv")
library(tidyverse)
temp = list.files(pattern="*_forward.tsv")

data = tibble(File = temp) %>%
  extract(File, "sample", remove = FALSE) %>%
  mutate(Data = lapply(File,read.csv,col_names=F)) %>%
  unnest(Data) %>%
  select(-File)

dat<-separate(data,X1,c("RNAME","TLEN","count"),sep = " ")
#names(dat)<-dat[1,]
miss<-dat[grep("|",dat$RNAME,fixed = T),]
miss<-miss[which(miss$sample!='ref'),]

temp=""
data=""

data1<-rbind(dat,miss)
data1$count<-as.numeric(data1$count)
dat1<-spread(data=unique(data1[,-3]),key = RNAME,value = count,fill=0)
dat2<-data.frame(unique(data1[,2:3]))
dat=""
dat3<-data.frame(t(dat1))
names(dat3)<-dat3[1,]
dat3$RNAME<-row.names(dat3)
dat4<-merge(dat2,dat3[-1,],by='RNAME')

#dat4[is.na(dat4)]<-0
write_tsv(dat4,"result/ref_count.tsv")