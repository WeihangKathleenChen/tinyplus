#! /usr/bin/R
library(tidyverse)
library(readr)
# setwd("C:/Users/kathleen/Desktop/Scoring")

# man sed | sed -i "s/ //g" ref_count.tsv
# man sed | sed -i "s/ //g" ref_filtered_count.txt
count<-read_delim("ref_count.tsv","\t") # total count table
ref_count<-read_delim("ref_count.tsv","\t",
                               col_types=rep(c(col_character(),col_double()),times=c(1,length(names(count[,-(1)])))))
count<-""
filtered_count<-read_delim("ref_filtered_count.txt","\t") # total readable table
col_type<-rep(c(col_character(),col_integer()),times=c(1,length(names(filtered_count[,-(1)]))))
ref_filtered_count<-read_delim("ref_filtered_count.txt","\t",col_names = T,
                               col_types=col_type)
filtered_count<-""
ref_filtered_count_processed<-read_delim("ref_count.out","\t") # lineage output
names(ref_filtered_count_processed)<-c("RNAME",names(ref_filtered_count_processed[,-1]))


####since no accession failed in finding its genus,
####filtered_count table is produced by binding instead of merging
filtered_count<-unique(merge(ref_filtered_count_processed,ref_filtered_count))


######## to sum counts by genus, all taxon info were removed except genus #####
genus_count<-as_tibble(filtered_count[,-(1:7)][,-(2:3)])%>%
  group_by(genus)%>%
  summarise_all(funs(sum))
# genus_count <- read_delim("C:/Users/kathleen/Desktop/genus_count.tsv",\t", escape_double = FALSE, trim_ws = TRUE)

family_count<-na.omit(as_tibble(subset(filtered_count,is.na(filtered_count$genus))[,-(1:6)][,-(2:4)])) %>%
  group_by(family)%>%
  summarise_all(funs(sum))

class_count <-na.omit(subset(subset(filtered_count,is.na(filtered_count$genus)),is.na(filtered_count$family))[,-(1:4)][-(2:6)]) %>%
  group_by(class)%>%
  summarise_all(funs(sum))

genus_list<-genus_count[,1] # list of all genus in alphabetic order

names(genus_count)<-c("name",names(genus_count[,-1]))
names(family_count)<-c("name",names(genus_count[,-1]))
names(class_count)<-c("name",names(genus_count[,-1]))
tot_count<-rbind(na.omit(genus_count),family_count,class_count)
write_tsv(tot_count,"tot_count.tsv")

########## sumy : summary table of count infomation per sample based on assigned genus
sumy<-as_tibble(names(genus_count[,-1])) # list all sample id as first column
names(sumy)<-"sample"
tot_counts<-as_tibble(colSums(ref_count[,-(1:2)]))
names(tot_counts)<-"tot_counts"
tot_counts$sample<-row.names(tot_counts)
summy<-merge(sumy,tot_counts)
summy$filter_count<-colSums(genus_count[,-(1)])
summy$assigned_count<-colSums(genus_count[,-1])
summy$Max_genus_count<-as.numeric(apply(as.matrix(na.omit(genus_count[,-1])),2,max))
summy<-na.omit(summy)
########### uses the count number to identidy the genus with maximum count after filteration by mismatch
y<-genus_list[(genus_count[,names(genus_count) %in% summy[1,1]])==as.numeric(summy[1,5]),1]

y<-vector("list",length(summy$sample))
for (i in 1:length(summy$sample)){
  y[i]<-genus_list[(genus_count[,names(genus_count) %in% summy[i,1]])==as.numeric(summy[i,5]),1]
}
summy$Max_genus<-as.character(y)
summy$MainGenus<-(summy$Max_genus_count/summy$tot_counts)
summy$annotated<-(summy$assigned_count/summy$tot_counts)
summy$genus_annoted<-colSums(tot_count[,-1]>0)
write_tsv(summy,"summy.tsv")

sumy<- summy[(summy$annotated > 0.25 & summy$MainGenus > 0.02),]
sum_Bacteroides<-sumy[which(sumy$Max_genus=='Bacteroides'),]
sum_Prevotella<-sumy[which(sumy$Max_genus=='Prevotella'),]
sum_Other<-sumy[which(sumy$Max_genus!='Prevotella'),][which(sumy$Max_genus!='Bacteroides'),]

############### to be continued #############
case_Bacteroides<-cbind(tot_count[,1],tot_count[,names(tot_count)%in%sum_Bacteroides$sample])
ifelse(length(names(case_Bacteroides)) > 1, write_tsv(case_Bacteroides,"case_Bacteroides.tsv"), print("no Bacteroides found"))
case_Prevotella<-cbind(tot_count[,1],tot_count[,names(tot_count)%in%sum_Prevotella$sample])
ifelse(length(names(case_Prevotella)) > 1, write_tsv(case_Prevotella,"case_Prevotellas.tsv"), print("no Prevotella found"))
case_Other<-cbind(tot_count[,1],tot_count[,names(tot_count)%in%sum_Other$sample])
ifelse(length(names(case_Other)) > 1, write_tsv(case_Other,"case_Other.tsv"), print("no Other found"))