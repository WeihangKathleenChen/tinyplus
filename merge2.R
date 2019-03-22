#! /usr/bin/R
################################# merge all sample count to 1 table##########
library(tidyverse)
library(readr)
temp = list.files(pattern="*_filtered_forward.txt")

data = tibble(File = temp) %>%
  extract(File, "sample", remove = FALSE) %>%
  mutate(Data = lapply(File,read_csv,col_names=F)) %>%
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
sample_id<-data.frame(dat1[,1])
dat2<-as_data_frame(unique(data1[,2:3]))
dat=""
dat3<-data.frame(t(dat1))
names(dat3)<-dat3[1,]
dat3$RNAME<-row.names(dat3)
dat4<-merge(dat2,dat3[-1,],by='RNAME')
names(dat4)<-c(names(dat4[,(1:2)]),as.character(dat1$sample))
write_tsv(dat4,"ref_filtered_count.txt")

# ################ load table ##########
# accession2tax<-read_table2("/home/kathleen/bacteria_gb.accession2taxid3",col_names=T)
# names<-read_table2("/home/analysisTemp/TinyPlus/reference/pcoa/control_chinese/names.tsv",col_names =T)
# ################ get genus from RNAME ###########################
# A<-merge(accession2tax,dat4,by.y='RNAME',by.x='accession.version')
# C<-read_table2("/home/analysisTemp/TinyPlus/tool/Reference_Taxonomy/names.tax",col_names=T)
# D<-merge(C,A,by.y='taxid',by.x='id')
# E<-(D[,-4])[,-1]

# ################ get genus with unsepecific name ###############

# ref_count_genus<-E
# ref_count_genus$name<-sub("[","",ref_count_genus$name,fixed=T)
# ref_count_genus$name<-sub("]","",ref_count_genus$name,fixed=T)
# ref_count_genus$name<-sub("'","",ref_count_genus$name,fixed=T)
# ref_count_genus$name<-sub("-like","",ref_count_genus$name,fixed=T)
# ref_count_genus$name<-sub("eubacterium","Eubacterium",ref_count_genus$name,fixed=T)

# dbc<-c(grep('uncultured',ref_count_genus$name),
#        grep('Uncultured',ref_count_genus$name),
#        grep('Candidatus',ref_count_genus$name),
#        grep('unidentified',ref_count_genus$name),
#        grep('UNVERIFIED',ref_count_genus$name),
#        grep('candidate',ref_count_genus$name),
#        grep('alph',ref_count_genus$name),
#        grep('beta',ref_count_genus$name),
#        grep('delta',ref_count_genus$name),
#        grep('epsilon',ref_count_genus$name),
#        grep('Clone',ref_count_genus$name),
#        grep('bacterium',ref_count_genus$name),
#        grep('human',ref_count_genus$name),
#        grep('-',ref_count_genus$name,fixed=T),
#        grep('.',ref_count_genus$name,fixed=T),
#        grep(':',ref_count_genus$name,fixed=T),
#        grep('Gram',ref_count_genus$name),
#        row.names(dat4[which(is.element(dat4$RNAME,E$accession.version)==F),]))
# dbc<-sort(unique(as.numeric(dbc)))


# skip_all<-data.frame(ref_count_genus[(dbc),2])
# names(skip_all)<-"RNAME"
# skip_ed<-merge(skip_all,(merge(merge(skip_all,names,by='RNAME'),dat4,by='RNAME')))
# write_tsv(skip_ed,"skip_ed.tsv")
# skip_ed2<-skip_ed%>%group_by(RNAME)%>%summarise_all(funs(first))
# D1<-skip_ed2[,-(3:4)][,-1]


# D3<-ref_count_genus[-dbc,-(2:3)]
# names(D3)<-names(D1)
# D<-rbind(D1,D3)%>%group_by(Genus)%>%summarise_all(funs(sum))

# ################## filter low count genus ####################################
# row.sums <- apply(D[,-1], 1, sum)
# D2<-D[row.sums>5,]

# write_tsv(D1,"D1.tsv")
# write_tsv(D3,"D3.tsv")
# write_tsv(D2,"ref_count_genus_filtered.tsv")