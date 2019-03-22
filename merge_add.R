#! /usr/bin/R

library(tidyverse)

ref_count<-read_table2("ref_count.tsv",col_names =T)
dat4<-ref_count
accession2tax<-read_table2("/home/kathleen/bacteria_gb.accession2taxid3",col_names=T)
A<-merge(accession2tax,dat4,by.y='RNAME',by.x='accession.version')
C<-read_table2("/home/analysisTemp/TinyPlus/tool/Reference_Taxonomy/names.tax",col_names=T)
D<-merge(C,A,by.y='taxid',by.x='id')
E<-(D[,-4])[,-1]

names<-read_table2("/home/analysisTemp/TinyPlus/reference/pcoa/control_chinese/names.tsv",col_names =F)
names(names)<-c("RNAME","Genus","Species")


ref_count_genus<-E
ref_count_genus$name<-sub("[","",ref_count_genus$name,fixed=T)
ref_count_genus$name<-sub("]","",ref_count_genus$name,fixed=T)
ref_count_genus$name<-sub("'","",ref_count_genus$name,fixed=T)
ref_count_genus$name<-sub("-like","",ref_count_genus$name,fixed=T)
ref_count_genus$name<-sub("eubacterium","Eubacterium",ref_count_genus$name,fixed=T)

dbc<-c(grep('uncultured',ref_count_genus$name),
       grep('Uncultured',ref_count_genus$name),
       grep('Candidatus',ref_count_genus$name),
       grep('unidentified',ref_count_genus$name),
       grep('UNVERIFIED',ref_count_genus$name),
       grep('candidate',ref_count_genus$name),
       grep('alph',ref_count_genus$name),
       grep('beta',ref_count_genus$name),
       grep('delta',ref_count_genus$name),
       grep('epsilon',ref_count_genus$name),
       grep('Clone',ref_count_genus$name),
       grep('bacterium',ref_count_genus$name),
       grep('human',ref_count_genus$name),
       grep('-',ref_count_genus$name,fixed=T),
       grep('.',ref_count_genus$name,fixed=T),
       grep(':',ref_count_genus$name,fixed=T),
       grep('Gram',ref_count_genus$name),
       row.names(dat4[which(is.element(dat4$RNAME,E$accession.version)==F),]))
dbc<-sort(unique(as.numeric(dbc)))

skip_all<-data.frame(ref_count_genus[(dbc),2])
names(skip_all)<-"RNAME"
skip_ed<-unique(merge(skip_all,(merge(merge(skip_all,names,by='RNAME'),dat4,by='RNAME'))))
skip_ed2<-skip_ed%>%group_by(RNAME)%>%summarise_all(funs(first))
D1<-skip_ed2[,-(3:4)][,-1]


D3<-ref_count_genus[-dbc,-(2:3)]
names(D3)<-names(D1)
D<-rbind(D1,D3)%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(D[,-1], 1, sum)
D1<-D[row.sums>5,]

write_tsv(D1,"ref_count_genus_filtered.tsv")