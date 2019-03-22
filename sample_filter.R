library(vegan)
library(tidyverse)


# rarecurve would indiate whether the sequecing depth is 
# enough to capture majority of the species in the environment
# "good" sample would plateau-off after XXXX base pair


# input 
# ref_count.tsv               from TinyRefDis.py
# names.txt                   from gitlab
# output
# Diabetic_filter.svg         filtering samples by eye
# Diabetic_filter.tsv         for keeping input for rarecurve



setwd("D:/mothur/curve_filter")
ref_count <- read_delim("D:/mothur/Diabetic/ref_count.tsv","\t", escape_double = FALSE, trim_ws = TRUE)
ref_count[,-(1:2)]<-as.numeric(as.matrix(ref_count[,-(1:2)]))
names <- read_table2("D:/mothur/R/names.txt", col_names = FALSE)[,1:3]
names(names)<-c("RNAME","Genus","Species")

# match the NC id with species namee
# names and ref_count has a common name of RNAME
A<-merge(names,ref_count)
# add all count from the same genus together
# the B is a table with genus name in the first row 
# with one sample per column
# -(3:4) removes species name and gene length
# -1 removes the NC is of the sequences
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))

# the B[,-1] contains all count data
# B[,1] is the genus name
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)

# eliminate genus that have total count lower than 2 and sample with less than 50 reads
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
Diabetic_raw<-data.frame(t(D[,-1]))
names(Diabetic_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(Diabetic_raw[,-1], 1, sum)
col.sums <- apply(Diabetic_raw[,-1], 2, sum)
C<-cbind(row.sums,Diabetic_raw)
Diabetic_filter<-C[which(C$row.sums>50),-1]
# Diabetic_filter contains count data only
# each row is a sample 
# each column is a genus


# clean up the made-up table for filtering and formatting
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""

sample_size=min(rowSums(Diabetic_filter))
sample_pool<-rarefy(Diabetic_filter,sample_size)

svg("Diabetic_filter.svg",width = 20, height=30)
rarecurve(Diabetic_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()



# for Parkinson case
ref_count <- read.delim("D:/mothur/Parkinson/PRJEB14928.count")
ref_count[,-(1:2)]<-as.numeric(as.matrix(ref_count[,-(1:2)]))
A<-merge(names,ref_count)
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
Parkinson_raw<-data.frame(t(D[,-1]))
names(Parkinson_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(Parkinson_raw[,-1], 1, sum)
col.sums <- apply(Parkinson_raw[,-1], 2, sum)
C<-cbind(row.sums,Parkinson_raw)
Parkinson_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(Diabetic_filter))
sample_pool<-rarefy(Diabetic_filter,sample_size)
svg("Parkinson_filter.svg",width = 20, height=30)
rarecurve(Parkinson_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()


# for Parkinson case
ref_count <- read.delim("D:/mothur/Parkinson/PRJEB14928.count")
ref_count[,-(1:2)]<-as.numeric(as.matrix(ref_count[,-(1:2)]))
A<-merge(names,ref_count)
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
Parkinson_raw<-data.frame(t(D[,-1]))
names(Parkinson_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(Parkinson_raw[,-1], 1, sum)
col.sums <- apply(Parkinson_raw[,-1], 2, sum)
C<-cbind(row.sums,Parkinson_raw)
Parkinson_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(Diabetic_filter))
sample_pool<-rarefy(Diabetic_filter,sample_size)
svg("Parkinson_filter.svg",width = 20, height=30)
rarecurve(Parkinson_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()


# for gallstone
ref_count <- read_delim("D:/mothur/gallstone/v1v2/ref_count.bk.tsv","\t", escape_double = FALSE, trim_ws = TRUE)
ref_count[,-(1:2)]<-as.numeric(as.matrix(ref_count[,-(1:2)]))
A<-merge(names,ref_count)
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
gallstone_raw<-data.frame(t(D[,-1]))
names(gallstone_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(gallstone_raw[,-1], 1, sum)
col.sums <- apply(gallstone_raw[,-1], 2, sum)
C<-cbind(row.sums,gallstone_raw)
gallstone_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(gallstone_filter))
sample_pool<-rarefy(gallstone_filter,sample_size)
svg("gallstone_filter.svg",width = 20, height=30)
rarecurve(gallstone_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()

# for autism
ref_count <- read_delim("D:/mothur/autistism/kid_study/PRJEB15418.count","\t", escape_double = FALSE, trim_ws = TRUE)
ref_count[,-(1:2)]<-as.numeric(as.matrix(ref_count[,-(1:2)]))
A<-merge(names,ref_count)
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
autism_raw<-data.frame(t(D[,-1]))
names(autism_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(autism_raw[,-1], 1, sum)
col.sums <- apply(autism_raw[,-1], 2, sum)
C<-cbind(row.sums,autism_raw)
autism_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(autism_filter))
sample_pool<-rarefy(autism_filter,sample_size)
svg("autism_filter.svg",width = 20, height=30)
rarecurve(autism_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()



# for obese
ref_count <- read_delim("D:/mothur/obese/ref_count.tsv","\t", escape_double = FALSE, trim_ws = TRUE)
ref_count[,-(1:2)]<-as.numeric(as.matrix(ref_count[,-(1:2)]))
A<-merge(names,ref_count)
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
obese_raw<-data.frame(t(D[,-1]))
names(obese_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(obese_raw[,-1], 1, sum)
col.sums <- apply(obese_raw[,-1], 2, sum)
C<-cbind(row.sums,obese_raw)
obese_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(obese_filter))
sample_pool<-rarefy(obese_filter,sample_size)
svg("obese_filter.svg",width = 20, height=30)
rarecurve(obese_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()

# for IBS case data
ref_count <- read_delim("D:/mothur/IBS/ref_count.tsv","\t", escape_double = FALSE, trim_ws = TRUE)

A<-merge(names,ref_count)
miss<-merge(A[grep(".",A$Genus,fixed=T),1:2],ref_count_genus[,],by.x='RNAME',by.y='accession.version')

A1<-A[which(is.element(A$RNAME,miss$RNAME)==F),-3]
A2<-miss[,-2]
names(A2)<-names(A1)
A<-rbind(A2,A1)

B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
IBS_raw<-data.frame(t(D[,-1]))
names(IBS_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(IBS_raw[,-1], 1, sum)
col.sums <- apply(IBS_raw[,-1], 2, sum)
C<-cbind(row.sums,IBS_raw)
IBS_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(IBS_filter))
sample_pool<-rarefy(IBS_filter,sample_size)
svg("IBS_filter.svg",width = 20, height=30)
rarecurve(IBS_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()


# for MS
PRJNA321051_count <- read_delim("D:/mothur/MS/PRJNA321051.count.s","\t", escape_double = FALSE, trim_ws = TRUE)
PRJNA321051 <- read_delim("D:/mothur/MS/PRJNA321051.count","\t", escape_double = FALSE, trim_ws = TRUE)
ref_count<-merge(PRJNA321051,PRJNA321051_count[,-2],all=T,by="RNAME")
A<-merge(names,ref_count)
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
MS_raw<-data.frame(t(D[,-1]))
names(MS_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(MS_raw[,-1], 1, sum)
col.sums <- apply(MS_raw[,-1], 2, sum)
C<-cbind(row.sums,MS_raw)
MS_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(MS_filter))
sample_pool<-rarefy(MS_filter,sample_size)
svg("MS_filter.svg",width = 20, height=30)
rarecurve(MS_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()


# for CRC
A<-merge(names,ref_count)
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
CRC_raw<-data.frame(t(D[,-1]))
names(CRC_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(CRC_raw[,-1], 1, sum)
col.sums <- apply(CRC_raw[,-1], 2, sum)
C<-cbind(row.sums,CRC_raw)
CRC_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(CRC_filter))
sample_pool<-rarefy(CRC_filter,sample_size)
svg("CRC_filter.svg",width = 20, height=30)
rarecurve(CRC_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()

# for ZSL control
ref_count <- read_table2("D:/mothur/ZSL/control_result/test.tsv")
names <- read_table2("C:/Users/kathleen/Downloads/16S_addin/new/head/processed_names.txt",col_names=F)[,1:3]

#SMicrobial <- read_table2("C:/Users/kathleen/Downloads/pathogen_addin/Addin_1/SMicrobial.head",col_names = FALSE)[,1:5]
#names_add<-SMicrobial[which(is.element(SMicrobial$X1,names$RNAME)==F),]
#write_delim(names_add,"C:/Users/kathleen/Downloads/pathogen_addin/Addin_1/add_name.head"," ")
names_add<- read_table2("C:/Users/kathleen/Downloads/pathogen_addin/Addin_1/add_name.head",col_names=F)[,1:3]
name2<-rbind(names_add,names)
names(name2)<-c("RNAME","Genus","Species")


ref_count[,-(1:2)]<-as.numeric(as.matrix(ref_count[,-(1:2)]))
A<-merge(name2,ref_count_IBS)
B<-A[,-(3:4)][,-1]%>%group_by(Genus)%>%summarise_all(funs(sum))
row.sums <- apply(B[,-1], 1, sum)
col.sums <- apply(B[,-1], 2, sum)
C<-cbind(row.sums,B)
D<-C[which(C$row.sums>2),-1]
ZSL_control_raw<-data.frame(t(D[,-1]))
names(ZSL_control_raw)<-as.character(unlist(D[,1]))
row.sums <- apply(ZSL_control_raw[,-1], 1, sum)
col.sums <- apply(ZSL_control_raw[,-1], 2, sum)
C<-cbind(row.sums,ZSL_control_raw)
ZSL_control_filter<-C[which(C$row.sums>50),-1]
A<-""
B<-""
C<-""
D<-""
row.sums<-""
col.sums <-""
sample_size=min(rowSums(ZSL_control_filter))
svg("D:/mothur/curve_filter/ZSL_control_filter.svg",width = 20, height=30)
rarecurve(ZSL_control_filter,step=20,sample=sample_size,col="blue",cex=0.6)
dev.off()



# remove specified column/sample
# not needed for parkinson and obese data

# MS
rm<-c("SRR3501989",
  "SRR3502010",
  "SRR3581134",
  "SRR3581133",
  "SRR3581132",
  "SRR3501915",
  "SRR3501912",
  "SRR3501936",
  "SRR3501925",
  "SRR3501997")
genus_filtered<-MS_filter[!rownames(MS_filter) %in% rm, ]
write.csv(genus_filtered,file ="D:/mothur/MS/genus_count.csv")
case_selected<-ref_count[,!names(ref_count) %in% rm]
write.csv(case_selected,file ="D:/mothur/MS/ref_count_filter.csv")



# diabetics
Diabetic_filter
rm<-c("SRR1276219",
      "SRR1276229",
      "SRR1276234",
      "SRR1276245",
      "SRR1276246",
      "SRR1276252",
      "SRR1276255")
case_selected<-ref_count[,!names(ref_count) %in% rm]
write_tsv(case_selected,path = "ref_count_filter.tsv",col_names = T)

# gallstone
rm<-c("SRR946345",
      "SRR946355",
      "SRR946354",
      "SRR946269",
      "SRR946326",
      "SRR946124",
      "SRR1276255")
ref_count <- read_delim("D:/mothur/gallstone/v1v2/ref_count.bk.tsv","\t", escape_double = FALSE, trim_ws = TRUE)
case_selected<-ref_count[,!names(ref_count) %in% rm]
write_tsv(case_selected,path = "ref_count_filter.tsv",col_names = T)

# autism
rm<-c("ERR1637179",
      "ERR1637196",
      "ERR1637204",
      "ERR1637213",
      "ERR1637216",
      "ERR1637185")
ref_count <- read_delim("D:/mothur/autistism/kid_study/PRJEB15418.count","\t", escape_double = FALSE, trim_ws = TRUE)
case_selected<-ref_count[,!names(ref_count) %in% rm]
write_tsv(case_selected,path = "ref_count_filter.tsv",col_names = T)

# IBS
rm<-c("SRR5577168")
ref_count <- read_delim("D:/mothur/IBS/PRJNA386442.count", "\t", escape_double = FALSE, trim_ws = TRUE)
case_selected<-ref_count[,!names(ref_count) %in% rm]
write_tsv(case_selected,path = "ref_count_filter.tsv",col_names = T)

# zsl control

rm<-c("ERR2017436",
      "ERR2017444",
      "ERR2017588",
      "ERR2017589",
      "ERR2017590",
      "ERR2017591",
      "ERR2017592",
      "ERR2017593",
      "ERR2017594",
      "ERR2017595",
      "ERR2017596",
      "ERR2017597",
      "ERR2017598",
      "ERR2017599",
      "ERR2017523",
      "ERR2017443",
      "ERR2017437",
      "ERR2017559",
      "ERR2017453",
      "ERR2017450",
      "ERR2017427",
      "ERR2017449",
      "ERR2017421",
      "ERR2017459",
      "ERR2017455",
      "ERR2017413",
      "ERR2017456",
      "ERR2017571",
      "ERR2017411",
      "ERR2017423",
      "ERR2017420",
      "ERR2017474",
      "ERR2017430",
      "ERR2017566",
      "ERR2017440",
      "ERR2017432",
      "ERR2017473",
      "ERR2017416",
      "ERR2017469",
      "ERR2017580",
      "ERR2017520",
      "ERR2017435",
      "ERR2017464",
      "ERR2017563",
      "ERR2017506",
      "ERR2017460",
      "ERR2017420",
      "ERR2017415",
      "ERR2017433",
      "ERR2017414",
      "ERR2017482",
      "ERR2017564",
      "ERR2017495",
      "ERR2017422",
      "ERR2017533",
      "ERR2017452",
      "ERR2017451",
      "ERR2017454",
      "ERR2017461",
      "ERR2017483",
      "ERR2017467",
      "ERR2017554",
      "ERR2017558",
      "ERR2017568",
      "ERR2017563")
ref_count <- read_table2("D:/mothur/ZSL/control_result/test.tsv")
case_selected<-ref_count[,!names(ref_count) %in% rm]
write_tsv(case_selected,path = "D:/mothur/ZSL/control_result/ref_count_filter.tsv",col_names = T)
