#!bin/bash

bowtie_index=$1
data_path=$2
output=$2/result
cd $data_path
mkdir -p result
mkdir -p raw_data
echo $data_path

for i in $(ls *_1.fastq.gz | rev | cut -c 12- | rev | uniq)

do
# keep only the mapped reads
bowtie2 -p 8 --rg-id ${i} --very-sensitive-local -x $bowtie_index -1 ${i}_1.fastq.gz -2 ${i}_2.fastq.gz | samtools view -Shub -F 4 | samtools sort -o ${i}.sorted.bam 

# count only the forward strand as the effective read for bacteria
samtools view -b -F 0x40 ${i}.sorted.bam > ${i}.forward.bam
samtools index ${i}.forward.bam
samtools idxstats ${i}.forward.bam | tr ' ' '\t' | awk '{if ($3!=0) {print $1, $2, $3} }' > ${i}_forward.tsv

# extract fastq file back from bam
# samtools sort -n -o ${i}.sorted_by_name.bam ${i}.sorted.bam
# samtools fixmate ${i}.sorted_by_name.bam ${i}.fixed.bam
# bamToFastq -i ${i}.fixed.bam -fq ${i}_16S.fq
# rm ${i}.sorted_by_name.bam ${i}.fixed.bam

# filter out reads with high mismatch
samtools view -h ${i}.sorted.bam | grep -v 'NM:i:[4-9]' | grep -v 'NM:i:[1-3][0-9]' | samtools view -b -F 0x40 > ${i}_filtered.forward.bam
samtools index ${i}_filtered.forward.bam
samtools idxstats ${i}_filtered.forward.bam | tr ' ' '\t' | awk '{if ($3!=0) {print $1, $2, $3} }' > ${i}_filtered_forward.txt

# move the prcoessed files to a subdirectory named after sample ID under the working directory

# mv ${i}_1.fastq.gz raw_data
# mv ${i}_2.fastq.gz raw_data


# copy the forward count file to another subdirectory called result 
# an R script called merge.R will merge all reads from different samples in /result 
cp $data_path/${i}_forward.tsv $output
cp $data_path/${i}_filtered_forward.txt $output

done

# for single end fastq data
for i in $(ls $data_path | grep [0-9].fastq.gz | rev | cut -c 10- | rev | uniq)

do

bowtie2 -p 8 --rg-id ${i} --very-sensitive-local -x $bowtie_index -U ${i}.fastq.gz | samtools view -Shub -F 4 | samtools sort -o ${i}.sorted.bam 

# we only count the forward strand as the effective read for bacteria for now.
samtools view -b -F 0x40 ${i}.sorted.bam > ${i}.forward.bam
samtools index ${i}.forward.bam
samtools idxstats ${i}.forward.bam | tr ' ' '\t' | awk '{if ($3!=0) {print $1, $2, $3} }' > ${i}_forward.tsv

# extract fastq file back from bam
# samtools sort -n -o ${i}.sorted_by_name.bam ${i}.sorted.bam
# samtools fixmate ${i}.sorted_by_name.bam ${i}.fixed.bam
# bamToFastq -i ${i}.fixed.bam -fq ${i}_16S.fq
# rm ${i}.sorted_by_name.bam ${i}.fixed.bam

# filter by mismatch
samtools view -h ${i}.sorted.bam | grep -v 'NM:i:[4-9]' | grep -v 'NM:i:[1-3][0-9]' | samtools view -b -F 0x40 > ${i}_filtered.forward.bam
samtools index ${i}_filtered.forward.bam
samtools idxstats ${i}_filtered.forward.bam | tr ' ' '\t' | awk '{if ($3!=0) {print $1, $2, $3} }' > ${i}_filtered_forward.txt

# move the prcoessed files to a subdirectory named after sample ID under the working directory
# mv ${i}.fastq.gz raw_data

# copy the forward count file to another subdirectory called result 
# an R script called merge.R will merge all reads from different samples in /result 
cp $data_path/${i}_forward.tsv $output
cp $data_path/${i}_filtered_forward.txt $output

done

cd $output
Rscript /home/analysisTemp/TinyPlus/reference/pcoa/control_chinese/merge2.R &
Rscript /home/analysisTemp/TinyPlus/reference/pcoa/control_chinese/merge.R
python36 /home/analysisTemp/TinyPlus/tool/Reference_Taxonomy/tinyTax.py lineage -i ref_count.tsv -o ref_count.out -f ref_count.fail
man sed | sed -i "s/ //g" ref_count.tsv
man sed | sed -i "s/ //g" ref_filtered_count.txt
# Rscript /home/analysisTemp/TinyPlus/reference/pcoa/control_chinese/count_preprocess.R