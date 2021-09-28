#!/usr/bin/env bash 

## Generate simulated circRNA reads based on Glioblastoma circRNAs from circRNADb
## Barry Digby 

## 1. Download circRNADb database file.

wget http://reprod.njmu.edu.cn/circrnadb/doc/circRNA_dataset.zip && unzip circRNA_dataset.zip

## 2. Filter for glioblastoma, and reformat for CIRIsimulator.pl

grep "glioblastoma" circRNA_dataset/circRNA_dataset.txt | awk -v OFS="\t" '{print $2, $3, $4, 0, 0, 0, $9, 0, 0, 0, $8}' | sort > glioblastoma_reformat4CIRIsimulator.txt

## 3. Prep ucsc files

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz && gunzip refGene.txt.gz && mv refGene.txt hg19.txt

wget wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf && chmod +x genePredToGtf

cut -f2-11 hg19.txt | ./genePredToGtf file stdin hg19.gtf

mkdir -p ucsc 

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz && gunzip hg19.fa.gz

awk '/^>/ {F=substr($0, 2, length($0))".fa"; print >F;next;} {print >> F;}' < hg19.fa && mv *.fa ucsc/

## 4. Generate Sim reads 

perl ./CIRIsimulator.pl -1 pos_1.fastq -2 pos_2.fastq -O pos_circRNA.list -G ./hg19.gtf -DB ./glioblastoma_reformat4CIRIsimulator.txt -C 10 -LC 0 -R 1 -LR 1 -L 101 -E 1 -I 350 -D ./ucsc -CHR1 0 -M 50 2>&1 | tee gen_pos.out && \
gzip pos_1.fastq pos_2.fastq
