#!/usr/bin/env bash 

## Generate simulated circRNA reads based on Glioblastoma, Oligodendroma, Leukemia circRNAs from circRNADb
## Barry Digby 

## 1. Download circRNADb database file.

wget http://reprod.njmu.edu.cn/circrnadb/doc/circRNA_dataset.zip && unzip circRNA_dataset.zip

## 2. Filter for cancers, reformat for CIRIsimulator.pl

mkdir -p truth_sets

grep "glioblastoma" circRNA_dataset/circRNA_dataset.txt | awk -v OFS="\t" '{print $2, $3, $4, 0, 0, 0, $9, 0, 0, 0, $8}' | sort > truth_sets/glioblastoma.txt

grep "oligodendroma" circRNA_dataset/circRNA_dataset.txt | awk -v OFS="\t" '{print $2, $3, $4, 0, 0, 0, $9, 0, 0, 0, $8}' | sort > truth_sets/oligodendroma.txt

grep "leukemia" circRNA_dataset/circRNA_dataset.txt | awk -v OFS="\t" '{print $2, $3, $4, 0, 0, 0, $9, 0, 0, 0, $8}' | sort > truth_sets/leukemia.txt

## 3. Prep ucsc files

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz && gunzip refGene.txt.gz && mv refGene.txt hg19.txt

grep -vf non_canon_chrs.txt hg19.txt > tmp.txt && rm hg19.txt && mv tmp.txt hg19.txt

wget wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf && chmod +x genePredToGtf

cut -f2-11 hg19.txt | ./genePredToGtf file stdin hg19.gtf

mkdir -p ucsc 

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz && gunzip hg19.fa.gz

awk '/^>/ {F=substr($0, 2, length($0))".fa"; print >F;next;} {print >> F;}' < hg19.fa && mv *.fa ucsc/

while read -r line; do rm ucsc/${line}.fa; done < non_canon_chrs.txt

## 4. Generate Sim reads 

mkdir -p reads
mkdir -p logs

perl ./CIRIsimulator.pl -1 reads/glio_1.fastq -2 reads/glio_2.fastq -O logs/glio_circRNA.list -G ./hg19.gtf -DB truth_sets/glioblastoma.txt -C 10 -LC 0 -R 1 -LR 1 -L 101 -E 1 -I 350 -D ./ucsc -CHR1 0 -M 50 2>&1 | tee logs/glio.out && \
gzip reads/glio_1.fastq reads/glio_2.fastq

perl ./CIRIsimulator.pl -1 reads/oligo_1.fastq -2 reads/oligo_2.fastq -O logs/oligo_circRNA.list -G ./hg19.gtf -DB truth_sets/oligodendroma.txt -C 10 -LC 0 -R 1 -LR 1 -L 101 -E 1 -I 350 -D ./ucsc -CHR1 0 -M 50 2>&1 | tee logs/oligo.out && \
gzip reads/oligo_1.fastq reads/oligo_2.fastq

perl ./CIRIsimulator.pl -1 reads/leuk_1.fastq -2 reads/leuk_2.fastq -O logs/leuk_circRNA.list -G ./hg19.gtf -DB truth_sets/leukemia.txt -C 10 -LC 0 -R 1 -LR 1 -L 101 -E 1 -I 350 -D ./ucsc -CHR1 0 -M 50 2>&1 | tee logs/leuk.out && \
gzip reads/leuk_1.fastq reads/leuk_2.fastq
