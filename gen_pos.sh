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

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf

chmod 777 ./genePredToGtf

cut -f2-11 hg19.txt | ./genePredToGtf file stdin hg19.gtf

# sort to match fasta file chr10, chr11, chr12 etc etc..
sort -k1 -n hg19.gtf > tmp.gtf && rm hg19.gtf && mv tmp.gtf hg19.gtf
sort -k3 -n hg19.txt > tmp.txt && rm hg19.txt && mv tmp.txt hg19.txt


mkdir -p ucsc 

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz && gunzip hg19.fa.gz

awk '/^>/ {F=substr($0, 2, length($0))".fa"; print >F;next;} {print >> F;}' < hg19.fa && mv *.fa ucsc/

sleep 5

rm ucsc/hg19.fa

while read -r line; do rm -f ucsc/${line}.fa; done < non_canon_chrs.txt

cat ucsc/*.fa > ucsc/hg19.fa

## 4. Generate Sim reads 

mkdir -p reads
mkdir -p logs

perl ./CIRIsimulator.pl -1 reads/glio_1.fastq -2 reads/glio_2.fastq -O logs/glio_circRNA.list -G ./hg19.gtf -DB truth_sets/glioblastoma.txt -C 10 -LC 0 -R 1 -LR 1 -L 101 -E 1 -I 350 -D ./ucsc -CHR1 0 -M 50 2>&1 | tee logs/glio.out && \
gzip reads/glio_1.fastq reads/glio_2.fastq

perl ./CIRIsimulator.pl -1 reads/oligo_1.fastq -2 reads/oligo_2.fastq -O logs/oligo_circRNA.list -G ./hg19.gtf -DB truth_sets/oligodendroma.txt -C 10 -LC 0 -R 1 -LR 1 -L 101 -E 1 -I 350 -D ./ucsc -CHR1 0 -M 50 2>&1 | tee logs/oligo.out && \
gzip reads/oligo_1.fastq reads/oligo_2.fastq

perl ./CIRIsimulator.pl -1 reads/leuk_1.fastq -2 reads/leuk_2.fastq -O logs/leuk_circRNA.list -G ./hg19.gtf -DB truth_sets/leukemia.txt -C 10 -LC 0 -R 1 -LR 1 -L 101 -E 1 -I 350 -D ./ucsc -CHR1 0 -M 50 2>&1 | tee logs/leuk.out && \
gzip reads/leuk_1.fastq reads/leuk_2.fastq

## 5. Run Workflow

nextflow -bg run -r dev nf-core/circrna -profile singularity,lugh --input "reads/*_{1,2}.fastq.gz" --input_type "fastq" --module "circrna_discovery" --tool "circexplorer2, ciriquant, circrna_finder, dcc, find_circ, mapsplice, segemehl" --bsj_reads 0 --tool_filter 0 --outdir "simulated_analysis" --fasta "hg19.fa" --gtf "hg19.gtf" --circexplorer2_annotation "hg19.txt" --max_cpus 16 --max_memory '100.GB'

## note: you will have to extract each cancer type from the results and check against the truth sets (i.e dont count a circRNA exclusive to glioblastoma as a missed call for the other 2). 
## Generate coordinates from the truth sets and use these as members of a set for analysis in R. 
