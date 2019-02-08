#!/bin/bash

# input data, edit this for input reads and contig files 
read_file=../data/hiv5_simulate_reads.fa
contig_file=../data/hiv5_contigs.fa

# contig index output path
if [ ! -d "./index" ]; then
    mkdir index
fi

index_out=./index

OUT=Reads_align_on_contig

L=0.6
bowtie2-build -f $contig_file $index_out"/Contig_index"
contig_ref=$index_out"/Contig_index"
bowtie2 -x $contig_ref -f --score-min G,20,8 --local -t -p 8 -S $OUT"_$L.sam" $read_file

samtools view -hS -F 4 $OUT"_$L.sam" >$OUT"_"$L"_mapped.sam"
samtools view -bS $OUT"_"$L"_mapped.sam" > $OUT"_"$L"_mapped.bam"
samtools sort $OUT"_"$L"_mapped.bam" $OUT"_"$L"_mapped_sort"
samtools mpileup -f $contig_file $OUT"_"$L"_mapped_sort.bam" -D -u >$OUT"_"$L".bcf"

bcftools view $OUT"_"$L".bcf" >$OUT"_"$L".vcf"


