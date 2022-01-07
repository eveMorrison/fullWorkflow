import os, sys

input_dir = '/Users/evamorrison/Desktop/Eva/fullWorkflow'

os.chdir(input_dir)

# filter repeatmasker output

awk_filter_nano = "awk '{if (($7-$6) >= 250) print $0;}' nanoporeParis.fa.out | awk '{if($2<=20.0) print $0;}' > nanoporeParis_filter.fa.out"
awk_filter_160 = "awk '{if (($7-$6) >= 250) print $0;}' strain160all.fasta.out | awk '{if($2<=20.0) print $0;}' > strain160all_filter.fa.out"
awk_filter_9 = "awk '{if (($7-$6) >= 250) print $0;}' strain9all.fasta.out | awk '{if($2<=20.0) print $0;}' > strain9all_filter.fa.out"

os.system(awk_filter_nano)
os.system(awk_filter_160)
os.system(awk_filter_9)

# Find the reads containing the TE

awk_find_nano = '''awk '{if (($10) == "Paris") print $0;}' nanoporeParis_filter.fa.out | awk '{Nanopore[$5];} END{for (var in Nanopore) print var;} ' > nanopore_reads_Paris.txt'''

os.system(awk_find_nano)

# isolate the reads containing the TE

awk_extract_reads = "awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=(substr($2,1,36) in ids)} f' nanopore_reads_Paris.txt nanoporeParis.fa.masked > nanoporeParisOnly.fa.masked"

os.system(awk_extract_reads)

# isolate the TE locations in the assemblies

awk_find_TE_160 = '''awk '{if (($10) == "Paris") print $0;}' strain160all_filter.fa.out > strain160all.paris.fa.out'''
awk_find_TE_9 = '''awk '{if (($10) == "Paris") print $0;}' strain9all_filter.fa.out > strain9all.paris.fa.out'''

os.system(awk_find_TE_160)
os.system(awk_find_TE_9)

# Index the reference assembly sequences

bwa_index_160 = "bwa index strain160all.fasta.masked"
bwa_index_9 = "bwa index strain9all.fasta.masked"

# Normalize Nanaopore Reads to Assemblies

bwa_mem_160 = "bwa mem strain160all.fasta.masked nanoporeParisOnly.fa.masked > nanoporeParisNormalized160.sam"
bwa_mem_9 = "bwa mem strain9all.fasta.masked nanoporeParisOnly.fa.masked > nanoporeParisNormalized9.sam"

# filter sam files

samtools_view_160 = "samtools view -q 50 -F 0x800 -b nanoporeParisNormalized160.sam > nanoporeParisFilter160.bam"
samtools_view_9 = "samtools view -q 50 -F 0x800 -b nanoporeParisNormalized9.sam > nanoporeParisFilter9.bam"

# sort sam files for bedfile comparison

samtools_sort_160 = "samtools sort nanoporeParisFilter160.bam -o nanoporeParisSorted160.bam"
samtools_sort_9 = "samtools sort nanoporeParisFilter9.bam -o nanoporeParisSorted9.bam"

# convert bam file to a bed file

bedtools_bamtobed_160 = "bedtools bamtobed -i nanoporeParisSorted160.bam > nanoporeParis160.bed"
bedtools_bamtobed_9 = "bedtools bamtobed -i nanoporeParisSorted9.bam > nanoporeParis9.bed"

# convert assmebly files to bed files

awk_tobed_160 = '''awk '{OFS="\t"; if($6>1000) print $5,$6-1000, $7+1000,$10; else print $5,1,$7+1000,$10;}' strain160all.paris.fa.out > strain160all.paris.extended.bed'''
awk_tobed_9 = '''awk '{OFS="\t"; if($6>1000) print $5,$6-1000, $7+1000,$10; else print $5,1,$7+1000,$10;}' strain9all.paris.fa.out > strain9all.paris.extended.bed'''

# find inverse intersections between assembly and nanopore

bedtools_intersect_160 = "bedtools intersect -a nanoporeParis160.bed -b strain160all.paris.extended.bed -v > insertParis160.bed"
bedtools_intersect_9 = "bedtools intersect -a nanoporeParis9.bed -b strain9all.paris.extended.bed -v > insertParis9.bed"