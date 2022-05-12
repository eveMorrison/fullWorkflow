import os, sys
import fileinput

input_dir = '/home/e378m007/fullWorkflow'

os.chdir(input_dir)

# filter repeatmasker output

awk_filter_nano = "awk '{if (($7-$6) >= 250) print $0;}' first15all.fa.out | awk '{if($2<=20.0) print $0;}' > first15all_filter.fa.out"
awk_filter_160 = "awk '{if (($7-$6) >= 250) print $0;}' strain160all.fasta.out | awk '{if($2<=20.0) print $0;}' > strain160all_filter.fa.out"
awk_filter_9 = "awk '{if (($7-$6) >= 250) print $0;}' strain9all.fasta.out | awk '{if($2<=20.0) print $0;}' > strain9all_filter.fa.out"

os.system(awk_filter_nano)
os.system(awk_filter_160)
os.system(awk_filter_9)

#pull the list of TEs

pull_TEs = '''grep -e '>' refTElib.fasta | awk 'sub(/^>/, "")' > TE_list.txt'''

os.system(pull_TEs)

# Index the reference assembly sequences

bwa_index_160 = '''bwa index strain160all.fasta.masked'''
bwa_index_9 = '''bwa index strain9all.fasta.masked'''

os.system(bwa_index_160)
os.system(bwa_index_9)

for line in fileinput.input("TE_list.txt"):
    TE = line
    TE = TE.strip('\n') #strip the new line character at the end of the TE

    # Find the reads containing the TE

    awk_find_nano = '''awk '{if (($10) == "'''+TE+'''") print $0;}' first15all_filter.fa.out | awk '{Nanopore[$5];} END{for (var in Nanopore) print var;} ' > nanopore_reads_'''+TE+'''.txt'''
    os.system(awk_find_nano)

    # isolate the reads containing the TE

    awk_extract_reads = '''awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=(substr($2,1,36) in ids)} f' nanopore_reads_'''+TE+'''.txt first15all.fa.masked > first15all_'''+TE+'''.fa.masked'''
    remove = '''rm nanopore_reads_'''+TE+'''.txt'''

    os.system(awk_extract_reads)
    os.system(remove)

    # isolate the TE locations in the assemblies

    awk_find_TE_160 = '''awk '{if (($10) == "'''+TE+'''") print $0;}' strain160all_filter.fa.out > strain160all.'''+TE+'''.fa.out'''
    awk_find_TE_9 = '''awk '{if (($10) == "'''+TE+'''") print $0;}' strain9all_filter.fa.out > strain9all.'''+TE+'''.fa.out'''


    os.system(awk_find_TE_160)
    os.system(awk_find_TE_9)

    # Normalize Nanaopore Reads to Assemblies

    bwa_mem_160 = '''bwa mem strain160all.fasta.masked first15all_'''+TE+'''.fa.masked > first15all_Normalized160_'''+TE+'''.sam'''
    bwa_mem_9 = '''bwa mem strain9all.fasta.masked first15all_'''+TE+'''.fa.masked > first15all_Normalized9_'''+TE+'''.sam'''
    remove = '''rm first15all_'''+TE+'''.fa.masked'''

    os.system(bwa_mem_160)
    os.system(bwa_mem_9)
    os.system(remove)

    # filter sam files

    samtools_view_160 = '''samtools view -q 50 -F 0x800 -b first15all_Normalized160_'''+TE+'''.sam > first15all_Filter160_'''+TE+'''.bam'''
    samtools_view_9 = '''samtools view -q 50 -F 0x800 -b first15all_Normalized9_'''+TE+'''.sam > first15all_Filter9_'''+TE+'''.bam'''
    remove_160 = '''rm first15all_Normalized160_'''+TE+'''.sam'''
    remove_9 = '''rm first15all_Normalized9_'''+TE+'''.sam'''

    os.system(samtools_view_160)
    os.system(samtools_view_9)
    os.system(remove_160)
    os.system(remove_9)

    # sort sam files for bedfile comparison

    samtools_sort_160 = '''samtools sort first15all_Filter160_'''+TE+'''.bam -o first15all_Sorted160_'''+TE+'''.bam'''
    samtools_sort_9 = '''samtools sort first15all_Filter9_'''+TE+'''.bam -o first15all_Sorted9_'''+TE+'''.bam'''
    remove_160 = '''rm first15all_Filter160_'''+TE+'''.bam'''
    remove_9 = '''rm first15all_Filter9_'''+TE+'''.bam'''

    os.system(samtools_sort_160)
    os.system(samtools_sort_9)
    os.system(remove_160)
    os.system(remove_9)

    # convert bam file to a bed file

    bedtools_bamtobed_160 = '''bedtools bamtobed -i first15all_Sorted160_'''+TE+'''.bam > first15all_160_'''+TE+'''.bed'''
    bedtools_bamtobed_9 = '''bedtools bamtobed -i first15all_Sorted9_'''+TE+'''.bam > first15all_9_'''+TE+'''.bed'''
    remove_160 = '''rm first15all_Sorted160_'''+TE+'''.bam'''
    remove_9 = '''rm first15all_Sorted9_'''+TE+'''.bam'''

    os.system(bedtools_bamtobed_160)
    os.system(bedtools_bamtobed_9)
    os.system(remove_160)
    os.system(remove_9)

    # convert assmebly files to bed files

    awk_tobed_160 = '''awk '{OFS="\t"; if($6>1000) print $5,$6-1000, $7+1000,$10; else print $5,1,$7+1000,$10;}' strain160all.'''+TE+'''.fa.out > strain160all.'''+TE+'''.extended.bed'''
    awk_tobed_9 = '''awk '{OFS="\t"; if($6>1000) print $5,$6-1000, $7+1000,$10; else print $5,1,$7+1000,$10;}' strain9all.'''+TE+'''.fa.out > strain9all.'''+TE+'''.extended.bed'''
    remove_160 = '''rm strain160all.'''+TE+'''.fa.out'''
    remove_9 = '''rm strain9all.'''+TE+'''.fa.out'''

    os.system(awk_tobed_160)
    os.system(awk_tobed_9)
    os.system(remove_160)
    os.system(remove_9)

    # find inverse intersections between assembly and nanopore

    bedtools_intersect_160 = '''bedtools intersect -a first15all_160_'''+TE+'''.bed -b strain160all.'''+TE+'''.extended.bed -v > insert_'''+TE+'''_160.bed'''
    bedtools_intersect_9 = '''bedtools intersect -a first15all_9_'''+TE+'''.bed -b strain9all.'''+TE+'''.extended.bed -v > insert_'''+TE+'''_9.bed'''
    remove_160 = '''rm first15all_160_'''+TE+'''.bed'''
    remove_9 = '''rm first15all_9_'''+TE+'''.bed'''
    remove_160_extend = '''rm strain160all.'''+TE+'''.extended.bed'''
    remove_9_extend = '''rm strain9all.'''+TE+'''.extended.bed'''

    os.system(bedtools_intersect_160)
    os.system(bedtools_intersect_9)
    os.system(remove_160)
    os.system(remove_9)
    os.system(remove_160_extend)
    os.system(remove_9_extend)
