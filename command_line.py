import os, sys
import fileinput

input_dir = '/home/e378m007/fullWorkflow'

os.chdir(input_dir)

# filter repeatmasker output

awk_filter_nano = "awk '{if (($7-$6) >= 250) print $0;}' dys_60000.fa.out | awk '{if($2<=20.0) print $0;}' > dys_60000_filter.fa.out"
awk_filter_160 = "awk '{if (($7-$6) >= 250) print $0;}' strain160all.fasta.out | awk '{if($2<=20.0) print $0;}' > strain160all_filter.fa.out"
awk_filter_9 = "awk '{if (($7-$6) >= 250) print $0;}' strain9all.fasta.out | awk '{if($2<=20.0) print $0;}' > strain9all_filter.fa.out"

os.system(awk_filter_nano)
os.system(awk_filter_160)
os.system(awk_filter_9)

#pull the list of TEs

pull_TEs = "grep -e '>' miniTELibrary.fasta |awk 'sub(/^>/, "")' > TE_list.txt"

os.system(pull_TEs)

for line in fileinput.input("TE_list.txt"):
    TE = line

    # Find the reads containing the TE

    awk_find_nano = '''awk '{if (($10) == '''+TE+''') print $0;}' dys_60000_filter.fa.out | awk '{Nanopore[$5];} END{for (var in Nanopore) print var;} ' > nanopore_reads_'''+TE+'''.txt'''
    os.system(awk_find_nano)

    # isolate the reads containing the TE

    awk_extract_reads = "awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=(substr($2,1,36) in ids)} f' nanopore_reads_"+TE+".txt dys_60000.fa.masked > dys_60000_"+TE+".fa.masked"
    remove = "rm nanopore_reads_"+TE+".txt"

    os.system(awk_extract_reads)
    os.system(remove)

    # isolate the TE locations in the assemblies

    awk_find_TE_160 = '''awk '{if (($10) == "'''+TE+'''") print $0;}' strain160all_filter.fa.out > strain160all.'''+TE+'''.fa.out'''
    awk_find_TE_9 = '''awk '{if (($10) == "'''+TE+'''") print $0;}' strain9all_filter.fa.out > strain9all.'''+TE+'''.fa.out'''


    os.system(awk_find_TE_160)
    os.system(awk_find_TE_9)

    # Index the reference assembly sequences

    bwa_index_160 = "bwa index strain160all.fasta.masked"
    bwa_index_9 = "bwa index strain9all.fasta.masked"

    os.system(bwa_index_160)
    os.system(bwa_index_9)

    # Normalize Nanaopore Reads to Assemblies

    bwa_mem_160 = "bwa mem strain160all.fasta.masked dys_60000_"+TE+".fa.masked > dys_60000Normalized160.sam"
    bwa_mem_9 = "bwa mem strain9all.fasta.masked dys_60000_"+TE+".fa.masked > dys_60000Normalized9.sam"

    os.system(bwa_mem_160)
    os.system(bwa_mem_9)

    # filter sam files

    samtools_view_160 = "samtools view -q 50 -F 0x800 -b dys_60000_Normalized160.sam > dys_60000_Filter160.bam"
    samtools_view_9 = "samtools view -q 50 -F 0x800 -b dys_60000_Normalized9.sam > dys_60000_Filter9.bam"

    os.system(samtools_view_160)
    os.system(samtools_view_9)

    # sort sam files for bedfile comparison

    samtools_sort_160 = "samtools sort dys_60000_Filter160.bam -o dys_60000_Sorted160.bam"
    samtools_sort_9 = "samtools sort dys_60000_Filter9.bam -o dys_60000_Sorted9.bam"

    os.system(samtools_sort_160)
    os.system(samtools_sort_9)

    # convert bam file to a bed file

    bedtools_bamtobed_160 = "bedtools bamtobed -i dys_60000_Sorted160.bam > dys_60000_160.bed"
    bedtools_bamtobed_9 = "bedtools bamtobed -i dys_60000_Sorted9.bam > dys_60000_9.bed"

    os.system(bedtools_bamtobed_160)
    os.system(bedtools_bamtobed_9)

    # convert assmebly files to bed files

    awk_tobed_160 = '''awk '{OFS="\t"; if($6>1000) print $5,$6-1000, $7+1000,$10; else print $5,1,$7+1000,$10;}' strain160all.'''+TE+'''.fa.out > strain160all.'''+TE+'''.extended.bed'''
    awk_tobed_9 = '''awk '{OFS="\t"; if($6>1000) print $5,$6-1000, $7+1000,$10; else print $5,1,$7+1000,$10;}' strain9all.'''+TE+'''.fa.out > strain9all.'''+TE+'''.extended.bed'''

    os.system(awk_tobed_160)
    os.system(awk_tobed_9)

    # find inverse intersections between assembly and nanopore

    bedtools_intersect_160 = "bedtools intersect -a dys_60000_160.bed -b strain160all."+TE+".extended.bed -v > insert"+TE+"160.bed"
    bedtools_intersect_9 = "bedtools intersect -a dys_60000_9.bed -b strain9all."+TE+".extended.bed -v > insert"+TE+"9.bed"

    os.system(bedtools_intersect_160)
    os.system(bedtools_intersect_9)
