import os, sys
import fileinput

input_dir = '/Users/evamorrison/Desktop/Eva/fullWorkflow'

os.chdir(input_dir)

# filter repeatmasker output

awk_filter_nano = "awk '{if (($7-$6) >= 250) print $0;}' dys_60000.fa.out | awk '{if($2<=20.0) print $0;}' > dys_60000_filter.fa.out"
awk_filter_160 = "awk '{if (($7-$6) >= 250) print $0;}' strain160all.fasta.out | awk '{if($2<=20.0) print $0;}' > strain160all_filter.fa.out"
awk_filter_9 = "awk '{if (($7-$6) >= 250) print $0;}' strain9all.fasta.out | awk '{if($2<=20.0) print $0;}' > strain9all_filter.fa.out"

os.system(awk_filter_nano)
os.system(awk_filter_160)
os.system(awk_filter_9)

#pull the list of TEs

pull_TEs = '''grep -e '>' miniTELibrary.fasta | awk 'sub(/^>/, "")' > TE_list.txt'''

os.system(pull_TEs)

for line in fileinput.input("TE_list.txt"):
    TE = line
    TE = TE.strip('\n') #strip the new line character at the end of the TE

    # Find the reads containing the TE

    awk_find_nano = '''awk '{if (($10) == "'''+TE+'''") print $0;}' dys_60000_filter.fa.out | awk '{Nanopore[$5];} END{for (var in Nanopore) print var;} ' > nanopore_reads_'''+TE+'''.txt'''
    os.system(awk_find_nano)

    # isolate the reads containing the TE

    awk_extract_reads = '''awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=(substr($2,1,36) in ids)} f' nanopore_reads_'''+TE+'''.txt dys_60000.fa.masked > dys_60000_'''+TE+'''.fa.masked'''
    remove = '''rm nanopore_reads_'''+TE+'''.txt'''

    os.system(awk_extract_reads)
    os.system(remove)

    # isolate the TE locations in the assemblies

    awk_find_TE_160 = '''awk '{if (($10) == "'''+TE+'''") print $0;}' strain160all_filter.fa.out > strain160all.'''+TE+'''.fa.out'''
    awk_find_TE_9 = '''awk '{if (($10) == "'''+TE+'''") print $0;}' strain9all_filter.fa.out > strain9all.'''+TE+'''.fa.out'''

    # print(awk_find_TE_160)
    # print(awk_find_TE_9)

    os.system(awk_find_TE_160)
    os.system(awk_find_TE_9)

    