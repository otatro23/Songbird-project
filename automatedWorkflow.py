import subprocess
import time

startTime = time.time()

R1 = "PDR3_S3_L001_R1_001.fastq.gz"
R2 = "PDR3_S3_L001_R2_001.fastq.gz"


#Finding the number of raw reads, total bp data, and average coverage 
subprocess.call(["conda", "init", "--all"])
subprocess.call(["conda", "activate", "genomics"])
numRawReads = subprocess.check_output(["zgrep", "-c", "@A01", "raw_reads/"+R1])
numRawReads = int(numRawReads.strip())
totalbp = numRawReads*250*2
coverage = totalbp/numRawReads


#Exam raw reads quality with Fastqc
subprocess.call(["mkdir", "fastqc_raw-reads"])
subprocess.call(["fastqc", "raw_reads/"+R1, "raw_reads/"+R2, "-o", "fastqc_raw-reads"])


#Adapter and quality trimming with Trimmomatic
subprocess.call(["trim_scriptV2.sh", R1, R2])


#Examine trimmed read quality with Fastqc
subprocess.call(["mkdir", "fastqc_trimmed-reads"])
subprocess.call(["fastqc", "trimmed-reads/"+R1, "trimmed-reads/"+R2, "trimmed-reads/unpaired-"+R1, "trimmed-reads/unpaired-"+R2, "-o", "fastqc_trimmed-reads"])
numTrimmedReads = subprocess.check_output(["zgrep", "-c", "A01", "trimmed-reads/"+R1, "trimmed-reads/"+R2, "trimmed-reads/unpaired-"+R1, "trimmed-reads/unpaired-"+R2])


# Genome Assembly with SPAdes
subprocess.call("spades.py -1 trimmed-reads/"+R1+" -2 trimmed-reads/"+R2+" -s trimmed-reads/unpaired-"+R1+ "-s trimmed-reads/unpaired-"+R2+" -o spades_assembly_default -t 24", shell=True)
numContigs = subprocess.check_output(["grep", "-c", ">", "spades_assembly_default/contigs.fasta"])
numContigs = int(numContigs.strip())

#Genome structure assessment with QUAST
subprocess.call(["quast.py", "spades_assembly_default/contigs.fasta", "-o", "quast_results"])


# Genome content Assessment with BUSCO
subprocess.call("busco -i spades_assembly_default/contigs.fasta -m genome -o busco-results -l bacteria", shell=True)


# Genome Annotation with PROKKA
subprocess.call("prokka spades_assembly_default/contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200", shell=True)
subprocess.call("grep -o product=.* prokka_output/PROKKA_*.gff | sed s/product=//g | sort | uniq -c | sort -nr > protein_abundances.txt", shell=True)
sixteenS = subprocess.check_output("grep 16S prokka_output/*.ffn", shell=True)
sixteenS = str(sixteenS.strip())


# Organism ID with BLAST
subprocess.call("makeblastdb -in spades_assembly_default/contigs.fasta -dbtype nucl -out contigs_db", shell=True)
subprocess.call("blob_blast.sh spades_assembly_default/contigs.fasta", shell=True)


# Read Mapping with BWA and samtools
subprocess.call("bwa index spades_assembly_default/contigs.fasta", shell=True)
subprocess.call("bwa mem -t 24 spades_assembly_default/contigs.fasta trimmed-reads/" +R1+ " trimmed-reads/" +R2+ "  > raw_mapped.sam", shell=True)
subprocess.call("samtools view -@ 24 -Sb raw_mapped.sam | samtools sort -@ 24 - sorted_mapped", shell=True)
subprocess.call("samtools flagstat sorted_mapped.bam", shell=True)

subprocess.call("samtools flagstat sorted_mapped.bam", shell=True)
subprocess.call("samtools index sorted_mapped.bam", shell=True)
subprocess.call("bedtools genomecov -ibam sorted_mapped.bam > coverage.out", shell=True)
subprocess.call("gen_input_table.py --isbedfiles spades_assembly_default/contigs.fasta coverage.out > coverage_table.tsv", shell=True)


# Non target contig removal
subprocess.call("blobtools create -i spades_assembly_default/contigs.fasta -b sorted_mapped.bam -t contigs.fasta.vs.nt.cul5.1e5.megablast.out -o blob_out", shell=True)
subprocess.call("blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy", shell=True)
subprocess.call("blobtools plot -i blob_out.blobDB.json -r genus", shell=True)


# Filter the genome assembly
subprocess.call("mkdir filtered_assembly", shell=True)
subprocess.call("cp blob_taxonomy.blob_out.blobDB.table.txt filtered_assembly", shell=True)
subprocess.call("grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F '\t' '$2>500' | awk -F '\t' '$5 > 20' > filtered_assembly/list_of_contigs_to_keep_len500_cov20.txt", shell=True)
subprocess.call("filter_contigs_by_list.py spades_assembly_default/contigs.fasta filtered_assembly/list_of_contigs_to_keep_len500_cov20.txt filtered_assembly/filtered.fasta", shell=True)


# Report data
endTime = time.time()
elapsedTime = endTime - startTime
print(" Raw Reads: " + str(numRawReads) + "\n Total BP data: " + str(totalbp) + "\n Average coverage: " + str(coverage) + "\n Contigs: " + str(numContigs) + "\n 16S: " + str(sixteenS) + "\n Run time: " + str(round(elapsedTime)) + " seconds") 
                

                
                
