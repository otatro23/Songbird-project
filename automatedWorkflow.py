import subprocess
import time

startTime = time.time()
R1 = "PDR3_S3_L001_R1_001.fastq.gz"
R2 = "PDR3_S3_L001_R2_001.fastq.gz"


# Finding the number of raw reads, total bp data, and average coverage 
numRawReads = subprocess.check_output(["zgrep", "-c", "@A01", "raw_reads/"+R1])
numRawReads = int(numRawReads.strip())
totalbp = numRawReads*250*2
coverage = totalbp/numRawReads



# Exam raw reads quality with Fastqc
subprocess.call(["mkdir", "fastqc_raw-reads"])
subprocess.call(["fastqc", "raw_reads/"+R1, "raw_reads/"+R2, "-o", "fastqc_raw-reads"])


# Adapter and quality trimming with Trimmomatic
subprocess.call(["trim_scriptV2.sh", R1, R2])


# Examine trimmed read quality with Fastqc
subprocess.call(["mkdir", "fastqc_trimmed-reads"])
subprocess.call(["fastqc", "trimmed-reads/"+R1, "trimmed-reads/"+R2, "trimmed-reads/unpaired-"+R1, "trimmed-reads/unpaired-"+R2, "-o", "fastqc_trimmed-reads"])
numTrimmedReads = subprocess.check_output(["zgrep", "-c", "A01", "trimmed-reads/"+R1, "trimmed-reads/"+R2, "trimmed-reads/unpaired-"+R1, "trimmed-reads/unpaired-"+R2])


# Genome Assembly with SPAdes
subprocess.call("spades.py -1 trimmed-reads/"+R1+" -2 trimmed-reads/"+R2+" -s trimmed-reads/unpaired-"+R1+ "-s trimmed-reads/unpaired-"+R2+" -o spades_assembly_default -t 24", shell=True)
numContigs = subprocess.check_output(["grep", "-c", ">", "spades_assembly_default/contigs.fasta"])
numContigs = int(numContigs.strip())


# Report data
endTime = time.time()
elapsedTime = endTime - startTime
print(" Raw Reads: " + str(numRawReads) + "\n Total BP data: " + str(totalbp) + "\n Average coverage: " + str(coverage) + "\n Contigs: " + str(numContigs) + "\n Run time: " + str(round(elapsedTime)) + " seconds")
