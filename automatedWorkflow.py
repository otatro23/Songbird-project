import subprocess

def seqReadAssessment(rawReads):
    numRawReads = subprocess.check_output(['zgrep', '-c', '^>', rawReads]) 
    numRawReads = int(str(numRawReads)[2:4]) *2 #fixes numRawReads being b'24\n' and accoutns for the forward and backwards read
    totalbp = numRawReads * 250 * 2
    coverage = totalbp / 5000000 #approx.
    print(numRawReads)

def examineReadQuality(rawReadsR1, rawReadsR2):
    subprocess.call(['mkdir', 'fastqc_raw-reads'], cwd='/home/tatroo23/Songbird/Songbird-project')
    subprocess.call(['fastqc', 'rawReadsR1' 'rawReadsR2', '-o', 'fastqc_raw-reads'])

def trim(rawReadsR1, rawReadsR2):
    subprocess.call(['trim_scriptV2.sh', 'rawReadsR1', 'rawReadsR2'], cwd='/home/tatroo23/Songbird/Songbird-project')
    subprocess.call(['mkdir', 'trimmed_reads'])
    subprocess.call(['mv *fastq.gz', 'trimmed_reads/'])

def genomeAssembly():
    subprocess.call(['spades.py', '-1', 'trimmed_reads/paired_forward.fastq.gz', '-2 trimmed_reads/paired_reverse.fastq.gz', '-s', 'trimmed_reads/unpaired_forward.fastq.gz', '-s', 'trimmed_reads/unpaired_reverse.fastq.gz', '-o', 'spades_assembly_default'])
    subprocess.check_output(['grep', '-c', '^>', 'spades_assembly_default/contigs.fasta']) 
    subprocess.call(['trim_scriptV2.sh', 'rawReadsR1', 'rawReadsR2'], cwd='/home/tatroo23/Songbird/Songbird-project/spades_spades_assembly_default')
    subprocess.call(['mv', 'contigs.fasta', 'spades.log', '../'])
    subprocess.call(['rm', '-r', '*'])
    subprocess.call(['mv', '../contigs.fasta', '../spades.log', './'])
    subprocess.call(['quast.py', 'contigs.fasta', '-o', 'quast_results'])
    subprocess.call(['busco', '-i', 'contigs.fasta', '-m', 'genome', '-o', 'busco-results', '-l', 'bacteria'])
    subprocess.call(['prokka', 'contigs.fasta', '--outdir', 'prokka_output', '--cpus', '24', '--mincontiglen', '200'])
    subprocess.call(['grep', '-o', '"product=.*"', 'prokka_output/PROKKA_*.gff', '|', 'sed', ''s/product=//g'', '|', 'sort', '|', 'uniq', '-c' '|', 'sort', '-nr', '>', 'protein_abundances.txt'])
    subprocess.call(['extract_sequences', '"16S ribosomal RNA"', 'prokka_output/PROKKA_*.ffn', '>', '16S_sequence.fasta'])
    subprocess.call(['makeblastdb', '-in', 'contigs.fasta', '-dbtype', 'nucl', '-out', 'contigs_db'])
    subprocess.call(['blastn', '-query', '16S_sequence.fasta', '-db', 'contigs_db', '-out', '16S_vs_contigs_6.tsv', '-outfmt', '6'])
    subprocess.call(['blob_blast.sh', 'contigs.fasta'])
    subprocess.call(['bwa', 'index', '$fasta'])
    subprocess.call(['bwa', 'mem', '-t', '24', '$fasta', '$forward', '$reverse', '>', 'raw_mapped.sam'])
    subprocess.call(['samtools', 'view', '-@', '24', '-Sb',  'raw_mapped.sam',  '|', 'samtools sort', '-@', '24', '-', 'sorted_mapped'])
    subprocess.call(['samtools', 'flagstat', 'sorted_mapped.bam'])
    subprocess.call(['samtools', 'index', 'sorted_mapped.bam'])
    subprocess.call(['bedtools', 'genomecov', '-ibam', 'sorted_mapped.bam', '>', 'coverage.out'])
    subprocess.call(['gen_input_table.py',  '--isbedfiles', '$fasta', 'coverage.out', '>',  'coverage_table.tsv'])
    subprocess.call(['blobtools', 'create', '-i', 'contigs.fasta', '-b', 'sorted_mapped.bam', '-t', 'contigs.fasta.vs.nt.cul5.1e5.megablast.out', '-o', 'blob_out'])
    subprocess.call(['blobtools', 'view', '-i', 'blob_out.blobDB.json', '-r', 'all', '-o', 'blob_taxonomy'])
    subprocess.call(['blobtools', 'plot', '-i', 'blob_out.blobDB.json', '-r', 'genus'])
    subprocess.call(['mkdir', '~/mdibl-t3-2018-WGS/filtered_assembly'])
    subprocess.call([''], cwd='/home/tatroo23/Songbird/Songbird-project/spades_spades_assembly_default/~/mdibl-t3-2018-WGS/filtered_assembly')
    subprocess.call(['cp', '~/~/mdibl-t3-2018-WGS/blob_taxonomy.blob_out.blobDB.table.txt' './')]
    subprocess.call(['grep', '-v', 'blob_taxonomy.blob_out.blobDB.table.txt', '|', 'awk', "-F'\t'", "'$2 > 500'", '|', 'awk', "-F'\t'", "'$5 > 20'", '|', 'awk', "-F'\t'" "'{print $1}'" '>', 'list_of_contigs_to_keep_len500_cov20.txt'])
    subprocess.call(['filter_contigs_by_list.py', '~/mdibl-t3-WGS/spades_assembly/contigs.fasta list_of_contigs_to_keep_len500_cov20.txt', 'Streptomyces_A1277_filtered.fasta'])
    subprocess.call(['grep', '-f', 'list_to_keep.txt', 'blob_taxonomy.blob_out.blobDB.table.txt', '|', 'awk', '{w = w + $2; e = e + $5 * $2;} END {print e/w}'])
    subprocess.call(['wget', '"https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"')
    subprocess.call(['blastn', '-reward', '1', 'penalty', '-5', '-gapopen', '3', 'gapextend', '3', '-dust', 'yes', '-soft_masking', 'true', '-evalue', '700', 'searchsp', '1750000000000', '-query', 'filtered-SKBO6.fasta', '-subject', 'UniVec', '-outfmt', '6', '-out', 'genome_vs_univec.6'])


    

seqReadAssessment('/home/tatroo23/Songbird/Songbird-project/data.fna')
