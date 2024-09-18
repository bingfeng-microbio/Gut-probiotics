############Metagenome assemblies
##Assembly
megahit -t 10 -m 0.5 --min-contig-len 500 --k-step 10 --k-min 27 \
-1 /02.paired_data/$sampleid_paired_1.clean.fastq.gz \\
-2 /02.paired_data/$sampleid_paired_2.clean.fastq.gz \
-o /assembly/$sampleid

###Construct a non-redundant unigene set
##run prodigal
prodigal -i all.contigs.fa -a all.proteins.faa -d all.cds.fa -p meta
##Filter sequences whose length is less than 100bp
seqkit seq -m 100 -g all.cds.fa > all.cds.filter.fa
##The CDS were clustered to remove redundancy, and the similarity was 95%
cd-hit-est -i all.cds.filter.fa -o all.cds.filter.cdhit.fa -c 0.95 -aS 0.9 -g 1 -d 0 -G 0 -T 100 -M 0
##Extract the ID of CDS
sed -n "s/^>\(\S\+\).*$/\1/p" all.cds.filter.cdhit.fa > all.cds.filter.cdhit.id

seqtk subseq all.protein.fa all.cds.filter.cdhit.id > all.protein.filter.cdhit.fa

########Functional annotation
##ARG
diamond blastp --max-target-seqs 1 -p 100 -b12 -c1 -q all.contigs.fa -d /home/star/database/CARD3.1.0/card-3.1.0.dmnd -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq qlen --id 40 -o diamond_ARG_contigs.out --evalue 1e-5 --query-cover 80 --sensitive
##KEGG
diamond blastp --max-target-seqs 1 -p 100 -b12 -c1 -q all.protein.filter.cdhit.fa -d /home/star/database/kegg.dmnd -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq qlen -o diamond_kegg.out --evalue 1e-5 --sensitive
########Taxonomic annotation
diamond blastp -d /mnt/sda/nr/nr -q all.protein.filter.cdhit.fa -o unigene.daa -f 100 -e 1e-5 -p 256 -b12 -c1

# Phylogenetic analyse
phylophlan -i input_folder -d phylophlan --diversity high --accurate --min_num_markers 100 -f /home/user/phylophlan.cfg -o phylophlan --nproc 40 --verbose 2>&1 | tee tree.log


##c_AMP prediction
The c_AMP prediction codes can be found at https://github.com/mayuefine/c_AMPs-prediction