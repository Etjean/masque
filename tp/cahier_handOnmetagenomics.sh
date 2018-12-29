# WORKDIR
WORKDIR=./
# DATADIR
DATADIR=../test/data/
# MASQUEDIR
MASQUEDIR=../


################
# 16S analysis #
################
# Even samples from MOCK communities
#sample=SRR053818
#sample=SRR072220
#sample=SRR072239
# Staggered samples from MOCK communities
#sample=SRR072221
#sample=SRR072223
#sample=SRR072237
for sample in $DATADIR/*.gz
do
    sample_name=$(basename $sample | cut -f 1 -d '.')
    # Filtering quality
    $MASQUEDIR/vsearch_bin/bin/vsearch --fastq_filter $DATADIR/${sample_name}.fastq.gz --fastqout $WORKDIR/${sample_name}_filt.fastq --fastq_truncqual 16 --fastq_trunclen 250
    # Add sample
    $MASQUEDIR/fastq2fasta/fastq2fasta.py -i $WORKDIR/${sample_name}_filt.fastq -s $sample_name -o $WORKDIR/${sample_name}_filt.fasta
done
sample="mock_community"
cat *_filt.fasta > ${sample}_filt.fasta
# Dereplication
$MASQUEDIR/vsearch_bin/bin/vsearch --derep_fulllength $WORKDIR/${sample}_filt.fasta --output $WORKDIR/${sample}_drep.fasta --sizeout --strand both
# Remove singleton
$MASQUEDIR/vsearch_bin/bin/vsearch --sortbysize $WORKDIR/${sample}_drep.fasta --output $WORKDIR/${sample}_nosing.fasta  --minsize 2
# Remove chimera
$MASQUEDIR/vsearch_bin/bin/vsearch --uchime_denovo $WORKDIR/${sample}_nosing.fasta --chimeras $WORKDIR/${sample}_chim.fasta --nonchimeras $WORKDIR/${sample}_nochim.fasta
# Clustering
$MASQUEDIR/vsearch_bin/bin/vsearch --cluster_size $WORKDIR/${sample}_nochim.fasta --id 0.97 --centroids $WORKDIR/${sample}_otu.fasta --sizein --relabel OTU_  --strand both
# Mapping
$MASQUEDIR/vsearch_bin/bin/vsearch --usearch_global $WORKDIR/${sample}_filt.fasta --db $WORKDIR/${sample}_otu.fasta --id 0.97 --strand both --otutabout $WORKDIR/${sample}_otu_table.tsv --biomout $WORKDIR/${sample}.biom
# Annotation MOCK
$MASQUEDIR/vsearch_bin/bin/vsearch --usearch_global $WORKDIR/${sample}_otu.fasta --db $DATADIR/HMP_MOCK_v35_annotated.fasta --id 0.9  --blast6out $WORKDIR/${sample}_vs_mock.tsv --alnout $WORKDIR/${sample}_vs_mock_ali.txt --strand both
$MASQUEDIR/get_taxonomy/get_taxonomy.py -i $WORKDIR/${sample}_vs_mock.tsv -d $DATADIR/HMP_MOCK_v35_annotated.fasta -u $WORKDIR/${sample}_otu.fasta -o $WORKDIR/${sample}_annotation_mock.tsv -ob $WORKDIR/${sample}_annotation_mock_biom.txt
# Biom
biom add-metadata -i $WORKDIR/${sample}.biom -o $WORKDIR/${sample}_annotated.biom --observation-metadata-fp $WORKDIR/${sample}_annotation_mock_biom.txt --observation-header id,taxonomy --sc-separated taxonomy --sample-metadata-fp sample_type.txt --output-as-json
