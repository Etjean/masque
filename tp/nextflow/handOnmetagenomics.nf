
MASQUEDIR = '/home/etienne/masque/'
DATADIR = "$MASQUEDIR/test/data/"
WORKDIR = "$MASQUEDIR/tp/nextflow"
EXTENSION = '.fastq.gz'
SAMPLE = 'mock_community'
fq = Channel.fromPath("$DATADIR/*$EXTENSION")
db = Channel.value(file("$DATADIR/HMP_MOCK_v35_annotated.fasta"))
sample_type = Channel.fromPath("$WORKDIR/sample_type.txt")


// Filtering quality
process qualityFiltering {

	publishDir path: WORKDIR, mode: 'copy'

	input:
		file reads from fq

	output:
		file reads_filt into fq_filt
	
	script:
		sample = reads.name.take(reads.name.lastIndexOf(EXTENSION))
		reads_filt = sample+'_filt.fastq'
		"""
		$MASQUEDIR/vsearch_bin/bin/vsearch \
			--fastq_filter $reads \
			--fastqout $reads_filt \
			--fastq_truncqual 16 \
			--fastq_trunclen 250
		"""
}

// Add sample
process fastq2fasta {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file reads_filt from fq_filt
		
	output:
		file reads_fasta into fa_filt
		
	script:
		sample = reads_filt.name.take(reads_filt.name.lastIndexOf('_filt.fastq'))
		reads_fasta = sample+'_filt.fasta'
		"""
		$MASQUEDIR/fastq2fasta/fastq2fasta.py \
		-i $reads_filt \
		-s $sample \
		-o $reads_fasta
		"""
}

// Concatenate files
fa_filt
	.collectFile(name: SAMPLE+'_filt.fasta', storeDir: WORKDIR)
	.into {mock_filt1; mock_filt2}

// Dereplication
process dereplication {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file mock_filt1
	
	output:
		file SAMPLE+'_drep.fasta' into mock_drep
	
	script:
		"""
		$MASQUEDIR/vsearch_bin/bin/vsearch \
		--derep_fulllength $mock_filt1 \
		--output ${SAMPLE}_drep.fasta \
		--sizeout \
		--strand both
		"""
}

// Remove singleton
process removeSingleton {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file mock_drep
	
	output:
		file SAMPLE+'_nosing.fasta' into mock_nosing
	
	script:
		"""
		$MASQUEDIR/vsearch_bin/bin/vsearch \
		--sortbysize $mock_drep \
		--output ${SAMPLE}_nosing.fasta  \
		--minsize 2
		"""
}

// Remove chimera
process removeChimera {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file mock_nosing
	
	output:
		file SAMPLE+'_nochim.fasta' into mock_nochim
		file SAMPLE+'_chim.fasta'
	
	script:
		"""
		$MASQUEDIR/vsearch_bin/bin/vsearch \
		--uchime_denovo $mock_nosing \
		--chimeras ${SAMPLE}_chim.fasta \
		--nonchimeras ${SAMPLE}_nochim.fasta
		"""
}

// Clustering
process clustering {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file mock_nochim
	
	output:
		file SAMPLE+'_otu.fasta' into mock_otu
	
	script:
		"""
		$MASQUEDIR/vsearch_bin/bin/vsearch \
		--cluster_size $mock_nochim \
		--id 0.97 \
		--centroids ${SAMPLE}_otu.fasta \
		--sizein \
		--relabel OTU_  \
		--strand both
		"""
}

mock_otu.into {mock_otu1; mock_otu2; mock_otu3}

// Mapping
process mapping {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file mock_otu1
		file mock_filt2
		
	output:
		file SAMPLE+'_otu_table.tsv'
		file SAMPLE+'.biom' into mock_biom
	
	"""
	$MASQUEDIR/vsearch_bin/bin/vsearch \
	--usearch_global $mock_filt2 \
	--db $mock_otu1 \
	--id 0.97 \
	--strand both \
	--otutabout ${SAMPLE}_otu_table.tsv \
	--biomout ${SAMPLE}.biom
	"""
}

// Annotation MOCK
process annotation {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file mock_otu2
		file db
	
	output:
		file SAMPLE+'_vs_mock.tsv' into mock_vs_mock
		file SAMPLE+'_vs_mock_ali.txt'
	
	"""
	$MASQUEDIR/vsearch_bin/bin/vsearch \
	--usearch_global $mock_otu2 \
	--db $db \
	--id 0.9  \
	--blast6out ${SAMPLE}_vs_mock.tsv \
	--alnout ${SAMPLE}_vs_mock_ali.txt \
	--strand both
	"""
}

// Get taxonomy
process getTaxonomy {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file mock_vs_mock
		file db
		file mock_otu3
	
	output:
		file SAMPLE+'_annotation_mock.tsv'
		file SAMPLE+'_annotation_mock_biom.txt' into mock_annotation_biom
	
	"""
	$MASQUEDIR/get_taxonomy/get_taxonomy.py \
	-i $mock_vs_mock \
	-d $db \
	-u $mock_otu3 \
	-o ${SAMPLE}_annotation_mock.tsv \
	-ob ${SAMPLE}_annotation_mock_biom.txt
	"""
}

// Biom
process biom {
	
	publishDir path: WORKDIR, mode: 'copy'
	
	input:
		file mock_biom
		file mock_annotation_biom
		file sample_type
	
	output:
		file SAMPLE+'_annotated.biom'
		
	
	"""
	biom add-metadata \
	-i $mock_biom \
	-o ${SAMPLE}_annotated.biom \
	--observation-metadata-fp $mock_annotation_biom \
	--observation-header id,taxonomy \
	--sc-separated taxonomy \
	--sample-metadata-fp $sample_type \
	--output-as-json
	"""
}


