.PHONY: clean data lint env

################################################################################
# GLOBALS                                                                      #
################################################################################

# Snakemake
SNAKEMAKE_ARGS =

ENSEMBL_FASTA = data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa
ENSEMBL_GTF =  data/external/ensembl/Mus_musculus.GRCm38.76.gtf

ENSEMBL_FASTA_URL = ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
ENSEMBL_GTF_URL = ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/Mus_musculus.GRCm38.76.gtf.gz

BOWTIE_INDEX = data/interim/references/bowtie/Mus_musculus.GRCm38.dna.primary_assembly.1.bt2
TOPHAT_INDEX = data/interim/references/tophat/Mus_musculus.GRCm38.76.t2onc.tophat2.1.bt2
STAR_INDEX = data/interim/star/Mus_musculus.GRCm38.dna.primary_assembly.51bp


################################################################################
# COMMANDS                                                                     #
################################################################################

all:

env:
	conda env create -f envs/sb-screen.yml

clean:
	find . -name "*.pyc" -exec rm {} \;

lint:
	flake8 --exclude=lib/,bin/,docs/conf.py .


################################################################################
# PROJECT RULES                                                                #
################################################################################

build_tophat_index: $(TOPHAT_INDEX)

build_bowtie_index: $(BOWTIE_INDEX)

build_star_index: $(STAR_INDEX)

preload_star_index: $(STAR_INDEX)
	STAR --genomeLoad LoadAndExit --genomeDir $(STAR_INDEX)
	rm Log.out Log.progress.out Aligned.out.sam
	rm -rf _STARtmp

shear_splink: $(BOWTIE_INDEX) $(ENSEMBL_GTF)
	snakemake $(SNAKEMAKE_ARGS) -s pipelines/shear-splink.snake \
		--configfile ./configs/shear-splink.yml \
		--config samples=data/raw/sb/samples.txt \
				 raw_dir=data/raw \
		         interim_dir=data/interim/sb/shear_splink \
				 processed_dir=data/processed/sb/shear_splink \
				 log_dir=logs/sb/shear_splink

rnaseq_sb: $(STAR_INDEX) $(ENSEMBL_GTF)
	snakemake $(SNAKEMAKE_ARGS) -s pipelines/rnaseq.snake \
		--configfile ./configs/rnaseq.yml \
		--config samples=data/raw/sb/samples.rnaseq.txt \
		         interim_dir=data/interim/sb/rnaseq \
				 processed_dir=data/processed/sb/rnaseq \
				 log_dir=logs/sb/rnaseq \
				 qc_dir=qc/sb/rnaseq

rnaseq_kb1p: $(STAR_INDEX) $(ENSEMBL_GTF)
	snakemake $(SNAKEMAKE_ARGS) -s pipelines/rnaseq.snake \
		--configfile ./configs/rnaseq.yml \
		--config samples=data/raw/kb1p/samples.rnaseq.txt \
		         interim_dir=data/interim/kb1p \
				 processed_dir=data/processed/kb1p \
				 log_dir=logs/kb1p \
				 qc_dir=qc/kb1p

rnaseq_pten: $(STAR_INDEX) $(ENSEMBL_GTF)
	snakemake $(SNAKEMAKE_ARGS) -s pipelines/rnaseq.snake \
		--configfile ./configs/rnaseq.yml \
		--config samples=data/raw/pten/samples.rnaseq.txt \
		         interim_dir=data/interim/pten \
				 processed_dir=data/processed/pten \
				 log_dir=logs/pten \
				 qc_dir=qc/pten

nmf:
	Rscript --vanilla scripts/nmf.R \
		--counts data/processed/sb/rnaseq/gene_counts.txt \
		--norm_counts data/processed/sb/rnaseq/gene_counts.log2.txt \
		--output_dir data/processed/sb/nmf \
		--cores 20


################################################################################
# FILE RULES                                                                   #
################################################################################

$(ENSEMBL_FASTA):
	wget --quiet -P $(dir $(ENSEMBL_FASTA)) $(ENSEMBL_FASTA_URL)
	gunzip $(ENSEMBL_FASTA).gz

$(ENSEMBL_GTF):
	wget --quiet -P $(dir $(ENSEMBL_GTF)) $(ENSEMBL_GTF_URL)
	gunzip $(ENSEMBL_GTF).gz
	(grep ^"#" $(ENSEMBL_GTF); grep -v ^"#" $(ENSEMBL_GTF)| sort -k1,1 -k4,4n) \
		| bgzip > $(ENSEMBL_GTF).gz;
	tabix -p gff $(ENSEMBL_GTF).gz

$(BOWTIE_INDEX): $(ENSEMBL_FASTA)
	bowtie2-build $(ENSEMBL_FASTA) $(basename $(basename $(BOWTIE_INDEX)))

$(TOPHAT_INDEX): $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	mkdir -p $(dir $(TOPHAT_INDEX))
	bowtie2-build $(ENSEMBL_FASTA) $(basename $(basename $(TOPHAT_INDEX)))
	tophat2 -G $(ENSEMBL_GTF) \
		--transcriptome-index=$(basename $(basename $(TOPHAT_INDEX))).transcriptome \
		$(basename $(basename $(TOPHAT_INDEX)))
	rm -rf ./tophat_out

$(STAR_INDEX): $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	mkdir -p $(STAR_INDEX)
	STAR --runMode genomeGenerate --genomeDir $(STAR_INDEX) \
		--genomeFastaFiles $(ENSEMBL_FASTA) \
		--sjdbGTFfile $(ENSEMBL_GTF) \
		--sjdbOverhang 50
