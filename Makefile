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
STAR_INDEX = data/interim/references/star/Mus_musculus.GRCm38.dna.primary_assembly.51bp

FIREBROWSE_RNASEQ = data/external/tcga-brca-firebrowse/rnaseqv2.normalized.txt
FIREBROWSE_CNV = data/external/tcga-brca-firebrowse/gistic.all_thresholded.by_genes.txt

TGCA_2012_PAM50 = data/external/tcga-brca-2012/BRCA.547.PAM50.SigClust.Subtypes.txt

TCGA_ILC_RNASEQ = data/external/tcga-ilc-2015/ILC_rnaseqv2_RSEM_genes_normalized_data_F2.txt
TCGA_ILC_SUBTYPES = data/external/tcga-ilc-2015/mmc9.xlsx

SHEAR_SPLINK_URL = https://ndownloader.figshare.com/articles/4765111?private_link=282f6ac0a014d81e137f
SHEAR_SPLINK_READS = data/interim/sb/shear_splink/reads

PROCESSED_FREEZE_URL = https://ndownloader.figshare.com/files/8297795?private_link=dd2515b13a5d022eba4d

################################################################################
# COMMANDS                                                                     #
################################################################################

all:

env:
	conda env create -f environment.yml

clean:
	find . -name "*.pyc" -exec rm {} \;

lint:
	flake8 --exclude=lib/,bin/,docs/conf.py .

docs:
	(cd docs && make html)

################################################################################
# PROJECT RULES                                                                #
################################################################################

build_bowtie_index: $(BOWTIE_INDEX)

build_star_index: $(STAR_INDEX)

preload_star_index: $(STAR_INDEX)
	STAR --genomeLoad LoadAndExit --genomeDir $(STAR_INDEX)
	rm Log.out Log.progress.out Aligned.out.sam
	rm -rf _STARtmp

download_shear_splink: $(SHEAR_SPLINK_READS)

download_freeze: tmp
	if [ -d "data/processed" ]; then $(error Data directory (data/processed) already exists); fi
	wget -O tmp/processed_freeze.tar.gz $(PROCESSED_FREEZE_URL)
	mkdir -p data/processed
	tar -xvf tmp/processed_freeze.tar.gz -C data/processed
	rm tmp/processed_freeze.tar.gz


shear_splink: $(BOWTIE_INDEX) $(ENSEMBL_GTF) $(SHEAR_SPLINK_READS)
	snakemake $(SNAKEMAKE_ARGS) -s pipelines/shear-splink.snake \
		--configfile ./configs/shear-splink.yml \
		--config samples=data/raw/sb/samples.txt \
				 reads_dir=data/interim/sb/shear_splink/reads \
		         interim_dir=data/interim/sb/shear_splink/subset \
				 processed_dir=data/processed/sb/shear_splink/subset \
				 log_dir=logs/sb/shear_splink/subset

shear_splink_full: $(BOWTIE_INDEX) $(ENSEMBL_GTF)
	snakemake $(SNAKEMAKE_ARGS) -s pipelines/shear-splink.snake \
		--configfile ./configs/shear-splink.yml \
		--config samples=data/raw/sb/samples.txt \
				 reads_dir=data/interim/sb/shear_splink/reads \
		         interim_dir=data/interim/sb/shear_splink/full \
				 processed_dir=data/processed/sb/shear_splink/full \
				 log_dir=logs/sb/shear_splink/full \
				 all_samples=True \
				 per_strain=False

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

tcga_brca_firebrowse: $(FIREBROWSE_RNASEQ) $(FIREBROWSE_CNV)

tcga_brca_2012: $(TGCA_2012_PAM50)

tcga_ilc_2015: $(TCGA_ILC_SUBTYPES) $(TCGA_ILC_RNASEQ)


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

$(STAR_INDEX): $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	mkdir -p $(STAR_INDEX)
	STAR --runMode genomeGenerate --genomeDir $(STAR_INDEX) \
		--genomeFastaFiles $(ENSEMBL_FASTA) \
		--sjdbGTFfile $(ENSEMBL_GTF) \
		--sjdbOverhang 50

$(FIREBROWSE_RNASEQ): tmp
	cd tmp && tar -xvf gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz
	mkdir -p $(dir $(FIREBROWSE_RNASEQ))
	python scripts/tcga_extract_counts.py --symbols tmp/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt $(FIREBROWSE_RNASEQ)
	rm tmp/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz
	rm -rf tmp/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0

$(FIREBROWSE_CNV): tmp
	wget -P tmp http://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/BRCA-TP/20160128/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz
	cd tmp && tar -xvf gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz
	mkdir -p $(dir $(FIREBROWSE_CNV))
	cp tmp/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt $(FIREBROWSE_CNV)
	rm tmp/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz
	rm -rf tmp/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/

$(TGCA_2012_PAM50): tmp
	wget -P tmp https://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt
	mkdir -p $(dir $(TGCA_2012_PAM50))
	mv tmp/BRCA.547.PAM50.SigClust.Subtypes.txt $(TGCA_2012_PAM50)

$(TCGA_ILC_RNASEQ): tmp
	wget -P tmp http://cbio.mskcc.org/cancergenomics/tcga/brca_tcga/ilc/ILC_rnaseqv2_RSEM_genes_normalized_data_F2.txt
	mkdir -p $(dir $(TCGA_ILC_RNASEQ))
	mv tmp/ILC_rnaseqv2_RSEM_genes_normalized_data_F2.txt $(TCGA_ILC_RNASEQ)

$(TCGA_ILC_SUBTYPES): tmp
	wget -P tmp http://www.cell.com/cms/attachment/2062298194/2064035165/mmc9.xlsx
	mkdir -p $(dir $(TCGA_ILC_SUBTYPES))
	mv tmp/mmc9.xlsx $(TCGA_ILC_SUBTYPES)

$(SHEAR_SPLINK_READS): tmp
	wget -O tmp/shearsplink_sb.zip $(SHEAR_SPLINK_URL)
	unzip -d $(SHEAR_SPLINK_READS) tmp/shearsplink_sb.zip
	rm tmp/shearsplink_sb.zip

tmp:
	mkdir -p tmp