# Download APA reproducible package
```
mkdir -p $HOME/projects
cd $HOME/projects
git clone https://github.com/hwanglab/apa_2019.git
```
- In this guideline, the package base directory is "$HOME/projects/apa_atingLab2019". In the following steps, it is necessary to modify scripts to be matched with your local setup.
 
# Prerequisite
1. python 2.7 
	```
	cd apa_atingLab2019
	pip install --user -r python_requirements.txt #Alternatively, use conda environment to install packages listed in the requirements.txt
	```

1. R 3.4.3 packages
	```
	handy (refer to https://github.com/jeffbhasin/handy for installation)
	goldmine
	argparse
	BSgenome.Hsapiens.UCSC.hg19
	corrplot
	data.table
	DEXSeq
	e1071
	edgeR
	ggbio
	ggplot2
	gridExtra
	IlluminaHumanMethylation450kanno.ilmn12.hg19
	matrixStats
	openxlsx
	parallel
	reshape
	RPMM
	Rsamtools
	Rsubread
	stringr
	ggseqlogo
	JASPAR2018
	TFBSTools
	UpSetR
	VennDiagram
	wateRmelon
	```
1. programs (the binary files should be in $PATH)
	```
	bowtie2
	mosdepth
	parallel
	picard
	samtools
	sra-tools
	STAR
	twoBitToFa
	wget
	wigToBigWig
	xz
	```
1. add R library path to your bash shell environment file ($HOME/.bashrc or $HOME/.bash_profile),
	```
	export R_UTIL_APA=$HOME/projects/apa_atingLab2019/resource/libs
	export PADD_GIT=$HOME/projects/apa_atingLab2019/01_polyAseq/paddle-git
	``` 
# PolyA-seq data processing

1. Configure program paths and resource files
	```bash
	bash ./01_resource_prep.sh
	```
1. Download PolyA-seq FASTQ files.  Note that, as of 05/21/2019, the SRAs, 
	- SRP083252 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86178)
	- SRP083254 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86180) 
	
	is not public. Secure tokens are available only reviewers to access processed data. Any raw data is not available yet until it becomes public.
	```
	cd 01_polyAseq
	bash ./02_get_fastq.sh
	```

1. Open the polyA-seq pipeline configuration file (`./configs/program.conf`) and edit the following two section aligned with your local setup,
	```bash
	[bowtie2-build]
	bin_path = bowtie2-build
	ref_hg19_path = ~/projects/apa_atingLab2019/resource/ref/hg19_PhiX.fa
	ref_hg19_prefix = ~/projects/apa_atingLab2019/resource/ref/hg19_PhiX
	chrom_size = ~/projects/apa_atingLab2019/resource/ref/ChromSizes.hg19.txt
	
	[ensembl]
	biotype_table = ~/projects/apa_atingLab2019/resource/ensembl/ensembl_gene_biotypes.csv
	gtf_file = ~/projects/apa_atingLab2019/resource/ensembl/Homo_sapiens.GRCh37.87.gtf
	```

1. Open the script (`03_polyaseq.sh`) and modify the package base directory  
	```
	export PPD=$HOME/projects/apa_atingLab2019/01_polyAseq
	```
1. Running polyA-seq processing pipeline
	```
	bash ./03_polyaseq.sh
	```
1. Transcription factor binding analysis in the genomic regions in between polyA sites at each 546 APA gene
	```
	Rscript ./04_enrich_perms_100k_select.r
	```
1. Plot sequence logos for the highly enriched TFs
	```
	Rscript ./05_plot_seqlogos.r
	```

# RNA-seq data processing  
1. Contact authors to download RNA-seq FASTQ files.
	
	```
	cd apa_atingLab2019/02_mRNAseq
	bash ./01_get_fastq.sh
	```

1. Run STAR 
	```
	bash ./02_run_star.sh
	```

1. Collect a basic alignment statistics 
	```
	Rscript ./03_get_alignment_stats.r
	```

1. Analyzing the APA regulatory gene expressions. If the number of CPU's available is 8, 
	```
	bash ./04_apafactorexp.sh 8
	```
	
# ChIP-seq data processing

1. We use ENCODE TF and Histone ChIP-Seq processing pipeine from Kundaje lab (e.g., chipseq.bds.20180726_175025_830). If you don't have the pipeline, then, 
	1. Download the pipeline at https://github.com/kundajelab/chipseq_pipeline
	1. Follow the installation instruction. The installation directories used in this guideline are,
		- `$HOME/apps/chipseq/chip-seq-pipeline2` # ENCODE pipeline path
		- `$HOME/.bds` # bds path
		- `$HOME/apps/miniconda3/bin` #miniconda3 path
		- `/opt/bds_pipeline_genome_data/aquas_chipseq_species.conf` #species configuration file
		- Make sure that the installation is successful
		
	1. Come back to the APA project directory, change to ChIP-seq working directory, and download ChIP/MBD-seq FASTQ files.
	
	- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131606
	- https://www.ncbi.nlm.nih.gov/sra/SRP001414
	 
	 As of 05/21/2019, note that ChIP-seq raw fastq files are not in public.
		```
		cd $HOME/projects/apa_atingLab2019/03_chipseq
		bash ./00_get_fastq.sh
		``` 
	1. Open the bash script (`chipseq_pipeline_arg.sh`) to modify the paths matched with your local setup.
		```
		progd="$HOME/apps/chipseqs/TF_chipseq_pipeline"
		species_conf="/opt/bds_pipeline_genome_data/aquas_chipseq_species.conf"
		
		export PATH=$HOME/.local/bin
		export PATH=$HOME/.bds:$PATH
		export PATH=$HOME/apps/miniconda3/bin:$PATH
		export PATH=/usr/lib64/ccache:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:$PATH
		```

1. Launch the ENCODE ChIP-seq pipeline. Note that it will take for a while to complete the job.
	```
	bash ./01_chipseq_by_encode.sh 2>&1 | tee logs/01_chipseq_by_encode_$(date '+%Y-%m-%d-%H').log
	``` 

1. Process MBD-seq FASTQ files.
	```
	bash ./02_mbdseq_by_bt2.sh 2>&1 | tee logs/02_mbdseq_by_bt2_$(date '+%Y-%m-%d-%H').log
	``` 

1. Generate bigWig files for ChIP-seq alignment files for visualization. If you have a ftp server to host bigWig files to import to UCSC genomebrowser, open the script and modify `ftp://account:passwd@hostname/apa_atingLab2019/chipseq_bigwigs`.
	```
	bash ./03_chipseq_tags_to_bigwig.sh
	bash ./04_input_tags_to_bigwig.sh
	```
	
1. Collect some basic alignment statistics
	```
	bash ./05_get_align_stats.sh
	```

1. Perform a differential ChIP-seq binding site analysis using MANorm.
	```
	bash ./06_manorm_diff_chipseq.sh
	Rscript ./06a_merge_manorm_report.r
	Rscript ./06b_store_manorm_sites.r
	```
1. Visualization with ChIP-seq differential binding sites. Note that b03_nmf.sh requires that MATLAB is installed and should be in $PATH.
	```
	Rscript ./b01_get_apa_intby_chipseq.r
	Rscript ./b02_uniq_xcvg_reg.r
	bash ./b03_nmf.sh
	Rscript ./b04_gen_heatmap_filtw.r
	Rscript ./b05_cluster_analysis.r
	```

# TCGA
TCGA data access permission is required to complete this analysis. Consider the following scripts to figure out which procedures/parameters were used in the paper for reference. Contact us for more help.

1. Change to TCGA working directory and download a preprocessed data (TCGA methylation beta values and pA usage ratio on APA genes of interest).
	```
	cd $HOME/projects/apa_atingLab2019/04_tcga
	bash ./01_resource_prep.sh
	less -S 546goi_all/github_Supplementary_Table_S6.tsv #table
	```

1. Access TCGA RNA-Seq BAM files and 4500 Infinium methylation array. Compute predicted polyA usage and normalized methylation level (beta value). Visualize the correlation with adjusted P-value in scatter plots and gene track along with ChIP-seq binding site read pileup from the cell line model designed.
	```
	./02a_tcga_analy_apa_paur.r
	./02b_tcga_analy_apa_paur_plot.sh
	```
1. To summarize/visualize the correlation between polyA usage predicted from mRNA-seq and methylation level,
	```
	cd ./03_tcga_summary
	Rscript ./01_collapse_corr_mpts_max_corr_per_cohort.r
	Rscrtip ./02_summary_corr_table.r
	```

1. To generate a box plot of polyA usage and methylation level in HEATR2 per tumor stage.
	```
	cd ../04_tcga_clinical_info
	bash ./runme.sh
	```

# Reference
