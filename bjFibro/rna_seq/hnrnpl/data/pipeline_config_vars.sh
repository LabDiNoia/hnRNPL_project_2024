#!/bin/bash
set -e
set -u
set -o pipefail
export SCRIPT_CUR_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/scripts/
export DATA_INFO_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/data/
export DATA_INFO_FILE=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/data//data_info.txt
export REFR_GROUP_NAME=CTLi
export TEST_GROUP_NAME=hnRNPLi
export SE_or_PE=PE
export SEQ_STRANDED=No
export SPIKEIN_GENOME=none
export GENOME=hg38
export SJDB_OVHG_VAL=49
export GENOME_DIR=/lustre07/scratch/subrampg/binf_analyses_data//common/star_indexes/hg38/sjdbOH_49/
export GTF_FILE=/lustre07/scratch/subrampg/binf_analyses_data/common/hg38/gencode.v43.primary_assembly.annotation.gtf
export GFF3_FILE=/lustre07/scratch/subrampg/binf_analyses_data/common/hg38/gencode.v43.primary_assembly.annotation.gff3
export CHROM_SIZE_FILE=/lustre07/scratch/subrampg/binf_analyses_data/common/hg38/GRCh38.primary_assembly.genome.fa.fai
export RAW_FASTQ_DIR=/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/bjFibro/rna_seq/hnrnpl/data/raw_fastq/
export FILT_FASTQ_DIR=/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/bjFibro/rna_seq/hnrnpl/data/filteredFastq/
export FILT_FASTQ_SYML_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/data/filteredFastq/
export ALIGNED_BAM_DIR=/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/bjFibro/rna_seq/hnrnpl/data/alignedBAM/
export ALIGNED_BAM_SYML_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/data/alignedBAM/
export BW_DIR=/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/bw_files/
export BW_SYML_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/bw_files/
export BAM_PROG_FILE=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/scripts//slurm_outputs/BAM_progress_file.txt
export BAM_PROG_2P_FILE=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/scripts//slurm_outputs/BAM_progress_2Pass_file.txt
export NORM_FACTOR_FILE=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/scripts/../data/normFactor_file.txt
export FC_DSQ_DIR=/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/fcounts_deseq/
export FC_DSQ_SYML_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/fcounts_deseq/
export RMATS_DIR=/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/rMATS/
export RMATS_SYML_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/rMATS/
export BW_DIR=/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/bw_files/
export BW_SYML_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/bw_files/
export QC_DIR=/lustre07/scratch/subrampg/binf_analyses_data/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/qc_stats/
export QC_SYML_DIR=/home/subrampg/binf_analyses/hnrnpl_project/bjFibro/rna_seq/hnrnpl/results/qc_stats/
export SIRT1_RMATS_EVENT=No
export SIRT1_RMATS_EVENT=No
export SIRT1_RMATS_EVENT=Yes
