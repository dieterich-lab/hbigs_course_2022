# Copyright (C) 2019 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Auto-generated Makefile
# Date:				2019-10-21 13:56:40.811369
# Script Git revision:		12c1d48
# Generating script call:	['/home//tjakobi/repos/tjakobi/rnaseq-pipeline/create_makefile.py', '-l', '@new_names.csv', '-o', '@old_names.csv', '-g', '/biodb/', '-d', 'raw_data', '-g', '/beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star', '-G', '/beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf', '-R', '/beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor', '-r', '/beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38_96_repeatmasker.gtf', '-a', '/beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38_96.fa', '-u', 'course.Rdata', '-p']
# Generating machine:		t470-tjakobi


# 17.03.2022 Modified by Etienne Boileau for HBIGS 2022

################## MAKEFILE STARTS BELOW #######################


.PHONY: clean check_clean SYMLINK FASTQC FLEXBAR RRNA MAPPING PREPARE_SUBREAD SUBREAD

all: SYMLINK FASTQC FLEXBAR BOWTIE STAR SUBREAD MULTIQC
# Rules for renaming all files to user-identifiable names

SYMLINK: CREATE_RAW_DIR raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1.fastq.gz raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1.fastq.gz raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1.fastq.gz raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1.fastq.gz raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2.fastq.gz raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2.fastq.gz raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2.fastq.gz raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2.fastq.gz

CREATE_RAW_DIR: 
	mkdir -pv raw_data
	
raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1.fastq.gz: 
	 ln -s /pub/hbigs_course_2022/raw-data/downsampled/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1.fastq.gz raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1.fastq.gz

raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1.fastq.gz: 
	 ln -s /pub/hbigs_course_2022/raw-data/downsampled/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1.fastq.gz raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1.fastq.gz

raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1.fastq.gz: 
	 ln -s /pub/hbigs_course_2022/raw-data/downsampled/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1.fastq.gz raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1.fastq.gz

raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1.fastq.gz: 
	 ln -s /pub/hbigs_course_2022/raw-data/downsampled/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1.fastq.gz raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1.fastq.gz

raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2.fastq.gz: 
	 ln -s /pub/hbigs_course_2022/raw-data/downsampled/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2.fastq.gz raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2.fastq.gz

raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2.fastq.gz: 
	 ln -s /pub/hbigs_course_2022/raw-data/downsampled/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2.fastq.gz raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2.fastq.gz

raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2.fastq.gz: 
	 ln -s /pub/hbigs_course_2022/raw-data/downsampled/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2.fastq.gz raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2.fastq.gz

raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2.fastq.gz: 
	 ln -s /pub/hbigs_course_2022/raw-data/downsampled/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2.fastq.gz raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2.fastq.gz


# Rules for QC-ing all files with FastQC

FASTQC: CREATE_FASTQC_DIR fastqc/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1_fastqc.zip fastqc/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1_fastqc.zip fastqc/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1_fastqc.zip fastqc/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1_fastqc.zip fastqc/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2_fastqc.zip fastqc/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2_fastqc.zip fastqc/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2_fastqc.zip fastqc/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2_fastqc.zip

CREATE_FASTQC_DIR: 
	mkdir -pv fastqc

fastqc/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1_fastqc.zip: raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1.fastq.gz
	$(eval FASTQC_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS:= $(FASTQC_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS):$(shell sbatch --parsable slurm_fastqc.sh raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1.fastq.gz fastqc/ ))
	@echo "Creating FASTQC job call for SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS -> $(FASTQC_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS)"

fastqc/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1_fastqc.zip: raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1.fastq.gz
	$(eval FASTQC_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS:= $(FASTQC_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS):$(shell sbatch --parsable slurm_fastqc.sh raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1.fastq.gz fastqc/ ))
	@echo "Creating FASTQC job call for SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS -> $(FASTQC_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS)"

fastqc/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1_fastqc.zip: raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1.fastq.gz
	$(eval FASTQC_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF:= $(FASTQC_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF):$(shell sbatch --parsable slurm_fastqc.sh raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1.fastq.gz fastqc/ ))
	@echo "Creating FASTQC job call for SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF -> $(FASTQC_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF)"

fastqc/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1_fastqc.zip: raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1.fastq.gz
	$(eval FASTQC_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF:= $(FASTQC_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF):$(shell sbatch --parsable slurm_fastqc.sh raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1.fastq.gz fastqc/ ))
	@echo "Creating FASTQC job call for SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF -> $(FASTQC_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF)"

fastqc/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2_fastqc.zip: raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2.fastq.gz
	$(eval FASTQC_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS:= $(FASTQC_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS):$(shell sbatch --parsable slurm_fastqc.sh raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2.fastq.gz fastqc/ ))
	@echo "Creating FASTQC job call for SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS -> $(FASTQC_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS)"

fastqc/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2_fastqc.zip: raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2.fastq.gz
	$(eval FASTQC_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS:= $(FASTQC_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS):$(shell sbatch --parsable slurm_fastqc.sh raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2.fastq.gz fastqc/ ))
	@echo "Creating FASTQC job call for SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS -> $(FASTQC_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS)"

fastqc/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2_fastqc.zip: raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2.fastq.gz
	$(eval FASTQC_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF:= $(FASTQC_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF):$(shell sbatch --parsable slurm_fastqc.sh raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2.fastq.gz fastqc/ ))
	@echo "Creating FASTQC job call for SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF -> $(FASTQC_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF)"

fastqc/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2_fastqc.zip: raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2.fastq.gz
	$(eval FASTQC_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF:= $(FASTQC_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF):$(shell sbatch --parsable slurm_fastqc.sh raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2.fastq.gz fastqc/ ))
	@echo "Creating FASTQC job call for SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF -> $(FASTQC_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF)"


# Rules for removing adapter via flexbar

FLEXBAR: FASTQC CREATE_FLEXBAR_DIR flexbar/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz flexbar/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz flexbar/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz flexbar/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz
CREATE_FLEXBAR_DIR: 
	mkdir -pv flexbar
flexbar/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz flexbar/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_2.fastq.gz: 
ifneq ("$(wildcard fastqc/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1_fastqc.zip)","")
	$(eval FLEXBAR_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS:= $(FLEXBAR_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS):$(shell sbatch --parsable slurm_flexbar_paired_make.sh raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1.fastq.gz raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2.fastq.gz flexbar/ SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS))
else
	$(eval FLEXBAR_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS:= $(FLEXBAR_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS):$(shell sbatch --dependency=afterok$(FASTQC_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS) --parsable slurm_flexbar_paired_make.sh raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R1.fastq.gz raw_data/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_R2.fastq.gz flexbar/ SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS))
endif
	@echo "Creating Flexbar job call for SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS -> $(FLEXBAR_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS)"
flexbar/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz flexbar/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_2.fastq.gz: 
ifneq ("$(wildcard fastqc/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1_fastqc.zip)","")
	$(eval FLEXBAR_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS:= $(FLEXBAR_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS):$(shell sbatch --parsable slurm_flexbar_paired_make.sh raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1.fastq.gz raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2.fastq.gz flexbar/ SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS))
else
	$(eval FLEXBAR_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS:= $(FLEXBAR_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS):$(shell sbatch --dependency=afterok$(FASTQC_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS) --parsable slurm_flexbar_paired_make.sh raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R1.fastq.gz raw_data/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_R2.fastq.gz flexbar/ SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS))
endif
	@echo "Creating Flexbar job call for SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS -> $(FLEXBAR_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS)"
flexbar/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz flexbar/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_2.fastq.gz: 
ifneq ("$(wildcard fastqc/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1_fastqc.zip)","")
	$(eval FLEXBAR_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF:= $(FLEXBAR_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF):$(shell sbatch --parsable slurm_flexbar_paired_make.sh raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1.fastq.gz raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2.fastq.gz flexbar/ SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF))
else
	$(eval FLEXBAR_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF:= $(FLEXBAR_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF):$(shell sbatch --dependency=afterok$(FASTQC_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF) --parsable slurm_flexbar_paired_make.sh raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R1.fastq.gz raw_data/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_R2.fastq.gz flexbar/ SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF))
endif
	@echo "Creating Flexbar job call for SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF -> $(FLEXBAR_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF)"
flexbar/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz flexbar/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_2.fastq.gz: 
ifneq ("$(wildcard fastqc/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1_fastqc.zip)","")
	$(eval FLEXBAR_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF:= $(FLEXBAR_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF):$(shell sbatch --parsable slurm_flexbar_paired_make.sh raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1.fastq.gz raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2.fastq.gz flexbar/ SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF))
else
	$(eval FLEXBAR_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF:= $(FLEXBAR_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF):$(shell sbatch --dependency=afterok$(FASTQC_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF) --parsable slurm_flexbar_paired_make.sh raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R1.fastq.gz raw_data/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_R2.fastq.gz flexbar/ SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF))
endif
	@echo "Creating Flexbar job call for SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF -> $(FLEXBAR_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF)"


# Rules for Bowtie2

BOWTIE: FLEXBAR FASTQC CREATE_BOWTIE_DIR rrna/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz rrna/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz rrna/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz rrna/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz
CREATE_BOWTIE_DIR: 
	mkdir -pv rrna
rrna/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz rrna/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_2.fastq.gz: 
ifneq ("$(wildcard flexbar/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz.lengthdist)","")
	$(eval BOWTIE_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS:= $(BOWTIE_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS):$(shell sbatch --parsable slurm_bowtie2_rRNA_filter_paired_make.sh /beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor flexbar/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz flexbar/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_2.fastq.gz rrna/ SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS))
else
	$(eval BOWTIE_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS:= $(BOWTIE_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS):$(shell sbatch --dependency=afterok$(FLEXBAR_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS) --parsable slurm_bowtie2_rRNA_filter_paired_make.sh /beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor flexbar/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz flexbar/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_2.fastq.gz rrna/ SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS))
endif
	@echo "Creating Bowtie2 job call for SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS -> $(BOWTIE_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS)"
rrna/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz rrna/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_2.fastq.gz: 
ifneq ("$(wildcard flexbar/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz.lengthdist)","")
	$(eval BOWTIE_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS:= $(BOWTIE_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS):$(shell sbatch --parsable slurm_bowtie2_rRNA_filter_paired_make.sh /beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor flexbar/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz flexbar/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_2.fastq.gz rrna/ SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS))
else
	$(eval BOWTIE_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS:= $(BOWTIE_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS):$(shell sbatch --dependency=afterok$(FLEXBAR_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS) --parsable slurm_bowtie2_rRNA_filter_paired_make.sh /beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor flexbar/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz flexbar/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_2.fastq.gz rrna/ SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS))
endif
	@echo "Creating Bowtie2 job call for SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS -> $(BOWTIE_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS)"
rrna/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz rrna/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_2.fastq.gz: 
ifneq ("$(wildcard flexbar/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz.lengthdist)","")
	$(eval BOWTIE_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF:= $(BOWTIE_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF):$(shell sbatch --parsable slurm_bowtie2_rRNA_filter_paired_make.sh /beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor flexbar/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz flexbar/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_2.fastq.gz rrna/ SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF))
else
	$(eval BOWTIE_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF:= $(BOWTIE_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF):$(shell sbatch --dependency=afterok$(FLEXBAR_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF) --parsable slurm_bowtie2_rRNA_filter_paired_make.sh /beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor flexbar/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz flexbar/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_2.fastq.gz rrna/ SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF))
endif
	@echo "Creating Bowtie2 job call for SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF -> $(BOWTIE_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF)"
rrna/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz rrna/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_2.fastq.gz: 
ifneq ("$(wildcard flexbar/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz.lengthdist)","")
	$(eval BOWTIE_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF:= $(BOWTIE_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF):$(shell sbatch --parsable slurm_bowtie2_rRNA_filter_paired_make.sh /beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor flexbar/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz flexbar/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_2.fastq.gz rrna/ SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF))
else
	$(eval BOWTIE_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF:= $(BOWTIE_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF):$(shell sbatch --dependency=afterok$(FLEXBAR_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF) --parsable slurm_bowtie2_rRNA_filter_paired_make.sh /beegfs/biodb/genomes/homo_sapiens/rRNA_cluster/rRNA_precursor flexbar/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz flexbar/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_2.fastq.gz rrna/ SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF))
endif
	@echo "Creating Bowtie2 job call for SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF -> $(BOWTIE_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF)"


# Rules for STAR

STAR: BOWTIE FLEXBAR FASTQC CREATE_STAR_DIR star/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS/Aligned.noS.bam star/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS/Aligned.noS.bam star/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF/Aligned.noS.bam star/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF/Aligned.noS.bam
CREATE_STAR_DIR: 
	mkdir -pv star
star/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS/Aligned.noS.bam:
ifneq ("$(wildcard rrna/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz )","")
	$(eval STAR_COLLECT:= $(STAR_COLLECT):$(shell sbatch --parsable slurm_circtools_detect_paired_mapping_make.sh /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star/ rrna/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz rrna/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_2.fastq.gz star/ SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf))
else
	$(eval STAR_COLLECT:= $(STAR_COLLECT):$(shell sbatch --dependency=afterok$(BOWTIE_SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS) --parsable slurm_circtools_detect_paired_mapping_make.sh /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star/ rrna/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_1.fastq.gz rrna/SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS_2.fastq.gz star/ SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf))
endif
	@echo "Creating STAR job call for SRR7451179__Rep_2__MP60-2__mRNA_Seq__PBS -> $(STAR_COLLECT)"
star/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS/Aligned.noS.bam:
ifneq ("$(wildcard rrna/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz )","")
	$(eval STAR_COLLECT:= $(STAR_COLLECT):$(shell sbatch --parsable slurm_circtools_detect_paired_mapping_make.sh /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star/ rrna/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz rrna/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_2.fastq.gz star/ SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf))
else
	$(eval STAR_COLLECT:= $(STAR_COLLECT):$(shell sbatch --dependency=afterok$(BOWTIE_SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS) --parsable slurm_circtools_detect_paired_mapping_make.sh /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star/ rrna/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_1.fastq.gz rrna/SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS_2.fastq.gz star/ SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf))
endif
	@echo "Creating STAR job call for SRR7451187__Rep_1__MP60-1__mRNA_Seq__PBS -> $(STAR_COLLECT)"
star/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF/Aligned.noS.bam:
ifneq ("$(wildcard rrna/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz )","")
	$(eval STAR_COLLECT:= $(STAR_COLLECT):$(shell sbatch --parsable slurm_circtools_detect_paired_mapping_make.sh /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star/ rrna/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz rrna/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_2.fastq.gz star/ SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf))
else
	$(eval STAR_COLLECT:= $(STAR_COLLECT):$(shell sbatch --dependency=afterok$(BOWTIE_SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF) --parsable slurm_circtools_detect_paired_mapping_make.sh /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star/ rrna/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_1.fastq.gz rrna/SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF_2.fastq.gz star/ SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf))
endif
	@echo "Creating STAR job call for SRR7451182__Rep_1__ME60-1__mRNA_Seq__EGF -> $(STAR_COLLECT)"
star/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF/Aligned.noS.bam:
ifneq ("$(wildcard rrna/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz )","")
	$(eval STAR_COLLECT:= $(STAR_COLLECT):$(shell sbatch --parsable slurm_circtools_detect_paired_mapping_make.sh /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star/ rrna/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz rrna/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_2.fastq.gz star/ SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf))
else
	$(eval STAR_COLLECT:= $(STAR_COLLECT):$(shell sbatch --dependency=afterok$(BOWTIE_SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF) --parsable slurm_circtools_detect_paired_mapping_make.sh /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/star/ rrna/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_1.fastq.gz rrna/SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF_2.fastq.gz star/ SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf))
endif
	@echo "Creating STAR job call for SRR7451201__Rep_2__ME60-2__mRNA_Seq__EGF -> $(STAR_COLLECT)"


# Rules for MultiQC call

MULTIQC: BOWTIE FLEXBAR FASTQC STAR multiqc_report.html

multiqc_report.html: 
ifneq ("$(wildcard star/*/Log.final.out)","")
	$(eval MULTIQC:= $(shell sbatch --parsable slurm_multiqc_make.sh star/ flexbar/ rrna/ fastqc))
else
	$(eval MULTIQC:= $(shell sbatch --dependency=afterok$(STAR_COLLECT) --parsable slurm_multiqc_make.sh star/ flexbar/ rrna/ fastqc))
endif

	@echo "Creating main MultiQC job call: $(MULTIQC)"


# Rules for subread call

SUBREAD: PREPARE_SUBREAD subread/course.Rdata

subread/course.Rdata: 
ifneq ("$(wildcard subread/*.bam)","")
	$(eval SUBREAD:= $(shell sbatch --parsable slurm_subread.sh -b subread/ --gtf /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf --file subread/course.Rdata))
else
	$(eval SUBREAD:= $(shell sbatch --dependency=afterok:$(PRE_SUBREAD) --parsable slurm_subread.sh -b subread/ --gtf /beegfs/biodb/genomes/homo_sapiens/GRCh38_96/GRCh38.96.gtf  --file subread/course.Rdata))
endif
# 	sbatch --dependency=afterok:$(SUBREAD) slurm_pbzip2.sh subread/*.featureCounts
	@echo "Creating main SUBREAD job call -> $(SUBREAD)"

# Rules for prepare-subread call

PREPARE_SUBREAD: subread

subread: 
ifneq ("$(wildcard star/*/Log.final.out)","")
	$(eval PRE_SUBREAD:= $(shell sbatch --parsable slurm_prepare_subread.sh star/ subread/ ))
else
	$(eval PRE_SUBREAD:= $(shell sbatch --dependency=afterok$(STAR_COLLECT) --parsable slurm_prepare_subread.sh star/ subread/ ))
endif

	@echo "Creating PREPARE-SUBREAD job call -> $(PRE_SUBREAD)"


# Rules for clean up all files generated by the pipeline

check_clean:
	@echo "THE FOLLOWING DIRECTORIES WILL BE REMOVED RECURSIVELY:"
	@echo fastqc
	@echo flexbar
	@echo rrna
	@echo star
	@multiqc_data
	@echo slurm-*
	@( read -p "Are you sure? [y/N]: " sure && case "$$sure" in [yY]) true;; *) false;; esac )
clean: check_clean
	rm -rv fastqc
	rm -rv flexbar
	rm -rv rrna
	rm -rv star
	rm -rv multiqc_data
	rm -rv slurm-*


