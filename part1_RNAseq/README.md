# Introduction to long and short RNA-seq analysis 

## Part 1: From raw sequencing reads to count tables and beyond

In this part of the course we fill focus on raw data processing. We will start with raw RNA seq reads. First, quality control is carried out. Afterwards, the reads will be pre-processed and cleaned up. Subsequently, the principal mapping step onto the reference genome will be carried out.  This part will close with read-to-gene assignment and differential expression analysis.

### Workflow in detail

1. First quality assessment with [**MultiQC**](https://multiqc.info/)
1. Removing sequencing adapters with [**FlexBar**](https://github.com/seqan/flexbar)
1. Removal of ribosomal RNA with  [**bowtie2**](https://github.com/BenLangmead/bowtie2)
1. Principal read mapping with  [**STAR**](https://github.com/alexdobin/STAR)
2. Read-to-gene assignment with [**subread/featurecount**](http://subread.sourceforge.net/)


### Data

The data for the course is available at `/pub/hbigs_course_2019/part_1/`.

