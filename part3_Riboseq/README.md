# Introduction to Ribo-seq data analysis and the Rp-Bp workflow


## Part 1.1: Virtual environments and installation of packages hosted on GitHub

### Sections:
    1.1.1 Short introduction to Python virtual environments: how does a virtual
    environments work, installing and using a virtual environment.
    1.1.2 Installation of packages hosted on GitHub.
    1.1.3 Short introduction to Jupyter (Notebook) and the JupyterHub

### Questions & Objectives:
    - Why do we need a virtual environment?
    - Learn how to install and use a Python virtual environment.
    - Learn how to install a package hosted on GitHub, in particular 
    how to install the rpbp package.
    - Jupyter Notebook basics: dashboard, user interface, navigation, running
    code, etc. How to use Jupyter with virtual environments.
    - How to use the JupyterHub.
### After I will be able to:
    - create and use a Python virtual environment to install packages.
    - understand what is a Jupyter Notebook, and how to use it.

    
## Part 1.2: Introduction to ribosome-profiling (Ribo-seq) and the Rp-Bp workflow

### Sections:
    1.2.1 Overview of the rpbp pipeline: command line options and configuration file
    1.2.2 Very short introduction to the Slurm workload manager (slurm-magic), used to run the rpbp pipeline

    1.2.3 High-level introduction to Ribo-seq, de novo ORF discovery (elements of annotation, transcript 
    isoforms, CDS, UTRs, etc.), biological relevance of alternative translation events (including 
    translation from non-coding transcripts), and why we need "dedicated software" to analyse Ribo-seq data.
    1.2.4 The rpbp pipeline step-by-step:
        - Creating reference genome indices;
        - Running the pipeline: creating ORF profiles, predicting translated ORFs.

### Questions & Objectives:
    - What is the translatome? What are the uses of Ribo-seq.
    - Why do we need dedicated software to analyse Ribo-seq data?
    - What softwares are available to analyse Ribo-seq data?
    - Understand how to use the rpbp package (on my laptop, on the cluster using the Slurm workload manager).
    - Run the complete rpbp pipeline on a selected Ribo-seq dataset.

### After I will be able to:
    - understand how to analyse Ribo-seq data for ORF discovery;
    - run the rpbp package (only ORF profiles, or full pipeline).
    
    
# Part 1.3: Ribo-seq quality control and downstream analysis of the Rp-Bp results

## Sections:
    - Ribo-seq quality control how-to.
    - Preprocessing analysis: standard rpbp preprocessing and more.
    - Prediction analysis: understand the output of rpbp.
    - Short introduction to visualisation using matplolib.

## Questions & Objectives:
    - Learn how to identify "good quality" Ribo-seq data.
    - Learn how to visualise the results.
    - Understand and use the output of rpbp (ORF predictions).

## After I will be able to:
    - run the rpbp downstream analysis pipeline and assess the quality
    of the data;
    - use the ORF predictions for follow-up studies.
    
    
# Part 2: Introduction to transcriptome-wide analysis of translational efficiency

## Sections:
    - Translational efficiency, translational control, and biological
    relevance.
    - Bioinformatics perspective: software and models for TE analysis 
    (high-level introduction).

## Questions & Objectives:
    - What is translational control, how RNA- and Ribo-seq together 
    overcome limitations of classical expression studies (RNA-seq alone).
    - Where to start if I have matched RNA- and Ribo-seq data and a 
    simple experimental design (treated vs. control)?

## After I will be able to:
    - have a better understanding of how to quantify translational
    control, and how to use and visualise the results.    
    
    
# Sofware, packages:
    - Python3
    - Jupyter
    - rpbp and dependencies (incl. slurm, bowtie2, etc.)
    - slurm-magic     
    
