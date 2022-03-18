# HBIGS course 2021/22: Introduction to computational RNA biology
## Part 3 - Ribo-seq data analysis and the Rp-Bp workflow

## Description

This part of the lecture contains material for the processing and analysis of the *Ribo-seq* data from the example dataset.
The content of the lecture is described below in more detail with learning objectives. We will cover aspects of each section by means of
[Jupyter notebooks](https://jupyter.org/). These notebooks constitute the material for this part of the course.

## Dependencies

You don't need to install anything before the course. We will go through all the required steps of installing a Python 
virtual environment and required packages. We will *run* the notebooks in the [Dieterich Lab JupyterHub server](https://jupyter.dieterichlab.org), which
provides an access to our computational environments and resources.

    
## Course content

## Part 3.1: Virtual environments and installation of packages hosted on GitHub

### Sections:
    3.1.1 Short introduction to Python virtual environments: how does a virtual
    environments work, installing and using a virtual environment.
    3.1.2 Installation of packages hosted on GitHub.
    3.1.3 Short introduction to Jupyter (Notebook) and the JupyterHub

### Questions & Objectives:
    - Why do we need a virtual environment?
    - Learn how to install and use a Python virtual environment.
    - Learn how to install a package hosted on GitHub, in particular 
    how to install the Rp-Bp package.
    - Jupyter Notebook basics: dashboard, user interface, navigation, running
    code, etc. How to use Jupyter with virtual environments.
    - How to use the JupyterHub.
### After I will be able to:
    - Create and use a Python virtual environment to install packages.
    - Understand what is a Jupyter Notebook, and how to use it.

    
## Part 3.2: Introduction to ribosome-profiling (Ribo-seq) and the Rp-Bp workflow

### Sections:
    3.2.1 Overview of the Rp-Bp pipeline: command line options and configuration file
    3.2.2 Very short introduction to the Slurm workload manager (slurm-magic), used to run the Rp-Bp pipeline

    3.2.3 High-level introduction to Ribo-seq, de novo ORF discovery (elements of annotation, transcript 
    isoforms, CDS, UTRs, etc.), biological relevance of alternative translation events (including 
    translation from non-coding transcripts), and why we need "dedicated software" to analyse Ribo-seq data.
    3.2.4 The Rp-Bp pipeline step-by-step:
        - Creating reference genome indices.
        - Running the pipeline: creating ORF profiles, predicting translated ORFs.

### Questions & Objectives:
    - What is the translatome? What are the uses of Ribo-seq.
    - Why do we need dedicated software to analyse Ribo-seq data?
    - What softwares are available to analyse Ribo-seq data?
    - Understand how to use the Rp-Bp package (on my laptop, on the cluster using the Slurm workload manager).
    - Run the complete Rp-Bp pipeline on a selected Ribo-seq dataset.

### After I will be able to:
    - Understand how to analyse Ribo-seq data for ORF discovery.
    - Run the Rp-Bp package (only ORF profiles, or full pipeline).
    
    
## Part 3.3: Ribo-seq quality control and downstream analysis of the Rp-Bp results

## Sections:
    3.3.1 Ribo-seq quality control how-to.
    3.3.2 Preprocessing analysis: standard Rp-Bp preprocessing and more.
    3.3.3 Prediction analysis: understand the output of Rp-Bp, and short introduction to visualisation using matplolib.

## Questions & Objectives:
    - Learn how to identify "good quality" Ribo-seq data.
    - Learn how to visualise the results.
    - Understand and use the output of Rp-Bp (ORF predictions).

## After I will be able to:
    - Run the Rp-Bp downstream analysis pipeline and assess the quality
    of the data.
    - Use the ORF predictions for follow-up studies.
    
    
## Part 3.4: Introduction to transcriptome-wide analysis of translational efficiency

## Sections:
    3.4.1 Translational efficiency, translational control, and biological
    relevance.
    3.4.2 Bioinformatics perspective: software and models for TE analysis 
    (high-level introduction).

## Questions & Objectives:
    - What is translational control, how RNA- and Ribo-seq together 
    overcome limitations of classical expression studies (RNA-seq alone).
    - Where to start if I have matched RNA- and Ribo-seq data and a 
    simple experimental design (treated vs. control)?

## After I will be able to:
    - Have a better understanding of how to quantify translational
    control, and how to use and visualise the results.    

    
