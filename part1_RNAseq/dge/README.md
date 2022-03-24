# HBIGS course 2021/22: Introduction to computational RNA biology
## Part 1 - From raw sequencing reads to differential expression

## Description

This part of the lecture contains material for the Differential Gene Expression (DGE) analysis of the Gene Expression data generated from first part of the lecture.

## Prerequisites
For running the analysis, the participiants need to do the following:

* Go to [https://r.dieterichlab.org:49200/](https://r.dieterichlab.org:49200/)
* Add the Log-in credentials (i.e. **Username:** coursexx and **Password:** coursexx)
* Set the working directory to: ```setwd("/beegfs/pub/hbigs_course_2022/part1_RNAseq/dge/")```
* Add the libPath: ```.libPaths <- .libPaths( c( .libPaths(), "/biosw/....") )```
* Check if the added libPath is present in the list: ```.libPaths()```
* Open *dge_analysis.R* script.

From here on we continue with the analysis workflow.