# Introduction to computational RNA biology: splicing, translation, and networks

## Part 2: Computational methods for alternative splicing identification

## Updated instructions for installing dependencies

```bash
module load R/4.1.1_deb10
module load hbigs_course_2022_tbb

cd ~/hbigs_course_2022/part2_Splicing
python3 -m venv splicing
source splicing/bin/activate

pip3 install --upgrade pip setuptools wheel
pip install snakemake 
pip install jupyter
```

Start a R-Console
```bash
R
# register the kernel
IRkernel::installspec(name='ir4.1.1', displayname='R4.1')
```

We also need to modify the kernel spec
```
jupyter kernelspec list
```

Use the R kernel, *e.g.* `~/.local/share/jupyter/kernels/ir4.1.1/kernel.json`

```
# add this line to the R kernel spec
"env": {"R_LIBS":"/biosw/hbigs_course_2022_tbb/1.0.0/rlib"},
```

This should look like:

```
{
  "argv": ["/biosw/R/4.1.1_deb10/lib/R/bin/R", "--slave", "-e", "IRkernel::main()", "--args", "{connection_file}"],
  "env": {"R_LIBS":"/biosw/hbigs_course_2022_tbb/1.0.0/rlib"},
  "display_name": "R4.1",
  "language": "R"
}
```


## A perspective to differential splicing

Here we aim to link the biological concepts to computational analysis.
In addition, we define concepts the elements of gene annotation.

1. The central dogma of molecular biology
1. RNA processing and splicing
1. The split-read

### Using workflows managers

This section aims to introduce the world of workflow managers, specifically  [**SnakeMake**](https://github.com/snakemake/snakemake). We apply these to differential transcript usage. However, workflow managers can be applied to many other topics. See [**ARMOR<**]<https://github.com/csoneson/ARMOR> for an example.

1. Why do we need a workflow manager?
1. What is workflow portability?
1. Basic elements of a SnakeMake workflow

### Workflows for detecting alternative splicing

Here we show how we developed a framework comprising several workflows to execute and integrate alternative splicing events.

1. Baltica, a framework to execute and integrate alternative splicing methods

### Workflow to discover new transcript isoforms

1. How StringTie discover new transcript isoforms from RNA-seq data
1. A deep dive on a SnakeMake workflow for StringTie
