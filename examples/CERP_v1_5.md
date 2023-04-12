# Tutorial: CERP (v1.5)

In this tutorial we demonstrate how to run CERP (v1.5) to produce annotated marker peak tables.

#### Required inputs:

* Predefined ArchR project with annotations

### Recommended environment:

* Run CERP on a compute node from HPC to avoid h5 file issues.

The BICore hosted cerp docker is compatibile with CERP v1.5.

* (non-interactive) singularity exec docker://bicore/cerp:latest YOUR_SCRIPT.R
* (interactive) singularity shell docker://bicore/cerp:latest

#### Run CERP:
```R
## Load the relevant libraries
library(ArchR)
library(CERP)

## Setup ArchR, this is required now.
addArchRGenome("mm10")
addArchRThreads(threads=10) 

## Load in ArchR projection
archr.proj = loadArchRProject(project.name)

## ArchR requires no "-" in annotations
archr.proj$supertype = gsub("-", "_", archr.proj$supertype)

## First we call marker peaks using standard ArchR methods
archr.proj = peakCaller(archr.proj = archr.proj,
                        groupBy="supertype", 
                        dataset="CERP_TUTORIAL")

## Now, annotate marker peaks using Allen Institute format
archr.proj = peakAnnotation(archr.proj,
                            groupBy = "supertype", 
                            dataset = "CERP_TUTORIAL",
                            publish = FALSE, ## Not implemented yet but will "push" new peak tables onto MolGen Shiny.
                            ucsc.user = "AIBSMolGen", 
                            ucsc.session = "RELEVANT_UCSC_BROWSER_SESSION")
saveArchRProject(archr.proj)

## (Module: Gini index) Compute specificity of each marker peak
archr.proj = peakSpecificity(archr.proj,
                             groupBy="supertype",
                             max_sample_size=500)
saveArchRProject(archr.proj)
```