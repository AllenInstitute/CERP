# Tutorial: CERP (v1.5)

In this tutorial we demonstrate how to run CERP (v1.5) to produce annotated marker peak tables.

### Required inputs:

* Predefined ArchR project with annotations

### Recommended environment:

* Run CERP on a compute node from HPC to avoid h5 file issues.

The BICore hosted cerp docker is compatibile with CERP v1.5.

* (non-interactive) singularity exec docker://bicore/cerp:latest YOUR_SCRIPT.R
* (interactive) singularity shell docker://bicore/cerp:latest

IMPORTANT! Please do a local install of CERP within the bicore/cerp docker for now:

```{R}
install.packages("/home/nelson.johansen/scripts/R/github/CERP_0.1.tar.gz")
.libPaths()
library(CERP, lib.loc="/home/nelson.johansen/R/x86_64-pc-linux-gnu-library/4.2")
```

### Run CERP:
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
```

### Run CERP modules
#### Module: Gini index
Here we will show how to run CERP modules that extend the annotated marker peak table
```R
## (Module: Gini index) Compute specificity of each marker peak
archr.proj = peakSpecificity(archr.proj,
                             groupBy="supertype",
                             max_sample_size=500)
saveArchRProject(archr.proj)
```

#### Module: PeakRankR
Here we run PeakRankR to rank our marker peaks. First we load in the subclass annotated marker peaks from CERP:
```R
## Load in the annotated marker peak table from `peakAnnotation()`
marker.peak.table = read.table(file.path(getOutputDirectory(archr.proj), "MarkerPeaks/subclass/subclass_annotated_markerPeaks.tsv"), 
                                sep="\t", 
                                header=T)
archr.proj = peakRankeR_annotation(archr.proj, 
                                    groupBy="subclass", 
                                    marker.peak.table)
saveArchRProject(archr.proj)
```

### Outputs of CERP:

In the ArchR folder CERP produces:

* markerPeaks directory within the ArchR project folder with annotated peak tables.
