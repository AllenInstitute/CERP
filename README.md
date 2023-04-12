# cross-species Enhancer Ranking Pipeline (CERP)

## Introduction

![](https://github.com/AllenInstitute/peakCallingPipeline/blob/main/schematic.jpg)

## CERP outputs

The main CERP function `peakCaller` returns an ArchR object that contains the results of the peak calling and marker peak analysis.

### Output files

CERP creates the following files/directories:

1. `MarkerPeaks`: Directory under the ArchR project folder, contains marker peak .tsv for each annotation level.
2. `GroupBigWigs`: Directory under the ArchR project folder, contains bigwigs for each annotation level.

### ArchR project additions

CERP adds additional fields to: `archr.proj@projectMetadata` specifically:

1. `markerPeaks`: Stores the marker peaks for each subclass.
2. `loc_marker_table_"groupBy"`: Stores the location of the marker peak .tsv, where groupBy is the annotation level (class, subclass...) used for CERP.

## Setup pipeline

To setup the peak calling pipeline, start by cloning it in your server.

```
git clone https://github.com/AllenInstitute/CERP
```

Or install via a local file:
```
install.packages("/home/nelson.johansen/scripts/R/github/CERP_0.1.tar.gz")
```

Next, load the conda env that contains all the dependencies for this pipeline:
```
$ conda activate /allen/programs/celltypes/workgroups/hct/conda_envs/peak_caller
```

## Tutorials
[Basic CERP tutorial](https://github.com/AllenInstitute/peakCallingPipeline/blob/main/examples/CERP_v1_5.md)
