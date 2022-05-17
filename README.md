# peakCallingPipeline

## Introduction

![](https://github.com/AllenInstitute/peakCallingPipeline/blob/main/schematic.jpg)

## Setup pipeline

To run the pipeline load the conda env that contains all the dependencies and the peakCallingPipeline:
```
$ conda activate /allen/programs/celltypes/workgroups/hct/conda_envs/atac
```

## Example: R
```
archr.proj = peakCaller(archr.proj = archr.proj, archr.genome="hg38", groupBy="subclass")
```

## Manual install of pipeline

To setup the peak calling pipeline, start by cloning it in your server.

```
git clone https://github.com/AllenInstitute/peakCallingPipeline
```

OR

```
install.packages("/allen/programs/celltypes/workgroups/hct/atac_pipeline/peakCallingPipeline_0.1.tar.gz")
```
