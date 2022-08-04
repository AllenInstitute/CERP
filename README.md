# peakCallingPipeline

## Introduction

![](https://github.com/AllenInstitute/peakCallingPipeline/blob/main/schematic.jpg)

## Setup pipeline

To setup the peak calling pipeline, start by cloning it in your server.

```
git clone https://github.com/AllenInstitute/peakCallingPipeline
```

Next, load the conda env that contains all the dependencies for this pipeline:
```
$ conda activate peakCallingPipeline
```

## Example: R
```
archr.proj = peakCaller(archr.proj = archr.proj, archr.genome="hg38", groupBy="subclass", dataset="region-information")
```
