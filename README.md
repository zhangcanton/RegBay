# RegBay: Bayesian network based computational framework for regulatory element prediction.

### Usage message
```
Usage: Rscript regBay.R [module] -i metadata.txt [options]
Modules:
    -C, --peakcall  Peak calling/refining module
    -R, --train     Training module for RegBay regulatory element model
    -P, --predict   Prediction module for regulatory element using RegBay
I/O options:
    -i, --input     Path to input files (.txt format)
    -m, --model     Path to RegBay model (.rds format)
    -o, --outdir    Output dircetory (DEFAULT:current folder)
Analysis parameters:
    -p, --pval      P-value for peak calling (DEFAULT:0.05)
    -l, --cutoff    Cutoff length for broad peaks (DEFAULT:1000)
    -w, --window    Window size of peak calling (DEFAULT:10000)
    -f, --flank     Size of flanking regions (DEFAULT:500)
Other options:
    -t, --cores     Number of cores to use (DEFAULT:1)
    -h, --help      Show help message
    -v, --verbose   Show detailed message
```

## I/O files and formats
### Reuiqred input files
A comma delimited metadata file is required, which lists all input files.

- BAM files: ATAC, H3K27ac, H3K4me1, H3K4me3, input (empty or IgG, optional)
- Annotation sets: active promoter, inactive promoter, active enhancer, inactive enhancer
- (optional) ATAC-seq peaks

Example:

```
ATAC,sample.atac.bam
H3K27ac,sample.k27ac.bam
active_pmo,act_pmo.bed
```

### Output files
- Peak calling module: peaks.bed contains refined ATAC-seq peak regions.
- Training module: model.RDS contains a trained Bayesian network model.
- Predicting module: prediction.bed contains predicted probabilities for all regulatory element candidates.

## Parameters

### ATAC-seq peak calling / refinement
```
-l/cutoff
-p/pval
-w/window
```
### Training / prediction of RegBay model
```
-f/flank
```
## Contact and citation

