# RegBay: Bayesian network-based computational framework for regulatory element prediction.

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
## General workflow
- Train new prediction model
<br>    Peak calling/refinement (C) --> Model training (R)
```
Rscript regBay.R -C -i metadata.txt [options]
Rscript regBay.R -R -i metadata.txt [options]
```
- Predict regulatory elements (A pre-trained model is required.)
<br>    Peak calling/refinement (C) --> Predicting regulatory elements (P)
```
Rscript regBay.R -C -i metadata.txt [options]
Rscript regBay.R -P -i metadata.txt -m model.RDS [options]
```
## I/O files and formats
### Reuiqred input files
A comma delimited metadata file is required, which lists all input files.
- Peak calling/refinement module
    - BAM files: ATAC
    - (optional) ATAC-seq peaks
- Model training module
    - BAM files: ATAC, H3K27ac, H3K4me1, H3K4me3, input (empty or IgG, optional)
    - Annotation sets: active promoter, inactive promoter, active enhancer, inactive enhancer
    - (optional) ATAC-seq peaks
- Prediction module
    - BAM files: ATAC, H3K27ac, H3K4me1, H3K4me3, input (empty or IgG, optional)
    - (optional) ATAC-seq peaks

Example metadata.txt:
```
ATAC,sample.atac.bam
H3K27ac,sample.k27ac.bam
active_pmo,act_pmo.bed
```
### Output files
- Peak calling/refinement module
    - peaks.bed contains refined ATAC-seq peak regions.
- Training module
    - model.RDS contains a trained Bayesian network model.
- Prediction module
    - predict.promoter.bedGraph contains predicted promoters with predicted probabilities.
    - predict.enhancer.bedGraph contains predicted enhancers with predicted probabilities.
    - predictions.RDS contains predicted probabilities of regulatory elements for all candidate genomic regions.

## Parameters

### ATAC-seq peak calling / refinement
```
-l/--cutoff  <int>    Cutoff length for broad peaks (DEFAULT:1000)
```
Peaks wider than this value are classified as broad peaks and later refined into narrower peak units.
```
-w/--window  <int>    Window size of peak calling (DEFAULT:10000)
```
RegBay estimates the binomial negative distribution within each sliding genomic window. The dynamically estimated distribution is then used for peak calling.
```
-p/--pval  <float>    P-value for peak calling (DEFAULT:0.05)
```
Statistical significance threshold for peak calling based on binomial negative distribution.
### Training / prediction of RegBay model
```
-f/flank  <int>    Size of flanking regions (DEFAULT:500)
```
RegBay searches histone modification signals within the flanking regions. By default, this search range crosses approximate two histones.
## Contact

