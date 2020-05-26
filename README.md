# bsaR

A simple R packages for bulk segregation analysis.

## Installation

```r
install.packages("vcfR")
remotes::install_github("xuzhougeng/bsaR")
```

## Usage

The input file should be processed by the GATK, e.g "gatk-hc.final.snp.vcf.gz"

Step1: Create the BSA project from the vcf file

```r
library(vcfR)
library(bsaR)

# Calculate the allele frequency ------------------------------------------
vcf.file <- "gatk-hc.final.snp.vcf.gz"

bsa <- CreateBsaFromVcf(vcf.file = vcf.file)
```

Step2: Filter the BSA with the following standard

- only keep bi
- exclude the NA row
- filter the only depth loci

```r
bsa <- FilterMultiVariant(bsa)
bsa <- FilterNaVariant(bsa)

bsa <- CalcDepth(bsa)
bsa <- FilterLowDepth(bsa, depth = 40)
```

Step3:  calculte the frequency

```r
bsa <- CalcAltFreq(bsa)
```

Step3.5: filter the potential error loci with low frequency


Step4: Visualization

"sample1" and "sample2" should be  sample name in the VCF file

```r
# Visualization -----------------------------------------------------------

bsa <- CalcFreqByWindow(bsa, window.size = 5000)
# plot per contig
# plot for contig  ------------------------------------------------------

pdf("dotplot.pdf")
for (i in paste0("chr", 1:8)){
  plotWindowFreq(bsa, i, "sample1", "sample2",
                 window.size = 5000)
}
dev.off()
```
