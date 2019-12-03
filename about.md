###  Purpose

**AutobasisApp** is a **Shiny** web application that allows users to project their own GWAS summary statistics data onto a 13 Immune-mediated disease (IMD) basis.
The application uses a Bayesian shrinkage approach (contained in the `cupcake` package) to perform a disease-focused decomposition, allowing to summarise the complex patterns of IMD genetic risk onto a lower dimensional representation and evaluate IMD genetic architectures, identify and characterise disease-discriminating components, as well as suggest novel associations. This approach allows the characterisation of the genetic architecture of traits with modest sample size where traditional GWAS approaches are challenging. 

The application allows user interaction and creates interactive visualizations such as PCA and Forest (Delta) plots, as well as a downloadable table with the PCA results of the trait projected onto the basis.


### File requirements

In order to analyze your own data, you must provide a GWAS summary statistics file with the following features:

- TSV, CSV, TXT or compressed (i.e. GZ) formats.

- In **hg19/GRCh37** build.

- First line must include, at least, the following columns: 
  
    - **CHR** (Chromosome)
    - **POS** (Base position)
    - **REF** (Reference allele)
    - **ALT** (Alternative, or effect allele), 
    - **SE** (Standard Error of the log OR, BETA)
    - **BETA** (log OR) *or* **OR** (Odds Ratio)
    - **P** (P-value)
      
If OR is provided, BETA will be calculated automatically.
  
  
### Sample file

We included a filtered dataset that serves as a default dataset to showcase what should be expected when inputting your own data.
This dataset correspond to a GWAS of C-X-C motif chemokine 10, Interferon gamma-induced protein 10 (CXCL10, IP-10) levels, published by Aholla-Olli et al., 2017 ([10.1016/j.ajhg.2016.11.007](https://doi.org/10.1016/j.ajhg.2016.11.007)), and publicly available at [GWAS Catalogue](http://computationalmedicine.fi/data#Cytokine_GWAS).

**Note** that it is *essential* that your BETA is relative to the same alleles as in the original basis.  The example file gives the list of SNPs (defined by CHR and POS) and the REF and ALT alleles.

We will reject any files that do not meet this requirement with an error message.


### Dependencies

**AutobasisApp** has been developed using **R** and **Shiny** and is dependent on the following software and **R** packages:

|  |   |
--- | ----
**Software**   | 
R  | Language and environment for statistical computing and graphics
**R packages** |
annotSnpStats | Simple class to allow SnpMatrix objects to be manipulated together with SNP annotation files. Available at [Github](https://github.com/chr1swallace/annotSnpStats)
BiocManager | Provides tools for managing Bioconductor's package versioning and release system. Dependency of cupcake/annotSnpStats
cupcake | Provides a set of functions that makes GWAS summary statistics amenable to PCA. Available at [Github](https://github.com/ollyburren/cupcake)
data.table | Fast aggregation of large data, among other functions
DT | A Wrapper of the JavaScript Library 'DataTables'. It helps rendering the PCA table for visualization
dplyr  | A fast, consistent tool for working with data frame like objects, both in memory and out of memory
knitr  | Tool for dynamic report generation in R
R.utils | Utility functions useful when programming and developing R packages, dependency of cupcake and/or annotSnpStats packages
shiny  | Web Application Framework for R


### Citation

To cite **AutobasisApp** in publications, please use:

- Burren O (2020) "Characterisation of the genetic architecture of immune mediated disease through informed dimension reduction allows the separation of clinical subtypes" *Manuscript in preparation*


### Developer

**AutobasisApp** has been developed by

Guillermo Reales 

Department of Medicine  
University of Cambridge  
Cambridge, CB2 0AW, United Kingdom  
E-mail: gr440@cam.ac.uk  
URL: http://github.com/GRealesM  
Twitter: @GReales7  
