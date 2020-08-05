
# IMD Basis App

#### A Shiny App for projecting your own GWAS data onto the 13-Immune-mediated-trait basis

Updated: 05/08/2020

**Note:** This is a beta version. Some important features and
compatibility with some data formats and builds will be surely missing,
but will be implemented in future versions.

IMD Basis App can be accessed at: https://grealesm.shinyapps.io/IMDbasisApp/

##  Introduction

Understanding the aetiological relationships between multiple (>2) complex diseases can provide insights into shared pathogenic themes which can then be used to inform therapeutic intervention. However, this has proved challenging for a number of reasons. In our paper [**'Informed dimension reduction of clinically-related genome-wide association summary data characterises cross-trait axes of genetic risk'**](https://www.biorxiv.org/content/10.1101/2020.01.14.905869v3), we expand on these challenges and describe a novel approach *cupcake* that seeks to overcome them. In the paper we propose a statistical framework that allows the creation of a lower-dimensional basis that summarises the sharing of genetic risk across a set of clinically related diseases. Using publicly available GWAS summary stats we apply the framework to the study of immune-mediated disease (IMD) to create an IMD specific basis. This software allows a user to upload their own GWAS summary data and project onto this basis. We envisage four main scenarios where this might be useful:

1. You want to understand if a trait shares a common genetic risk component with IMD.
2. You suspect your trait or traits might have a genetic relationship/shared aetiology with a particular set of IMDs but you require robust evidence.
3. Your GWAS is clinically related to IMD but is of a modest sample size preventing elucidation of its genetic architecture.
4. You want to assess the genetic architecture of a trait in the general context of IMD, this can be useful, for example, in describing the differences between clinically related subtypes of disease.


## Methodology

A  detailed account of the methodology is available in the paper but here is a brief explanation aimed at broad audience. Assessing the shared genetic architecture for more than two traits is challenging, the number of variants assayed in a GWAS is often in the millions which quickly makes useful comparisons problematic. One approach is dimension reduction; here we make sure that effect sizes (odds ratios) for each disease are with respect to the same allele and then arrange them in a large matrix. We can then perform standard principal component analysis (PCA) on this matrix to create lower-dimensional summaries of the shared and distinct genetic architectures of the input diseases creating a object called a **basis**. Unfortunately this naive approach does not work, the relationships that we are interested in are obscured by other phenomena such as linkage disequillibrium, study heterogeneity (e.g. sample size)  and variant allele frequencies. To overcome this we adapt strategies from genetic fine mapping to train a Bayesian "lens" to focus our attention only on the disease-associated fraction of the genome. We can then apply the lens to matrix of disease variant effect sizes and use this perform PCA and successfully elucidate components of shared genetic risk. 

This lens can be applied to external GWAS datasets, and these reweighted effect sizes can be projected onto this **basis** and their location with basis-space observed. We developed a method to assess whether the location of a projected trait within basis-space is significantly different from what would be expected under the null, thus enabling the kinds of analyses detailed in the introduction.  


## How to use this tool

Given our focus on the IMD-associated fraction of the genome the basis only require a small subset of variants. The tool itself will attempt to filter uploaded files of GWAS summary statistics to include only these variants, however please note that files greater than 200Mb cannot be accommodated. To upload your own file click the **Upload** button in the *Upload your data* panel on the left hand size and follow the on screen dialogue. 

To improve accuracy, we have now included the possibility to choose the broad ancestry group of the population that was used in your GWAS (**Select population** menu), and use different LD matrices, tailored to the following populations: Admixed American, African, East Asian, European, and South Asian. We computed such matrices using data from 1000 Genomes Project Phase III. See [here](https://www.internationalgenome.org/faq/which-populations-are-part-your-study/) for more details about specific populations used to generate the matrices. Note that LD matrices are used to calculate delta P-values and confidence intervals, so matrix choice should not affect projection itself.

To label your trait on subsequent plots fill in the text box labelled **Trait name**.

Subject to  various checks that are carried out on the uploaded file, the panel **Uploaded data overview** will become populated.
  - **SNPs in the basis** reports the number of variants in the basis and is mainly a sanity check for us to make sure the you have access to the most up to date basis. 
  - **Uploaded SNPs overlapping basis** reports the number (and percentage) of variants that you have uploaded that overlap the basis. **Note that if this percentage is below 95% then results may be unreliable**. See the section on **missing data** for a further discussion of this.
  - **Overall P-value** reports whether over all components your trait is significantly different from baseline.
  - **Most significant component** reports on which component your trait appears to be the most significantly different from baseline.

After your data has been sucessfully projected you will have access to two further tabs on the right hand side. 

### Plot

This is a forest plot that shows for a principal component selected from the control panel on the left. The objective of this plot is for a user to be able to assess whether their traits location for a specific component significantly deviates from baseline and how these relate to traits used to create the basis. To create the baseline, we include a synthetic disease for which all variants have an effect size of zero. Using the example file for IP-10 cytokine concentration, we can see that for PC1 it lies to the left and therefore looks more similar to basis diseases of a more autoimmune rather than autoinflammtory nature. A button **Download plot** allows you to download the current plot as a PDF.

### Table

This table gives the full results for your trait across all components. 

- **PC** labels the principal component.
- **Var.Delta** measures the variance associated projection location with respect to the baseline of your trait onto this component.
- **Delta** shows the difference in location between the projected dataset and the baseline.
- **P** the significance of your projection location.
- **Trait** a label for your trait - see text box **Trait name**. 

For further details on these and their statistical derivations please see our paper. A button **Download table** allows you to download the current table as a csv file. The checkbox **Include basis traits in download** if selected means that **Delta** values for basis traits are also included. **Note** as these traits are used to create the basis they do not have a **Var.Delta** or **P** associated with them. This tabularised data can be useful for your own analyses and plots for example it can be used to create dendrograms of overall disease architectures. See this [vignette](https://chr1swallace.github.io/cupcake/articles/create-basis.html) for further details.  


###  GWAS summary file requirements

In order to analyze your own data, you must provide a GWAS summary statistics file with the following format:

- TSV, CSV, TXT or compressed (i.e. GZ) formats.

- Mapped to the **hg19/GRCh37** build.

- With the following colummns (same order is **not** a requirement), here **OR** is short for odds ratio and is relevant to case/control studies only. The basis actually uses the log of the OR or **BETA** so either can be supplied, with software performing the neccessary conversion. GWAS of quantitative traits can also be projected and by definition will use **BETA**.
  
    - **CHR** (Chromosome)
    - **POS** (Base position)
    - **REF** (Reference allele)
    - **ALT** (Alternative, or effect allele), 
    - **SE** (Standard Error of the log OR, BETA)
    - **BETA** (log OR) *or* **OR**
    - **P** (P-value)
      
**A side note on SNPs and file size**: This App caps the maximum uploaded file size at **200MB**, which might be too small for most current GWAS summary statistics datasets. Our method focuses on 566 SNPs for projection only, which means that if you filter your dataset to contain only those 566 SNPs you'll obtain the same results - and it will run faster! 
You can get the list of SNPs by downloading the example dataset, and a detailed explanation of how we chose those 566 SNPs in the publication (see Citation below).

### Aligning effect sizes 

It is important that the effect sizes (OR/BETA) for projected traits are with respect to the same allele. Currently we assume that the effect is with respect to the ALT allele. Thus prior to uploading please ensure that REF and ALT alleles for your study match the basis and that all effect sizes are with respect to the ALT allele. Further details on how to do this are available [here](https://chr1swallace.github.io/cupcake/articles/project-basis.html).

### Missing data

There is no requirement for a GWAS to be performed on a particular genotyping platform or to have undergone imputation and the variants used by the basis are those that are commonly available across the major platforms (in the paper we demonstrate this has no real effect on the resultant basis), however due to study design some variants may be missing. We provide an **uploaded file overview** panel that shows various metrics for the uploaded file. If you notice the **Uploaded SNPs overlapping basis** is below 95% resultant analysis may be compromised - see the paper for further details.      
  
## Example file

We included a filtered dataset that serves as a default dataset to showcase what should be expected when inputting your own data. You can download this example dataset [here](https://raw.githubusercontent.com/GRealesM/IMDbasisApp/master/data/Sample_dataset_B004_Ahola-Olli_27989323_1.tsv).
This dataset correspond to a GWAS of C-X-C motif chemokine 10, Interferon gamma-induced protein 10 (CXCL10, IP-10) levels, published by Aholla-Olli et al., 2017 ([10.1016/j.ajhg.2016.11.007](https://doi.org/10.1016/j.ajhg.2016.11.007)), and publicly available at [GWAS Catalog] (https://www.ebi.ac.uk/gwas/efotraits/EFO_0008056) or [here](http://computationalmedicine.fi/data#Cytokine_GWAS).

This is an example of dataset format that *should*
work:

| SNP\_ID    | CHR |       POS | REF | ALT |     BETA |     SE |      P |
| :--------- | --: | --------: | :-- | :-- | -------: | -----: | -----: |
| rs11190140 |  10 | 101291593 | T   | C   | \-0.0016 | 0.0233 | 0.9422 |
| rs11195128 |  10 | 112186148 | C   | T   | \-0.0049 | 0.0247 | 0.8415 |
| rs3814231  |  10 | 115481018 | C   | T   |   0.0007 | 0.0261 | 0.9782 |
| rs3793910  |  10 | 131562993 | G   | T   | \-0.0100 | 0.0421 | 0.7975 |
| rs10829131 |  10 |  27177245 | A   | G   |   0.0284 | 0.0350 | 0.4136 |
| rs11008080 |  10 |  30802799 | A   | C   |   0.0069 | 0.0246 | 0.7694 |

## Dependencies

**IMD basis App** has been developed using **R** and **Shiny** and is dependent on the following software and **R** packages:

|  |   |
--- | ----
**Software**   | 
R  | Language and environment for statistical computing and graphics
**R packages** |
BiocManager | Provides tools for managing Bioconductor's package versioning and release system. Dependency of cupcake/snpStats.
cupcake | Provides a set of functions that makes GWAS summary statistics amenable to PCA. Available at [Github](https://github.com/ollyburren/cupcake)
data.table | Fast aggregation of large data, among other functions
DT | A Wrapper of the JavaScript Library 'DataTables'. It helps rendering the PCA table for visualization
dplyr  | A fast, consistent tool for working with data frame like objects, both in memory and out of memory
knitr  | Tool for dynamic report generation in R
R.utils | Utility functions useful when programming and developing R packages, dependency of cupcake and/or annotSnpStats packages
shiny  | Web Application Framework for R

## Data treatment and privacy statement

IMD basis App can use user-supplied data, but it does not store or shares them with any third party, and it deletes all data after being closed by the user. It will do so too if the server automatically closes the connection due to idle time (usually 5 minutes). This is the default behaviour of Shiny Apps hosted in ShinyApps.io (See [here](https://docs.rstudio.com/shinyapps.io/Storage.html)). Users can also check our [source code](https://github.com/GRealesM/IMDbasisApp).
In case of further privacy concerns, we refer the user to perform projections locally using our [cupcake package](https://github.com/ollyburren/cupcake), which applies the same method as the IMD basis App.

## Citation

To cite this application, please use:

- Burren OS *et al.* (2020) "Informed dimension reduction of clinically-related genome-wide association summary data characterises cross-trait axes of genetic risk" *bioRxiv*. doi:[10.1101/2020.01.14.905869](https://www.biorxiv.org/content/10.1101/2020.01.14.905869v3).


## About

This software has been developed by Guillermo Reales (gr440 [at] cam [dot] ac [dot] uk) and Olly Burren (ob219 [at] cam [dot] ac [dot] uk) within the  [**Wallace
Group**](https://chr1swallace.github.io) and funded by the **Wellcome
Trust**.