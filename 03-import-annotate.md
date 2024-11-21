---
title: "Importing and annotating quantified data into R"
source: Rmd
teaching: 80
output:
  html_document:
    df_print: paged
exercises: 40
---





::::::::::::::::::::::::::::::::::::::: objectives
-   Learn how to import the quantifications into a SummarizedExperiment object.
-   Learn how to add additional gene annotations to the object.
::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions
-   How can one import quantified gene expression data into an object suitable for downstream statistical analysis in R?
-   What types of gene identifiers are typically used, and how are mappings between them done?
::::::::::::::::::::::::::::::::::::::::::::::::::

## Load packages

In this episode we will use some functions from add-on R packages. In order to use them, we need to load them from our `library`:


``` r
suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(hgu95av2.db)
    library(SummarizedExperiment)
})
```

If you get any error messages about `there is no package called 'XXXX'` it means you have not installed the package/s yet for this version of R. See the bottom of the [Summary and Setup](https://carpentries-incubator.github.io/bioc-rnaseq/index.html) to install all the necessary packages for this workshop. If you have to install, remember to re-run the `library` commands above to load them. 

## Load data


In the last episode, we used R to download 4 files from the internet and saved them on our computer. But we do not have these files loaded into R yet so that we can work with them. The original experimental design in [Blackmore et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28696309/) was fairly complex: 8 week old male and female C57BL/6 mice were collected at Day 0 (before influenza infection), Day 4 and Day 8 after influenza infection. From each mouse, cerebellum and spinal cord tissues were taken for RNA-seq. There were originally 4 mice per 'Sex x Time x Tissue' group, but a few were lost along the way resulting in a total of 45 samples. For this workshop, we are going to simplify the analysis by only using the 22 cerebellum samples. Expression quantification was done using STAR to align to the mouse genome and then counting reads that map to genes. In addition to the counts per gene per sample, we also need information on which sample belongs to which Sex/Time point/Replicate. And for the genes, it is helpful to have extra information called annotation.
Let's read in the data files that we downloaded in the last episode and start to explore them:


### Counts


``` r
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", 
                   row.names = 1)
dim(counts)
```

``` output
[1] 41786    22
```

``` r
# View(counts)
```

Genes are in rows and samples are in columns, so we have counts for 41,786 genes and 22 samples. The `View()` command has been commented out for the website, but running it will open a tab in RStudio that lets us look at the data and even sort the table by a particular column. However, the viewer cannot change the data inside the `counts` object, so we can only look, not permanently sort nor edit the entries. When finished, close the viewer using the X in the tab. 
It looks like the rownames are gene symbols and the column names are the GEO sample IDs, which are not very informative for telling us which sample is what.

### Sample annotations

Next read in the sample annotations. Because samples are in columns in the count matrix, we will name the object `coldata`:


``` r
coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv",
                    row.names = 1)
dim(coldata)
```

``` output
[1] 22 10
```

``` r
# View(coldata)
```

Now samples are in rows with the GEO sample IDs as the rownames, and we have 10 columns of information. The columns that are the most useful for this workshop are `geo_accession` (GEO sample IDs again), `sex` and `time`.

It is often the case that we will have multiple objects with important information, and we will need to cross reference and compare the objects to get the full picture.

### Gene annotations
The counts only have gene symbols, which while short and somewhat recognizable to the human brain, are not always good absolute identifiers for exactly what gene was measured. For this we need additional gene annotations that were provided by the authors. 
We will import this data using a new function, read.delim(). Why are we using a new function, and not read.csv() like we did above? Some files have contents, like commas or quotation marks, that mean we need to use specialised functions. This is explained in the next paragraph - skip this if you are brand new to R and don't yet need the detail. 

The `count` and `coldata` files were in comma separated value (.csv) format, but we cannot use that for our gene annotation file because the descriptions can contain commas that would prevent a .csv file from being read in correctly. Instead the gene annotation file is in tab separated value (.tsv) format. 
Likewise, the descriptions can contain the single quote `'` (e.g., 5'), which by default R assumes indicates a character entry. So we have to use a more generic function `read.delim()` with extra arguments to specify that we have tab-separated data (`sep = "\t"`) with no quotes used (`quote = ""`). 
We also put in other arguments to specify that the first row contains our column names (`header = TRUE`), the gene symbols that should be our `row.names` are in the 5th column (`row.names = 5`), and that NCBI's species-specific gene ID (i.e., ENTREZID) should be read in as character data even though they look like numbers (`colClasses` argument). 
You can look up this details on available arguments by simply entering the function name starting with question mark. (e.g., `?read.delim`)


``` r
rowranges <- read.delim("data/GSE96870_rowranges.tsv", 
                        sep = "\t", 
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE, 
                        quote = "", 
                        row.names = 5)
dim(rowranges)
```

``` output
[1] 41786     7
```

``` r
# View(rowranges)
```

For each of the 41,786 genes, we have the `seqnames` (e.g., chromosome number), 
`start` and `end` positions, `strand`, `ENTREZID`, gene product description 
(`product`) and the feature type (`gbkey`). These gene-level metadata are 
useful for the downstream analysis. For example, from the `gbkey` column, we
can check what types of genes and how many of them are in our dataset:


``` r
table(rowranges$gbkey)
```

``` output

     C_region     D_segment          exon     J_segment      misc_RNA 
           20            23          4008            94          1988 
         mRNA         ncRNA precursor_RNA          rRNA          tRNA 
        21198         12285          1187            35           413 
    V_segment 
          535 
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Challenge: Discuss the following points with your neighbor

1. How are the 3 objects `counts`, `coldata` and `rowranges` related to each other in terms of their rows and columns?
2. If you only wanted to analyse the mRNA genes, what would you have to do keep just those (generally speaking, not exact codes)?
3. If you decided the first two samples were outliers, what would you have to do to remove those (generally speaking, not exact codes)?
  
::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution

1. In `counts`, the rows are genes just like the rows in `rowranges`. The columns in `counts` are the samples, but this corresponds to the rows in `coldata`. 
2. I would have to remember subset both the rows of `counts` and the rows of `rowranges` to just the mRNA genes.
3. I would have to remember to subset both the columns of `counts` but the rows of `coldata` to exclude the first two samples.

:::::::::::::::::::::::::::::::::::

You can see how keeping related information in separate objects could easily lead to mis-matches between our counts, gene annotations and sample annotations. 
This is why Bioconductor has created a specialized S4 class called a `SummarizedExperiment`. The details of a `SummarizedExperiment` object are covered extensively at the end of the [Introduction to data analysis with R and Bioconductor](https://carpentries-incubator.github.io/bioc-intro/60-next-steps.html#next-steps) workshop. 
As a reminder, let's take a look at the figure below representing the anatomy of the `SummarizedExperiment` class:

<img src="https://uclouvain-cbio.github.io/WSBIM1322/figs/SE.svg" width="80%" style="display: block; margin: auto;" />

It is designed to hold any type of quantitative 'omics data (`assays`) along with linked sample annotations (`colData`) and feature annotations with (`rowRanges`) or without (`rowData`) chromosome, start and stop positions. Once these three tables are (correctly!) linked, subsetting either samples and/or features will correctly subset the `assay`, `colData` and `rowRanges`. Additionally, most Bioconductor packages are built around the same core data infrastructure so they will recognize and be able to manipulate `SummarizedExperiment` objects. Two of the most popular RNA-seq statistical analysis packages have their own extended S4 classes similar to a `SummarizedExperiment` with the additional slots for statistical results: [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseqdataset)'s `DESeqDataSet` and [edgeR](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/DGEList-class)'s `DGEList`. No matter which one you end up using for statistical analysis, you can start by putting your data in a `SummarizedExperiment`. 

## Assemble SummarizedExperiment
We will create a `SummarizedExperiment` from these objects:

- The `count` object will be saved in `assays` slot    
- The `coldata` object with sample information will be stored in `colData` slot (_**sample metadata**_)    
- The `rowranges` object describing the genes will be stored in `rowRanges` slot (_**features metadata**_)     

Before we put them together, you ABSOLUTELY MUST MAKE SURE THE SAMPLES AND GENES ARE IN THE SAME ORDER! Even though we saw that `count` and `coldata` had the same number of samples and `count` and `rowranges` had the same number of genes, we never explicitly checked to see if they were in the same order. One quick way to check:



``` r
all.equal(colnames(counts), rownames(coldata)) # samples
```

``` output
[1] TRUE
```

``` r
all.equal(rownames(counts), rownames(rowranges)) # genes
```

``` output
[1] TRUE
```

``` r
# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in coldata like this (which is fine
# to run even if the first was TRUE):

tempindex <- match(colnames(counts), rownames(coldata))
coldata <- coldata[tempindex, ]

# Check again:
all.equal(colnames(counts), rownames(coldata)) 
```

``` output
[1] TRUE
```

:::::::::::::::::::::::::::::::::::::::  challenge

If the features (i.e., genes) in the assay (e.g., `counts`) and the gene
annotation table (e.g., `rowranges`) are different, how can we fix them? 
Write the codes. 

:::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::: solution


``` r
tempindex <- match(rownames(counts), rownames(rowranges))
rowranges <- rowranges[tempindex, ]

all.equal(rownames(counts), rownames(rowranges)) 
```


:::::::::::::::::::::::::::::::::::


Once we have verified that samples and genes are in the same order, we can 
then create our `SummarizedExperiment` object.


``` r
# One final check:
stopifnot(rownames(rowranges) == rownames(counts), # features
          rownames(coldata) == colnames(counts)) # samples

se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    rowRanges = as(rowranges, "GRanges"),
    colData = coldata
)
```




A brief recap of how to access the various data slots in a `SummarizedExperiment` and how to make some manipulations:


``` r
# Access the counts
head(assay(se))
```

``` output
             GSM2545336 GSM2545337 GSM2545338 GSM2545339 GSM2545340 GSM2545341
Xkr4               1891       2410       2159       1980       1977       1945
LOC105243853          0          0          1          4          0          0
LOC105242387        204        121        110        120        172        173
LOC105242467         12          5          5          5          2          6
Rp1                   2          2          0          3          2          1
Sox17               251        239        218        220        261        232
             GSM2545342 GSM2545343 GSM2545344 GSM2545345 GSM2545346 GSM2545347
Xkr4               1757       2235       1779       1528       1644       1585
LOC105243853          1          3          3          0          1          3
LOC105242387        177        130        131        160        180        176
LOC105242467          3          2          2          2          1          2
Rp1                   3          1          1          2          2          2
Sox17               179        296        233        271        205        230
             GSM2545348 GSM2545349 GSM2545350 GSM2545351 GSM2545352 GSM2545353
Xkr4               2275       1881       2584       1837       1890       1910
LOC105243853          1          0          0          1          1          0
LOC105242387        161        154        124        221        272        214
LOC105242467          2          4          7          1          3          1
Rp1                   3          6          5          3          5          1
Sox17               302        286        325        201        267        322
             GSM2545354 GSM2545362 GSM2545363 GSM2545380
Xkr4               1771       2315       1645       1723
LOC105243853          0          1          0          1
LOC105242387        124        189        223        251
LOC105242467          4          2          1          4
Rp1                   3          3          1          0
Sox17               273        197        310        246
```

``` r
dim(assay(se))
```

``` output
[1] 41786    22
```

``` r
# The above works now because we only have one assay, "counts"
# But if there were more than one assay, we would have to specify
# which one like so:

head(assay(se, "counts"))
```

``` output
             GSM2545336 GSM2545337 GSM2545338 GSM2545339 GSM2545340 GSM2545341
Xkr4               1891       2410       2159       1980       1977       1945
LOC105243853          0          0          1          4          0          0
LOC105242387        204        121        110        120        172        173
LOC105242467         12          5          5          5          2          6
Rp1                   2          2          0          3          2          1
Sox17               251        239        218        220        261        232
             GSM2545342 GSM2545343 GSM2545344 GSM2545345 GSM2545346 GSM2545347
Xkr4               1757       2235       1779       1528       1644       1585
LOC105243853          1          3          3          0          1          3
LOC105242387        177        130        131        160        180        176
LOC105242467          3          2          2          2          1          2
Rp1                   3          1          1          2          2          2
Sox17               179        296        233        271        205        230
             GSM2545348 GSM2545349 GSM2545350 GSM2545351 GSM2545352 GSM2545353
Xkr4               2275       1881       2584       1837       1890       1910
LOC105243853          1          0          0          1          1          0
LOC105242387        161        154        124        221        272        214
LOC105242467          2          4          7          1          3          1
Rp1                   3          6          5          3          5          1
Sox17               302        286        325        201        267        322
             GSM2545354 GSM2545362 GSM2545363 GSM2545380
Xkr4               1771       2315       1645       1723
LOC105243853          0          1          0          1
LOC105242387        124        189        223        251
LOC105242467          4          2          1          4
Rp1                   3          3          1          0
Sox17               273        197        310        246
```

``` r
# Access the sample annotations
colData(se)
```

``` output
DataFrame with 22 rows and 10 columns
                     title geo_accession     organism         age         sex
               <character>   <character>  <character> <character> <character>
GSM2545336 CNS_RNA-seq_10C    GSM2545336 Mus musculus     8 weeks      Female
GSM2545337 CNS_RNA-seq_11C    GSM2545337 Mus musculus     8 weeks      Female
GSM2545338 CNS_RNA-seq_12C    GSM2545338 Mus musculus     8 weeks      Female
GSM2545339 CNS_RNA-seq_13C    GSM2545339 Mus musculus     8 weeks      Female
GSM2545340 CNS_RNA-seq_14C    GSM2545340 Mus musculus     8 weeks        Male
...                    ...           ...          ...         ...         ...
GSM2545353  CNS_RNA-seq_3C    GSM2545353 Mus musculus     8 weeks      Female
GSM2545354  CNS_RNA-seq_4C    GSM2545354 Mus musculus     8 weeks        Male
GSM2545362  CNS_RNA-seq_5C    GSM2545362 Mus musculus     8 weeks      Female
GSM2545363  CNS_RNA-seq_6C    GSM2545363 Mus musculus     8 weeks        Male
GSM2545380  CNS_RNA-seq_9C    GSM2545380 Mus musculus     8 weeks      Female
             infection      strain        time      tissue     mouse
           <character> <character> <character> <character> <integer>
GSM2545336  InfluenzaA     C57BL/6        Day8  Cerebellum        14
GSM2545337 NonInfected     C57BL/6        Day0  Cerebellum         9
GSM2545338 NonInfected     C57BL/6        Day0  Cerebellum        10
GSM2545339  InfluenzaA     C57BL/6        Day4  Cerebellum        15
GSM2545340  InfluenzaA     C57BL/6        Day4  Cerebellum        18
...                ...         ...         ...         ...       ...
GSM2545353 NonInfected     C57BL/6        Day0  Cerebellum         4
GSM2545354 NonInfected     C57BL/6        Day0  Cerebellum         2
GSM2545362  InfluenzaA     C57BL/6        Day4  Cerebellum        20
GSM2545363  InfluenzaA     C57BL/6        Day4  Cerebellum        12
GSM2545380  InfluenzaA     C57BL/6        Day8  Cerebellum        19
```

``` r
dim(colData(se))
```

``` output
[1] 22 10
```

``` r
# Access the gene annotations
head(rowData(se))
```

``` output
DataFrame with 6 rows and 3 columns
                ENTREZID                product       gbkey
             <character>            <character> <character>
Xkr4              497097 X Kell blood group p..        mRNA
LOC105243853   105243853 uncharacterized LOC1..       ncRNA
LOC105242387   105242387 uncharacterized LOC1..       ncRNA
LOC105242467   105242467 lipoxygenase homolog..        mRNA
Rp1                19888 retinitis pigmentosa..        mRNA
Sox17              20671 SRY (sex determining..        mRNA
```

``` r
dim(rowData(se))
```

``` output
[1] 41786     3
```

``` r
# Make better sample IDs that show sex, time and mouse ID:

se$Label <- paste(se$sex, se$time, se$mouse, sep = "_")
se$Label
```

``` output
 [1] "Female_Day8_14" "Female_Day0_9"  "Female_Day0_10" "Female_Day4_15"
 [5] "Male_Day4_18"   "Male_Day8_6"    "Female_Day8_5"  "Male_Day0_11"  
 [9] "Female_Day4_22" "Male_Day4_13"   "Male_Day8_23"   "Male_Day8_24"  
[13] "Female_Day0_8"  "Male_Day0_7"    "Male_Day4_1"    "Female_Day8_16"
[17] "Female_Day4_21" "Female_Day0_4"  "Male_Day0_2"    "Female_Day4_20"
[21] "Male_Day4_12"   "Female_Day8_19"
```

``` r
colnames(se) <- se$Label
# Now when we look at the se object, we have informative sample IDs:
head(assay(se))
```

``` output
             Female_Day8_14 Female_Day0_9 Female_Day0_10 Female_Day4_15
Xkr4                   1891          2410           2159           1980
LOC105243853              0             0              1              4
LOC105242387            204           121            110            120
LOC105242467             12             5              5              5
Rp1                       2             2              0              3
Sox17                   251           239            218            220
             Male_Day4_18 Male_Day8_6 Female_Day8_5 Male_Day0_11 Female_Day4_22
Xkr4                 1977        1945          1757         2235           1779
LOC105243853            0           0             1            3              3
LOC105242387          172         173           177          130            131
LOC105242467            2           6             3            2              2
Rp1                     2           1             3            1              1
Sox17                 261         232           179          296            233
             Male_Day4_13 Male_Day8_23 Male_Day8_24 Female_Day0_8 Male_Day0_7
Xkr4                 1528         1644         1585          2275        1881
LOC105243853            0            1            3             1           0
LOC105242387          160          180          176           161         154
LOC105242467            2            1            2             2           4
Rp1                     2            2            2             3           6
Sox17                 271          205          230           302         286
             Male_Day4_1 Female_Day8_16 Female_Day4_21 Female_Day0_4
Xkr4                2584           1837           1890          1910
LOC105243853           0              1              1             0
LOC105242387         124            221            272           214
LOC105242467           7              1              3             1
Rp1                    5              3              5             1
Sox17                325            201            267           322
             Male_Day0_2 Female_Day4_20 Male_Day4_12 Female_Day8_19
Xkr4                1771           2315         1645           1723
LOC105243853           0              1            0              1
LOC105242387         124            189          223            251
LOC105242467           4              2            1              4
Rp1                    3              3            1              0
Sox17                273            197          310            246
```

``` r
# Now we will group our samples based on sex and time. 
# First, create a variable called group:
se$Group <- paste(se$sex, se$time, sep = "_")
se$Group
```

``` output
 [1] "Female_Day8" "Female_Day0" "Female_Day0" "Female_Day4" "Male_Day4"  
 [6] "Male_Day8"   "Female_Day8" "Male_Day0"   "Female_Day4" "Male_Day4"  
[11] "Male_Day8"   "Male_Day8"   "Female_Day0" "Male_Day0"   "Male_Day4"  
[16] "Female_Day8" "Female_Day4" "Female_Day0" "Male_Day0"   "Female_Day4"
[21] "Male_Day4"   "Female_Day8"
```

``` r
# Now we want to re-order all samples so they are grouped. To do this, we must 
# a) tell R to treat the data as a factor (a special type of data) 
# b) manually specify the order of the groups:

se$Group <- factor(se$Group, levels = c("Female_Day0","Male_Day0", 
                                        "Female_Day4","Male_Day4",
                                        "Female_Day8","Male_Day8"))
se <- se[, order(se$Group)]
colData(se)
```

``` output
DataFrame with 22 rows and 12 columns
                         title geo_accession     organism         age
                   <character>   <character>  <character> <character>
Female_Day0_9  CNS_RNA-seq_11C    GSM2545337 Mus musculus     8 weeks
Female_Day0_10 CNS_RNA-seq_12C    GSM2545338 Mus musculus     8 weeks
Female_Day0_8  CNS_RNA-seq_27C    GSM2545348 Mus musculus     8 weeks
Female_Day0_4   CNS_RNA-seq_3C    GSM2545353 Mus musculus     8 weeks
Male_Day0_11   CNS_RNA-seq_20C    GSM2545343 Mus musculus     8 weeks
...                        ...           ...          ...         ...
Female_Day8_16  CNS_RNA-seq_2C    GSM2545351 Mus musculus     8 weeks
Female_Day8_19  CNS_RNA-seq_9C    GSM2545380 Mus musculus     8 weeks
Male_Day8_6    CNS_RNA-seq_17C    GSM2545341 Mus musculus     8 weeks
Male_Day8_23   CNS_RNA-seq_25C    GSM2545346 Mus musculus     8 weeks
Male_Day8_24   CNS_RNA-seq_26C    GSM2545347 Mus musculus     8 weeks
                       sex   infection      strain        time      tissue
               <character> <character> <character> <character> <character>
Female_Day0_9       Female NonInfected     C57BL/6        Day0  Cerebellum
Female_Day0_10      Female NonInfected     C57BL/6        Day0  Cerebellum
Female_Day0_8       Female NonInfected     C57BL/6        Day0  Cerebellum
Female_Day0_4       Female NonInfected     C57BL/6        Day0  Cerebellum
Male_Day0_11          Male NonInfected     C57BL/6        Day0  Cerebellum
...                    ...         ...         ...         ...         ...
Female_Day8_16      Female  InfluenzaA     C57BL/6        Day8  Cerebellum
Female_Day8_19      Female  InfluenzaA     C57BL/6        Day8  Cerebellum
Male_Day8_6           Male  InfluenzaA     C57BL/6        Day8  Cerebellum
Male_Day8_23          Male  InfluenzaA     C57BL/6        Day8  Cerebellum
Male_Day8_24          Male  InfluenzaA     C57BL/6        Day8  Cerebellum
                   mouse          Label       Group
               <integer>    <character>    <factor>
Female_Day0_9          9  Female_Day0_9 Female_Day0
Female_Day0_10        10 Female_Day0_10 Female_Day0
Female_Day0_8          8  Female_Day0_8 Female_Day0
Female_Day0_4          4  Female_Day0_4 Female_Day0
Male_Day0_11          11   Male_Day0_11 Male_Day0  
...                  ...            ...         ...
Female_Day8_16        16 Female_Day8_16 Female_Day8
Female_Day8_19        19 Female_Day8_19 Female_Day8
Male_Day8_6            6    Male_Day8_6 Male_Day8  
Male_Day8_23          23   Male_Day8_23 Male_Day8  
Male_Day8_24          24   Male_Day8_24 Male_Day8  
```

``` r
# Finally, also factor the Label column to keep in order in plots:

se$Label <- factor(se$Label, levels = se$Label)
```




## Save SummarizedExperiment

This was a bit of code and time to create our `SummarizedExperiment` object. We will need to keep using it throughout the workshop, so it can be useful to save it as an actual single file on our computer to read it back in to R's memory if we have to shut down RStudio. To save an R-specific file we can use the `saveRDS()` function and later read it back into R using the `readRDS()` function. 


``` r
saveRDS(se, "data/GSE96870_se.rds")
rm(se) # remove the object!
se <- readRDS("data/GSE96870_se.rds")
```





## Session info


``` r
sessionInfo()
```

``` output
R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] hgu95av2.db_3.13.0          org.Hs.eg.db_3.19.1        
 [3] org.Mm.eg.db_3.19.1         AnnotationDbi_1.66.0       
 [5] SummarizedExperiment_1.34.0 Biobase_2.64.0             
 [7] MatrixGenerics_1.16.0       matrixStats_1.4.1          
 [9] GenomicRanges_1.56.2        GenomeInfoDb_1.40.1        
[11] IRanges_2.38.1              S4Vectors_0.42.1           
[13] BiocGenerics_0.50.0         knitr_1.48                 

loaded via a namespace (and not attached):
 [1] Matrix_1.7-1            bit_4.5.0               jsonlite_1.8.9         
 [4] highr_0.11              compiler_4.4.2          BiocManager_1.30.25    
 [7] renv_1.0.11             crayon_1.5.3            blob_1.2.4             
[10] Biostrings_2.72.1       png_0.1-8               fastmap_1.2.0          
[13] yaml_2.3.10             lattice_0.22-6          R6_2.5.1               
[16] XVector_0.44.0          S4Arrays_1.4.1          DelayedArray_0.30.1    
[19] GenomeInfoDbData_1.2.12 DBI_1.2.3               rlang_1.1.4            
[22] KEGGREST_1.44.1         cachem_1.1.0            xfun_0.49              
[25] bit64_4.5.2             memoise_2.0.1           SparseArray_1.4.8      
[28] RSQLite_2.3.7           cli_3.6.3               zlibbioc_1.50.0        
[31] grid_4.4.2              vctrs_0.6.5             evaluate_1.0.1         
[34] abind_1.4-8             httr_1.4.7              pkgconfig_2.0.3        
[37] tools_4.4.2             UCSC.utils_1.0.0       
```

::: keypoints
-   Depending on the gene expression quantification tool used, there are different ways (often distributed in Bioconductor packages) to read the output into a `SummarizedExperiment` or `DGEList` object for further processing in R.
-   Stable gene identifiers such as Ensembl or Entrez IDs should preferably be used as the main identifiers throughout an RNA-seq analysis, with gene symbols added for easier interpretation.
:::
