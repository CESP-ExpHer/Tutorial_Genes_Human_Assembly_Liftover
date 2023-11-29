# Tutorial for Liftover Genes Positions (from hg19 to hg38)
*Created by:* Yazdan Asgari<br>
*Creation date:* 29 Nov 2023<br>
*Update:* 29 Nov 2023<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>

This repository contains an R script that performs the liftover of Genes positions data from hg19 to hg38 coordinates using the `rtracklayer`, `GenomicRanges`, and `vroom` packages. The liftover is achieved through the use of a chain file provided by UCSC.

## Prerequisites
Before running the script, make sure you have the required R packages installed:

```R
install.packages(c("rtracklayer", "GenomicRanges", "vroom"))
```

**IMPORTANT NOTE:** This script needs some column names to be available in your Genes data (*"chr"*, *"start"*, *"end"*, and *"gene_name"*). *"chr"* column is in **character format** (example: *"chr1"*). If their names are different in your input data, do not forget to change them in the script before using it.

## Data Processing Steps
1. Reading Genes data
```R
library(rtracklayer)
library(GenomicRanges)
library(vroom)

path_genes <- "/path/to/your/Genes/"
genes_hg19 <- vroom(file = paste0(path_genes, "Genes_hg19.txt"))
head(genes_hg19)

genes_hg19_sel <- genes_hg19
# IMPORTANT NOTE: We saw that sometimes if you did not add "strand" column info on your data, some start positions will add with +1 bp in the created hg38 file !!! So, we added this column to the data as well.
genes_hg19_sel["strand"] <- "+"

# Check data types of genes_hg19_sel
str(genes_hg19_sel)

# Use this section ONLY if you have "NA" data in "start" column
# Check for missing values in "start" column. If so, it then removes them
#sum(is.na(genes_hg19_sel$start))
#genes_hg19_sel <- genes_hg19_sel[!is.na(genes_hg19_sel$start),]
#sum(is.na(genes_hg19_sel$start))

```

2. Importing chain file
```R
# It was downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/
chain <- import.chain("hg19ToHg38.over.chain")
```
You could also download the "chain file" from the [0_files](/0_files) folder.
<br>
<br>
3. Creating GRanges object
```R
gr <- makeGRangesFromDataFrame(genes_hg19_sel, ignore.strand = F, seqnames.field = "chr", start.field = "start", end.field = "end")
```

4. Lifting over to hg38
```R
hg38 <- liftOver(gr, chain)
hg38_df <- as.data.frame(hg38)
```

5. Matching and updating data
```R
genes_hg19_sel$rownames_col <- as.numeric(rownames(genes_hg19_sel))
head(genes_hg19_sel)
match_pos <- match(genes_hg19_sel$rownames_col, hg38_df$group)

# Add the corresponding information from hg38_df to genes_hg19_sel
genes_hg19_sel$group[!is.na(match_pos)] <- hg38_df$group[match_pos[!is.na(match_pos)]]
genes_hg19_sel$seqnames[!is.na(match_pos)] <- hg38_df$seqnames[match_pos[!is.na(match_pos)]]
genes_hg19_sel$bp_hg38[!is.na(match_pos)] <- hg38_df$start[match_pos[!is.na(match_pos)]]
head(genes_hg19_sel)
```

6. Checking for missing values and writing output
```R
# Check for missing values in the "bp_hg38" column of "genes_hg19_sel" Genes file, then remove them
sum(is.na(genes_hg19_sel$start_hg38))
genes_hg19_sel <- genes_hg19_sel[!is.na(genes_hg19_sel$start_hg38),]
sum(is.na(genes_hg19_sel$start_hg38))

sum(is.na(genes_hg19_sel$end_hg38))
genes_hg19_sel <- genes_hg19_sel[!is.na(genes_hg19_sel$end_hg38),]
sum(is.na(genes_hg19_sel$end_hg38))

# Write the output to a file
vroom_write(genes_hg19_sel, file = "genes_hg38.txt", delim = "\t", quote = "none", col_names = TRUE)
```

**Notes**
- Ensure that the hg19 to hg38 chain file is present in the project directory.
- Verify that the required R packages are installed before running the script.

# Tutorial for Liftover Genes Positions (from hg38 to hg19)
Just simply use the exact previous procedure, but this time use *"hg38ToHg19.over.chain"* chain file (available in the [0_files](/0_files) folder) and switch the variable names *"hg19"* and *"hg38"* in the script.
