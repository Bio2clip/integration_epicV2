# Data Preparation and Preprocessing Script

# This script reads samplesheet files and performs preprocessing steps and combining all 3 Illumina microarray types for methylation data analysis.

# -----------------------------------------
# Libraries
# -----------------------------------------
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(modelTsne)
library(ComplexHeatmap)
library(dplyr)

# -----------------------------------------
# Data Reading and Processing
# -----------------------------------------

# Read the samplesheet file
targets <- read.metharray.sheet("D:/IDAT/")

# Read 450K samples
rgSet <- read.metharray.exp("D:/IDAT/", targets = targets %>% filter(CHIP == "450K"), force = TRUE, verbose = TRUE)

# Read EPIC samples
rgSetEPIC <- read.metharray.exp("D:/IDAT/", targets = targets %>% filter(CHIP == "EPIC"), force = TRUE, verbose = TRUE)

# Merge both arrays into the 450K format
rgSet <- combineArrays(rgSetEPIC, rgSet, outType = "IlluminaHumanMethylation450k")

# Prepare metadata in the same order as the final dataframe merging all CHIP types
targets_EPIC <- targets %>% filter(CHIP == "EPIC")
targets_450K <- targets %>% filter(CHIP == "450K")
targets_EPICV2 <- targets %>% filter(CHIP == "EPICV2")
targets_FULL <- rbind(targets_EPIC, targets_450K, targets_EPICV2)

# Perform quality control
detP <- detectionP(rgSet)
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[, keep]
detP <- detP[, keep]

# Perform preprocessing (Noob method)
mSetSq <- preprocessNoob(rgSet)

# Filter based on detection p-values
detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep, ]

# Free up memory
mSetSq <- NULL
rgSet <- NULL

# Remove CHRXY regions using the IlluminaHumanMethylation450k annotation  
annEPIC <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% c("chrX", "chrY")])
mSetSqFlt <- mSetSqFlt[keep, ]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)

# Get betas and m values 
MVals <- getM(mSetSqFlt)
BValsC <- getBeta(mSetSqFlt)


# V2 Data Handling

# Purpose: This section handles the EPICv2 data, transforming it into a format compatible with the 450K platform.


# Let's handle the v2 data now.
# Currently, I need to store the EPICv2 samples in a specific directory because Sesame only accepts folders and does not consider the Samplesheet.

# Read the folder where the EPICv2 samples are stored using prep = QCDPB
# For more information, see: https://bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/sesame.html

# -----------------------------------------
# Data Preparation
# -----------------------------------------

betas <- openSesame("d:/idat/208181020018_R08C01")

# Transform EPICV2 into a 450K dataframe
betas <- mLiftOver(x = betas, target_platform = "HM450")
betas <- as.data.frame(betas)
#need to be fixed somehow :
colnames(betas) <-'208181020018_R08C01'

# Ensure both datasets have the same CpGs
betas <- as.data.frame(betas) %>%
  filter(row.names(betas) %in% rownames(BValsC))

# Combine both dataframes: the one with EPIC and 450K, and the EPICV2 transformed into 450K
betas <- cbind(BValsC, betas)


# -----------------------------------------
# Perfome T-sne
# -----------------------------------------

# Select the 25K most variable CpGs
meth <- selectVariables(data = na.omit(betas), method = "SD", no.variables = 25000, threads = 16)

# t-SNE wrapper
tsne <- modelTsne(meth, perplexity = 5, group = targets_FULL$Sample_GROUP, dim = 2)
plotTsne(model = tsne, method = "plotly")
