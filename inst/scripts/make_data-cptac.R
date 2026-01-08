

library("MsDataHub")
peptideFile <- cptac_peptides.txt()
peptides <- read.delim(peptideFile)

## Subset the data to make a small dataset
allProteins <- unique(peptides$Proteins)
upsProteins <- grep("ups", allProteins, value = TRUE)
yeastProteins <- grep("ups", allProteins, value = TRUE, invert = TRUE)
set.seed(1234)
keepProteins <- c(sample(upsProteins, 10), sample(yeastProteins, 50))
peptides <- peptides[peptides$Proteins %in% keepProteins, ]

## Create sample annotations
quantCols <- grep("Intensity[.]", names(peptides), value = TRUE)
coldata <- data.frame(quantCols = quantCols)
coldata$lab <- rep(rep(paste0("lab", 1:3), each = 3), 5)
coldata$condition <- gsub("Intensity..(.)_.*", "\\1", quantCols)
concentrations <- c(A = 0.25, B = 0.74, C = 2.22, D = 6.67, E = 20)
coldata$concentration <- concentrations[coldata$condition]

## Remove unnecessary feature annotations
keepAnnot <- c(
    "Sequence", "Proteins", "Charges", "PEP", "Reverse",
    "Potential.contaminant", quantCols
)
peptides <- peptides[, keepAnnot]

## Create the QFeatures object
cptac <- readQFeatures(
    peptides, coldata, name = "peptides", fnames = "Sequence"
)

## Store data
save(cptac, file = file.path("data/cptac.rda"),
     compress = "xz", compression_level = 9)
