library(arrow)
library(tidyverse)
source('/Users/jacobkang/VisiumHD/analysis/AuxFunctions.R')

SampleData <- GenerateSampleData('/Users/jacobkang/VisiumHD/PatientCRC1/outs', size='008um')
bcs <- SampleData$bcs
MetaData <- read_parquet('/Users/jacobkang/Outputs/MetaData/P1CRC_Metadata.parquet')
bcs <- left_join(bcs, MetaData %>% select(barcode, UnsupervisedL2), by='barcode')

cat('Non-NA UnsupervisedL2 in bcs:', sum(!is.na(bcs$UnsupervisedL2)), '\n')
cat('First few values:\n')
print(head(na.omit(bcs$UnsupervisedL2), 5))
cat('barcode match count:', sum(bcs$barcode %in% MetaData$barcode), '\n')
cat('Sample bcs barcodes:\n')
print(head(bcs$barcode, 3))
cat('Sample MetaData barcodes:\n')
print(head(MetaData$barcode, 3))