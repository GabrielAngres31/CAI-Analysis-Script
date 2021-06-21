
Dataset <- 
  read.table("C:/Users/gjang/Documents/GitHub/CAI-Analysis-Script/PARL0_04262021.csv",
   header=TRUE, stringsAsFactors=TRUE, sep=",", na.strings="NA", dec=".", 
  strip.white=TRUE)
editDataset(Dataset)
Dataset <- within(Dataset, {
  accession <- as.factor(accession)
})
Boxplot(fresh_weight~accession, data=Dataset, id=list(method="y"))
Dataset$D_div_W <- with(Dataset, diameter/ width)
Boxplot(D_div_W~accession, data=Dataset, id=list(method="y"))
Boxplot(D_div_W~accession, data=Dataset, id=list(method="y"))

