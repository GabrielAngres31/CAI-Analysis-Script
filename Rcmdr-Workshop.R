
load("C:/Users/gjang/Documents/GitHub/CAI-Analysis-Script/PARL0_09092021.csv")
editDataset(SWITCHBOARD.csvMAINFILE)
SWITCHBOARD.csvMAINFILE <- within(SWITCHBOARD.csvMAINFILE, {
  accession <- as.factor(accession)
})
Boxplot(diameter~accession, data=SWITCHBOARD.csvMAINFILE, 
  id=list(method="y"))
Boxplot(fresh_weight~accession, data=SWITCHBOARD.csvMAINFILE, 
  id=list(method="y"))
Boxplot(height~accession, data=SWITCHBOARD.csvMAINFILE, id=list(method="y"))
Boxplot(thickness~accession, data=SWITCHBOARD.csvMAINFILE, 
  id=list(method="y"))
Boxplot(width~accession, data=SWITCHBOARD.csvMAINFILE, id=list(method="y"))
SWITCHBOARD.csvMAINFILE$D_div_W <- with(SWITCHBOARD.csvMAINFILE, 
  diameter/width)
Boxplot(D_div_W~accession, data=SWITCHBOARD.csvMAINFILE, 
  id=list(method="y"))

