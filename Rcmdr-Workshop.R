boxplot(accession~fresh_weight,data=PARL, main="Box PLott",
   xlab="Accession", ylab="fresh_weight")


editDataset(PARL)

Boxplot(fresh_weight~accession, data=PARL, id=list(method="y"))
Boxplot(height~accession, data=PARL, id=list(method="y"))
Boxplot(width~accession, data=PARL, id=list(method="y"))
Boxplot(diameter~accession, data=PARL, id=list(method="y"))
Boxplot(thickness~accession, data=PARL, id=list(method="y"))


