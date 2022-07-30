library("xlsx")
library("dplyr")
library("tibble")
library("rsq")
library("qwraps2")
library("ggplot2")
library("car")
library("dunn.test")
library("Rcmdr")
library("hash")

setwd("C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script\\SAS_Datasets")
df <- read.csv(file = 'PARL0C_95th.csv')
df95 <- read.csv(file = 'PARL0C_95th.csv')
df90 <- read.csv(file = 'PARL0C_90th.csv')

padeApprox <- function(height, width) {
  h <- ((height-width)/(height+width))^2
  return (pi*(height+width)*(64-3*h^2)/(64-16*h))
}

df$padeEllip <- padeApprox(df$height, df$width)
df$theoEllip <- pi * (df$height/2) * (df$width/2)

print(
  summary(
    lm(
      Area ~ padeEllip, df[which(df$Area != "N/A"),]
    )
  )
)

plot(Area ~ padeEllip, df[which(df$Area != "N/A"),])
text(Area ~ padeEllip, labels=accession, data = df[which(df$Area != "N/A"),], cex=0.7, font=2)


