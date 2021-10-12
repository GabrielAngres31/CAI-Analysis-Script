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
library("multcompView")

library("HyperG")




value = "width"

setwd("C:\\Users\\gjang\\Documents\\GitHub\\CAI-Analysis-Script")

SWITCHBOARD.csvMAINFILE <- read.csv(file = 'PARL0_09092021.csv')

tukeyAnalyzer <- function(dataset, threshold, measure) {
   
   
   generate_label_df <- function(tukey_in, variable){
      # Extract labels and factor levels from Tukey post-hoc
      Tukey.levels <- tukey_in[[variable]][,4]
      Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
      
      #I need to put the labels in the same order as in the boxplot :
      Tukey.labels$accession=rownames(Tukey.labels)
      Tukey.labels=Tukey.labels[order(Tukey.labels$accession) , ]
      return(Tukey.labels)
   }
   
   model=lm( dataset[[measure]] ~ factor(dataset$accession) )
   ANOVA=aov(model)
   TUKEY <- TukeyHSD(x=ANOVA, 'factor(dataset$accession)', conf.level=0.95)
   LABELS <- generate_label_df(TUKEY, "factor(dataset$accession)")
   
   tukeyMeanMaker <- function(tukey_results) {
      png(
         paste0(
            SWITCHBOARD.DIRECTORY,
            "TukeyPlots\\Pairwise_Mean_Comparisons\\",
            measure,
            "--",
            threshold,
            "th.png"
         )
      )
      plot(tukey_results, las=1 , col="brown")
      dev.off()
   }
   
   tukeyHypergraphMaker <- function(tukey_labels) {
      
      groupElements <- c("a")
      for (grouptags in tukey_labels$Letters) {
         for (char in strsplit(grouptags, "")) {
            groupElements <- sort(unique(append(groupElements, char)))
         }
      }
      
      collectGroup <- function(tukeyResults, targetgroup) {
         basket <- c()
         for (accession in row.names(tukeyResults)) {
            if (grepl(targetgroup, tukeyResults[accession, 1])) {
               basket <- append(basket, accession)
            }
         }
         return(basket)
      }
      
      hypergroups <- as.list(rep(".", length(groupElements))) 
      names(hypergroups) <- as.vector(groupElements)
      print(hypergroups)
      for (group in groupElements) {
         hypergroups[[group]] <- collectGroup(tukey_labels, group)
      }
      print("Point 5")
      print(hypergroups)
      
      generated_hypergraph <- hypergraph_from_edgelist(hypergroups)
      png(
         paste0(
            SWITCHBOARD.DIRECTORY,
            "TukeyPlots\\Representative_Hypergraphs\\",
            measure,
            "--",
            threshold,
            "th.png"
         )
      )
      plot(generated_hypergraph)
      dev.off() 
   }
   
   tukeyMeanMaker(TUKEY)
   tukeyHypergraphMaker(LABELS)
   
}


tukeyAnalyzer(SWITCHBOARD.csvMAINFILE, "100", "width")




h <- hypergraph_from_edgelist(list(
    c(246, 319, 326, 580, 582, 854),
    c(319, 326, 585),
    c(325, 390, 839),
    c(326, 572, 584, 585),
    c(572, 584, 585, 845),
    c(242)
  )
)
plot(h)

# 
model=lm( SWITCHBOARD.csvMAINFILE[[value]] ~ factor(SWITCHBOARD.csvMAINFILE$accession) )
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'factor(SWITCHBOARD.csvMAINFILE$accession)', conf.level=0.95)
# plot(TUKEY , las=1 , col="brown")
# #png()
# #dev.off()
# LABELS <- generate_label_df(TUKEY , "factor(SWITCHBOARD.csvMAINFILE$accession)")





# ##############
# # I need to group the treatments that are not different each other together.

# #print(value)
# #print(LABELS)
# print(LABELS)
# 
# 



# 
# # A panel of colors to draw each group with the same color :
# my_colors <- c( 
#   rgb(143,199,74,maxColorValue = 255),
#   rgb(242,104,34,maxColorValue = 255), 
#   rgb(111,145,202,maxColorValue = 255)
# )
# 
# # Draw the basic boxplot
# a <- boxplot(SWITCHBOARD.csvMAINFILE[[value]] ~ SWITCHBOARD.csvMAINFILE$accession ,
#              ylim=c(min(SWITCHBOARD.csvMAINFILE$thickness) ,
#                     1.1*max(SWITCHBOARD.csvMAINFILE$thickness)) ,
#              col=my_colors[as.factor(LABELS[,1])] ,
#              ylab=value ,
#              main="")
# 
# # I want to write the letter over each box. Over is how high I want to write it.
# over <- 0.1*max( a$stats[nrow(a$stats),] )
# 
# #Add the labels
# text( c(1:nlevels(SWITCHBOARD.csvMAINFILE$accession)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )