
########################
# Creating Heat maps   #
########################

# Import the library
library(gplots)

# First Heat Map: Contain (1) or not (0) TFBSs in TADs

  #Choose colors for the heat map
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 1000)

  #Import the data
table<-read.table("MatrixTFBS_YesNo.txt",sep="\t")
table<-as.matrix(table)
class(table)<-"numeric"

  #Get Heat map
heatmap.2(table,col=my_palette)


# Second Heat Map: Quantitative calcualtion of TFBSs in TADs
tableq<-read.table("MatrixTFBS_Quantitative.txt",sep="\t")
tableq<-as.matrix(tableq)
class(tableq)<-"numeric"

  #Choose colors for the heat map
my_palette <- colorRampPalette(c("green", "red", "blue"))(n = 1000)

  #Get Heat map
heatmap.2(tableq,col=my_palette)
