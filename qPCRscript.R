#writing a script to analyze qPCR data

#install and load libraries and functions
install.packages("data.table")
install.packages("ggplot2")
install.packages("Hmisc")
library(data.table)
library(Hmisc)
library(ggplot2)


MakeDCQ <- function( MyData, AvgCq=AvgCq, Ref, GOI="GOI" ) {
  input <- MyData[ MyData$GeneID == GOI ,]
  DeltaCQ <- numeric( length=nrow( input ) )
  DeltaSD <- numeric( length=nrow( input ) )
  
  for( i in 1:nlevels( input$Sample )) {
    thissample <- levels( input$Sample )[i]
    thisGOI <- input$AvgCq[ input$Sample == thissample ]
    thisREF <- MyData$AvgCq[ MyData$GeneID == Ref & MyData$Sample == thissample ]
    GOISD <- input$SD[ input$Sample == thissample ]
    RefSD <- MyData$SD[ MyData$GeneID==Ref & MyData$Sample == thissample]
    DeltaCQ[i] <- thisGOI - thisREF
    DeltaSD[i] <- sqrt(GOISD^2+RefSD^2)
  }
  return( data.frame( cbind( input, DeltaCQ, DeltaSD ) ))
}

#upload data and check, #make sure to change file path for current analysis 
qPCR<-read.csv("C:/Users/Savvas Constantinou/Downloads/sjc 3-3-17 EO, scn4aa Bact, rps11.csv",stringsAsFactors=F)
head(qPCR)
unique(qPCR$Gene)
unique(qPCR$Sample)

#calculate average Cq and SD for each sample: make a dataframe to put it into
AvgCq = data.frame(Gene=character(0), GeneID=character(0), Sample=character(0), SampleID=character(0), AvgCq=numeric(0), SD=numeric(0) )
for(i in unique(qPCR$Gene)) {
  currentGene = qPCR[qPCR$Gene==i,]
  for(j in unique(currentGene$Sample)) {
    currentSample = currentGene[currentGene$Sample==j,]
    AvgCq = rbind(AvgCq, data.frame(Gene=i, GeneID=unique(currentSample$GeneID), Sample=j, SampleID=unique(currentSample$SampleID), AvgCq=mean(currentSample$Cq), SD=sd(currentSample$Cq)))
  }
}
print(AvgCq)
MyData=AvgCq

#calculate deltaCq values and delta SD
DCq1 <- MakeDCQ( AvgCq, AvgCq, "Ref1" )
DCq1
#can use below code if wanting the second Ref gene
#DCq1 <- MakeDCQ( AvgCq, AvgCq, "Ref2" )

#calculate DDCq. 
#specify control DCq value
control <- DCq1[DCq1$SampleID=="control", 7]
#Calculate DDCq = DCq - DCq of control
DDCq<-DCq1$DeltaCQ - control
#add DDCq to data frame with gene, gene ID, sample, Sample ID
DDCq<-cbind(DCq1, DDCq)
#calculate Relative Quantity and add to dataframe
RelQuant<- 2^(-DDCq$DDCq)
DDCq<-cbind(DDCq, RelQuant)
#Calculate Upper Error and add to dataframe
UpperError <- (2^-(DDCq$DDCq-DDCq$DeltaSD))-DDCq$RelQuant
LowerError <- DDCq$RelQuant-(2^-(DDCq$DDCq+DDCq$DeltaSD))
DDCq<-cbind(DDCq, UpperError, LowerError)

#Set limits of error bars.
limits <- aes(ymax=DDCq$RelQuant + DDCq$UpperError,ymin=DDCq$RelQuant - DDCq$LowerError)

#the below code restructures the graph to put in a different order. Remove for other datasets
DDCq$Sample <-factor(DDCq$Sample, levels=c("MOC","MOL","MOH"))

#plot the data
p <- ggplot(data=DDCq, aes(x=Sample, y=RelQuant, fill=Sample))
p + geom_col(position=position_dodge()) +
  #change the size and color of text to make more legable (changed for poster, defaults may be good for papers)
  theme(title=element_text(size=20), axis.text = element_text(size=18,color="black"), axis.title = element_text(size=18), 
        legend.title = element_text(size=14),legend.text = element_text(size=14)) +
  #add errorbars
  geom_errorbar( aes(ymax=RelQuant + UpperError,ymin=RelQuant - LowerError), width=0.2, position=position_dodge(0.9))+
  #add labels
  ggtitle(expression(paste("Relative Abundance of ", italic("scn4aa"), "(compared to B-actin) after MO treatment"))) +
  xlab("Fish") +
  ylab("Relative Expression") +
  #change title of legend
  scale_fill_hue(name="MO concentration 
  (mg MO/kg fish)", # Legend label, use darker colors
                 breaks=c("MOC", "MOL","MOH"),
                 labels=c("Control MO; 20.0", "scn4aa SB; 12.5 (low)","scn4aa SB; 20.0 (high)"))
