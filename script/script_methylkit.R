# ---------------------------
##
## Script name: script-methylkit-pipe.R
##
## Purpose of script:
##
## Author: Abdeljalil
##
## Date Created: 2020-02
##
## V 1.0.0
## Email: abdeljalil.senhajirachik@cnrgh.fr
##
## ---------------------------

#======================#
##### Load library #####
#======================#

#library(methylKit)

#=========================#




##########################################################################################

#====================#
##### logs_file: #####
#====================#

#Diverts R output to a connection:
log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

#=================#
##### Library #####
#=================#

library(methylKit)


#===================#
##### Functions #####
#===================#

# Extract and merge dataframe from multiple methread object
#@ methReadObj      : a methylRawList object
mergeMethylkitOutput <- function(methReadObj){
  
 # Loop over each samples processed by methread
  for(numList in 1:length(methReadObj)){
    
    # extract the values
    temp <- getData(methReadObj[[numList]])
    
    # rename specific columns (coverage, numCs and numTs) by addind samples id
    colnames(temp)[5:7] <- paste0(getSampleID(myobj[[numList]]),".",colnames(temp)[5:7])
    
    # If first one initiate the merged tabled else merge with the global dataframe
    if(numList == 1){
      merged <- temp
    } else {
      merged <- merge(merged,temp,all = T,by = c("chr","start","end","strand"))
    }
    
  }
 
  # return the dataframe
  return(merged)
   
}

#==========================#
##### Global variables #####
#==========================#

# Get working context from snakemake rule parameters
context <- snakemake@params[["context"]]
print(context)

# Get context directory from snakemake rule parameters
directory <- snakemake@params[["path"]]
print(directory)
# Get the context files of each samples
files <- list.files(directory)

# Named each samples using file name
sample_id <- gsub("\\..*","",files)

# Turn samples names and files as list for methylkit input
files <- as.list(paste0(directory,files))
sample_id <- as.list(sample_id)

table <- snakemake@output[["out"]]

# print(files)
# print(sample_id)

#============================#
##### Methylkit analysis #####
#============================#

# Read methylation using methylkit function methRead
# Here argument pipeline is used to precise the column of the input file since methylratio.py output is not natively recognize
# All samples are set to group 0 (treatment) since we don't intend to do differential analyses here
# minimum of coverage is set to 1 to include all positions (can be filtered out after using the final file)
myobj <- methRead(location = files, sample.id = sample_id, assembly = "populus tricharpa v3.1", mincov = 1, context = context,treatment = rep(0, length(files)),
                pipeline = list(fraction=TRUE, chr.col=1, start.col=2, end.col=2, coverage.col=6, strand.col=3, freqC.col=5 ))  

# Concatenate all samples tables into one unique table
finalFrame <- mergeMethylkitOutput(myobj)

#Write the final table as a csv2 file
write.csv2(finalFrame,file = table,)

# head(myobj)

# plots for statistcs and coverage simple :

pdf(file = snakemake@output[["plot"]])
getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
dev.off()




