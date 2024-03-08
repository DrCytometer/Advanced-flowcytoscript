# get channels script
# run this to get the channel names for the flowcytoscript parameter file

library(dplyr)
library( RColorBrewer )

# source parameters

param.filename <- "./flowcytoscript_parameter.r"
source( param.filename )

flow.csv <- read.csv(list.files(data.dir, "\\.csv$", full.names = TRUE)[1])

channels <- data.frame(colnames(flow.csv))
colnames(channels)[1] <- "name" 

for(i in 1:dim(channels)[1] ){
  channels$out1[i] = paste0("\"", channels$name[i], "\", #", channels$name[i], "\n")
  channels$out2[i] = paste0("\"", channels$name[i], "\" = \"", channels$name[i], "\",\n") 
}

cat('"', paste0(channels$name, collapse = '","'), '"', sep = "")

# copy the previous output into the channels.of.interest and delete unnecessary markers

channels.of.interest = c("Ki67","CCR9",
                         "CCR7","Ly.6C","FR4","CTLA.4","CD103",
                         "CD62L","GITR","CXCR3","CD44","PD.1","RORgT",
                         "CD127","ICOS","CCR6","CD38","Blimp1","CD69","KLRG1",
                         "T.bet","GATA.3","Helios","Neuropilin",
                         "CD25","IRF4")


# add the output to fcs.channel in the parameter file
temp <- channels$out1[channels$name %in% channels.of.interest]
temp[length(temp)] <- gsub(',', '', temp[length(temp)])
cat(paste(temp, collapse = ""))

# add the output to fcs.channel.label in the parameter file
temp <- channels$out2[channels$name %in% channels.of.interest]
temp[length(temp)] <- gsub(',', '', temp[length(temp)])
cat(paste(temp, collapse = ""))


# save the flowcytoscript_parameter.r file
# and run the advanced_flowcytoscript.r script