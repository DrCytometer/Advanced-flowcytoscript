# get channels script
# run this to get the channel names for the flowcytoscript parameter file

library(flowCore)
library(dplyr)
library( RColorBrewer )

# source parameters

param.filename <- "./flowcytoscript_parameter.r"
source( param.filename )

# need to modify for fcs or csv options
# add option to generate fcs.condition.label

flowFrame <- read.FCS(list.files(data.dir, "\\.fcs$", full.names = TRUE)[1], truncate_max_range = FALSE)

channels <- data.frame(name = unname(pData(parameters(flowFrame))$name), 
                      desc = unname(pData(parameters(flowFrame))$desc))

for(i in 1:dim(channels)[1] ){
  channels$out1[i] = paste0("\"", channels$name[i], "\", #", channels$desc[i], "\n")
  channels$out2[i] = paste0("\"", channels$name[i], "\" = \"", channels$desc[i], "\",\n") 
}

descs.filtered <- channels$desc[!is.na(channels$desc) & channels$desc != '-'] 
channels.filtered <- dplyr::filter(channels, desc %in% descs.filtered)
cat('"', paste0(channels.filtered$desc, collapse = '","'), '"', sep = "")

# copy the previous output into the channels.of.interest and delete unnecessary markers

channels.of.interest = c("Foxp3","T-bet","IRF4","Ki67","CCR9","CD95","Ly-6C","CD103",
                         "CTLA-4","CD62L","CXCR6","CCR2",
                         "CD44","ICOS","RORgT","PD-1","CXCR3","TNFRII",
                         "CD69","CD25","ST2","GATA-3","Neuropilin",
                         "KLRG1","Helios")


# add the output to fcs.channel in the parameter file
temp <- channels$out1[channels$desc %in% channels.of.interest]
temp[length(temp)] <- gsub(',', '', temp[length(temp)])
cat(paste(temp, collapse = ""))

# add the output to fcs.channel.label in the parameter file
temp <- channels$out2[channels$desc %in% channels.of.interest]
temp[length(temp)] <- gsub(',', '', temp[length(temp)])
cat(paste(temp, collapse = ""))


# save the flowcytoscript_parameter.r file
# and run the advanced_flowcytoscript.r script