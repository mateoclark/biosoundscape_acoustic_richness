# concatenates multiple acoustic index CSV files
# Matthew Clark, November 1, 2022

workspace = "E:/archive/project/biosoundscape/AcousticIndicesCampaign2"

#out file name
filename <- "E:/archive/project/biosoundscape/bioscape_acoustic_indices_campaign1_240203.csv"

setwd(dir=workspace)

# load in files in directory
files = list.files(workspace,pattern=".csv", ignore.case=TRUE)

# load CSV files in a single table
table = data.frame()
for (i in 1:length(files)){
  print(i)
  data = read.csv(files[i])
  table = rbind(data, table)
}

# write out file
write.csv(table, file = filename,row.names = FALSE) 
