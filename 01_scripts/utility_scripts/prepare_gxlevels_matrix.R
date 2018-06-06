# Prepare input from eXpress to input into edgeR
setwd("~/Documents/bernatchez/01_Sfon_projects/04_Sfon_eQTL/Simple_reads_to_counts_2017-04-05/05_gx_levels")

# Identify file names
files <- list.files(
  #path = "05_gx_levels", 
  pattern = "*_results.xprs")

# temp for troubleshooting
#files <- files[1:2]


# Read in files
current.name <- NULL;
expr.list <- list()

for(i in 1:length(files)){
  
  # obtain the name of the current file
  current.name <- files[i]
  print(current.name)
  current.name <- gsub(x = current.name, pattern = "_trimmed.fastq.gz.bowtie2.sorted.bam_results.xprs",replacement = "")
  print(current.name)
  
  # put the current file into a list
  expr.list[[i]] <- read.table(file = files[i], header = F, sep = "\t", )[-1,c(2,8)]
  colnames(expr.list[[i]]) <- c(paste(current.name,".transcript", sep = ""), 
                                paste(current.name,".eff.counts", sep = ""))
  
  # sort the current list accession
  sort.column.name <- colnames(expr.list[[i]])[1]
  
  expr.list[[i]] <- 
    expr.list[[i]][with(expr.list[[i]], order(expr.list[[i]][sort.column.name])), ]
  
}

str(expr.list)

# this results in a list of all necessary info


# Build final matrix file
# create variable with first column as the transcript ID
out.matrix <- expr.list[[1]][1]

for(f in 1:length(expr.list)){
  out.matrix <- cbind(out.matrix, expr.list[[f]][2])
}

head(out.matrix)

colnames(out.matrix)[1] <- "transcript.id"

write.csv(x = out.matrix, file = "out.matrix.csv", row.names = F)
