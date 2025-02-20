### Input: Rscript AGREPAssess_##.R [/path/to/indexfile] [MetaFile] [nThreads]
# Run in a folder containing all of the *.fa files resulting from AGREP/ugrep

# Indexfile is a tab-separated file containing [Designated Name for Probe Sequence] [RNAdesignation] [list of coordinates]
  # RNAdesignation: RT (Read-through) TS (Trans-spliced) (inter/intra)
  # For coordinates input, use the following format: "chr17:30700000-30800000;chr9:15000000-17000000"
    # Gene 1 and gene 2 coordinates are separated by a semicolon. Gene 2 coordinates are only required for TS RNA designation.

# MetaFile is a file of metadata that you want for each sample attached to counts. 
# Ex: SRR           Tissue  CellType  Treatment
# Ex: SRR123456789  Testes  Spermatid Treated  
  # The column names and any information you put on this table will be carried into the final outputs.
  # Suggestions for Tissue type, treatment, age, etc.
  # The first column must include ALL filenames used in the AGREP without the file extension.
    # This is required for the current implementation of MetaCounts.tsv, which matches the chimera name
    # And file/source name to count coincidences in a large matrix.

# nThreads needs to be specified, since it also must be provided in the batch submission.
  # If specified, this script will run pBLAT with that many threads. If not, normal BLAT, defaulting to n=1 threads.
  # nThreads > 8 will yield diminishing returns.


library(gtools)
library(stringr)
library(dplyr)

# Take in and Parse Inputs
vars.tmp <- commandArgs(trailingOnly = TRUE)
split.vars <- unlist(strsplit(vars.tmp, " "))
indexfile <- split.vars[1]


if (length(split.vars)>=2){
  MetaFile <- split.vars[2]
} 

if (length(split.vars)>=3){
  threads <- as.numeric(split.vars[3])
} else {
  threads <-1
}

indextable <- read.table(indexfile,stringsAsFactors = FALSE,sep = "\t")
IndexInputs<- ncol(indextable)
seqName <- indextable[,1]
if (IndexInputs > 1){
  RNAdesignation <- indextable[,2]
  coords <- indextable[,3]
}

# Load in Metatable. Contains RunID in leftmost column, and any other information user wants attached.
# Output file will keep provided header.
# Ex: SRR           Tissue  CellType  Treatment
# Ex: SRR123456789  Testes  Spermatid Treated
if (exists("MetaFile")){
Metatable <- read.table(MetaFile,stringsAsFactors = FALSE,sep = "\t",header = TRUE)
MetaInputs<- ncol(Metatable)
SRRnames <- Metatable[,1]
}

# Make folder to house results
foldername<- "AA_output"
if (file.exists("AA_output")){ 
  # if folder exists, remove contents
  # else, make folder.
system("rm -r AA_output")
system(paste("mkdir",foldername))
} else {
  system(paste("mkdir",foldername))
}

# Define RT and TS filenames
filename_RT <- paste0(foldername,"/","PosHits_RT.fa")
filename_TS <- paste0(foldername,"/","PosHits_TS.fa")

# Combine all files into one
system(paste0("awk '(NR == 1) || (FNR > 1)' *fastq*.fa > ",foldername,"/all.fa"))
# Sort all RT reads into RT.fa, TS into TS.fa
system(paste0("grep -A1 _RT_ ",foldername,"/all.fa > ", filename_RT))
system(paste0("grep -A1 _TS_ ",foldername,"/all.fa > ", filename_TS))
# Remove the large all.fa
system(paste0("rm ",foldername,"/all.fa"))

# Define function to read BLAT .psl outputs in simple format
read.psl <- function(x) {
  # Ensure x is a character string of .psl file path
  psltable <- try(read.table(x,header = FALSE, sep = "\t", stringsAsFactors = FALSE, skip = 5),silent=TRUE)
  if(inherits(psltable, "try-error")){
    return(NULL)
  }
  else{
    colnames(psltable) <- c("match","mismatch","rep.match","N's","QgapCount","QgapBases","TgapCount","TgapBases","strand","Qname","Qsize","Qstart","Qend","Tname","Tsize","Tstart","Tend","BlockCount","blockSizes","qStarts","tStarts")
    return(psltable)
  }
}

# Set RT and TS filenames
out_filename_RT <- paste(foldername,"/","RT_out.psl",sep = "")
reps_filename_RT <- paste(foldername,"/","RT_reps.psl",sep = "")
psr_filename_RT <- paste(foldername,"/","RT_out.psr",sep = "")

out_filename_TS <- paste(foldername,"/","TS_out.psl",sep = "")
reps_filename_TS <- paste(foldername,"/","TS_reps.psl",sep = "")
psr_filename_TS <- paste(foldername,"/","TS_out.psr",sep = "")

## Submit BLAT command to console
# **pBLAT** submission to mimic default web blat: 
# Note that these are local files are publicly available: https://genome.ucsc.edu/FAQ/FAQblat.html#blat3
# pBLAT is available here: https://github.com/icebert/pblat

# New commands without patches:
system(paste("/scratch/jme5fe/pblat-2.5-2/pblat /sfs/weka/scratch/jme5fe/genomes/hg38.2bit", filename_RT,paste0("-threads=",threads),"-stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -ooc=/scratch/jme5fe/blatSuite36/11.ooc", out_filename_RT,sep = " "))
system(paste("/scratch/jme5fe/pblat-2.5-2/pblat /sfs/weka/scratch/jme5fe/genomes/hg38.2bit", filename_TS,paste0("-threads=",threads),"-stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -ooc=/scratch/jme5fe/blatSuite36/11.ooc", out_filename_TS,sep = " "))

## Pick most likely true hit/result from BLAT output
if (file.exists(out_filename_RT)){
system(paste("/project/hlilab/justin/software/blatSuite36/pslReps -singleHit",out_filename_RT, reps_filename_RT, psr_filename_RT,sep = " "))
}
if (file.exists(out_filename_TS)){
system(paste("/project/hlilab/justin/software/blatSuite36/pslReps",out_filename_TS, reps_filename_TS, psr_filename_TS,sep = " "))
}

# Parse Coords into Coord Matrix
coord_mat<-as.data.frame(matrix(data = NA, nrow = length(coords),ncol = 7))
colnames(coord_mat)<-c("SeqName","Chrom1","Start1","End1","Chrom2","Start2","End2")
coord_mat$SeqName<-as.matrix(seqName)

# Coords for RT
coord_mat$Chrom1[RNAdesignation=="RT"]<-gsub(":.*","",coords[RNAdesignation=="RT"])
coord_mat$Start1[RNAdesignation=="RT"]<-as.numeric(gsub("-.*","",gsub(".*:","",coords[RNAdesignation=="RT"])))
coord_mat$End1[RNAdesignation=="RT"]<-as.numeric(gsub(".*-","",coords[RNAdesignation=="RT"]))
coord_mat$Chrom2[RNAdesignation=="RT"]<-NA
coord_mat$Start2[RNAdesignation=="RT"]<-NA
coord_mat$End2[RNAdesignation=="RT"]<-NA

# Coords for TS
TS_coord1<-gsub(";.*","",coords[RNAdesignation=="TS"])
TS_coord2<-gsub(".*;","",coords[RNAdesignation=="TS"])
coord_mat$Chrom1[RNAdesignation=="TS"]<-gsub(":.*","",TS_coord1)
coord_mat$Start1[RNAdesignation=="TS"]<-as.numeric(gsub("-.*","",gsub(".*:","",TS_coord1)))
coord_mat$End1[RNAdesignation=="TS"]<-as.numeric(gsub(".*-","",TS_coord1))
coord_mat$Chrom2[RNAdesignation=="TS"]<-gsub(":.*","",TS_coord2)
coord_mat$Start2[RNAdesignation=="TS"]<-as.numeric(gsub("-.*","",gsub(".*:","",TS_coord2)))
coord_mat$End2[RNAdesignation=="TS"]<-as.numeric(gsub(".*-","",TS_coord2))

# Read in .psl filtered BLAT results
if (file.exists(reps_filename_RT)){
RT_table <- read.psl(reps_filename_RT)
} else {
RT_table<-NULL
}
if (file.exists(reps_filename_TS)){
TS_table <- read.psl(reps_filename_TS)
} else {
TS_table <- NULL
}

## Transform into usable format.
# SeqName | Qname | Qstart1 | Qend1 | Tchr1 | Tstart1 | Tend1 | Qstart2 | Qend2 | Tchr2 | Tstart2 | Tend2 
# Create RT_table. Include second coordinate as NA to match w/ TS dimensions.
if (!is.null(RT_table)){
RT_table$SeqName <- gsub("_(RT|TS).*","",RT_table$Qname)
RT_table<-cbind(RT_table[,c(22,10,12:14,16,17)],matrix(NA,nrow = nrow(RT_table),ncol = 5))
colnames(RT_table)<-c("SeqName","Qname","Qstart1","Qend1","Tname1","Tstart1","Tend1","Qstart2","Qend2","Tname2","Tstart2","Tend2")
}

# Create TS_table.
if (!is.null(TS_table)){
TS_unique<-unique(TS_table$Qname)
TS_table_mod<-matrix(NA,nrow = length(TS_unique),ncol = 11)
for (i in 1:length(TS_unique)){
  # Pull out Qname + relevant columns
  TS_temp<-TS_table[TS_table$Qname==TS_unique[i],c(10,12:14,16:17)]
  TS_temp<-TS_temp[order(TS_temp$Qstart),]
  TS_table_mod[i,]<-as.matrix(unlist(c(TS_temp[1,],TS_temp[2,-1])))
}
TS_table<-data.frame(TS_table_mod)
rm(TS_table_mod)
colnames(TS_table)<-c("Qname","Qstart1","Qend1","Tname1","Tstart1","Tend1","Qstart2","Qend2","Tname2","Tstart2","Tend2")
TS_table$SeqName <- gsub("_(RT|TS).*","",TS_table$Qname)
TS_table<-TS_table[,c(12,1:11)]
}

## Next steps: compare Tchr/start/end/ to coords, assign hits/no hits
## and create: counts for each, table of positive alignments, .fa for pos hits only.

# For each entry in RT & TS tables (each result for BLAT sequence)...
# Compare to coords for its "index" or searched sequence.
# Assign 1) TruePos -- all match, 2) PartPos -- match+mismatch, 3) Neg -- NoMatch

# RT_assess
if (!is.null(dim(RT_table))){ # If RT table is not empty
RT_assess<-as.data.frame(matrix(data = NA, nrow = nrow(RT_table),ncol = 4))
colnames(RT_assess)<-c("Chr1","Tstart1","Tend1","Assess") 
for (i in 1: nrow(RT_table)){
  RT_assess$Chr1[i]<-RT_table$Tname1[i]==coord_mat$Chrom1[coord_mat$SeqName==RT_table$SeqName[i]]
  RT_assess$Tstart1[i]<-as.numeric(as.character(RT_table$Tstart1[i]))>=coord_mat$Start1[coord_mat$SeqName==RT_table$SeqName[i]]
  RT_assess$Tend1[i]<-as.numeric(as.character(RT_table$Tend1[i]))<=coord_mat$End1[coord_mat$SeqName==RT_table$SeqName[i]]
}
}

# TS_assess
if (!is.null(dim(TS_table))){ # If TS table is not empty
  TS_assess<-as.data.frame(matrix(data = NA, nrow = nrow(TS_table),ncol = 7))
  colnames(TS_assess)<-c("Chr1","Tstart1","Tend1","Chr2","Tstart2","Tend2","Assess") 
  for (i in 1: nrow(TS_table)){
    TS_assess$Chr1[i]<-TS_table$Tname1[i]==coord_mat$Chrom1[coord_mat$SeqName==TS_table$SeqName[i]]
    TS_assess$Tstart1[i]<- as.numeric(as.character(TS_table$Tstart1[i]))>=coord_mat$Start1[coord_mat$SeqName==TS_table$SeqName[i]]
    TS_assess$Tend1[i]<- as.numeric(as.character(TS_table$Tend1[i]))<=coord_mat$End1[coord_mat$SeqName==TS_table$SeqName[i]]
    TS_assess$Chr2[i]<-TS_table$Tname2[i]==coord_mat$Chrom2[coord_mat$SeqName==TS_table$SeqName[i]]
    TS_assess$Tstart2[i]<- as.numeric(as.character(TS_table$Tstart2[i]))>=coord_mat$Start2[coord_mat$SeqName==TS_table$SeqName[i]]
    TS_assess$Tend2[i]<- as.numeric(as.character(TS_table$Tend2[i]))<=coord_mat$End2[coord_mat$SeqName==TS_table$SeqName[i]]
  }
}

# Assess_outputs
PartPos <- function(x) {
  # Ensure x is a character string of .psl file path
  x<-as.data.frame(x)
  out <- any(all(x[1:3,]),all(x[4:6,])) & any(!all(x[1:3,]),!all(x[4:6,]))
  return(out)
}

if (exists("RT_assess")){
RT_assess$Assess[apply(!RT_assess[,1:3],1,any)]<-"FalsePos"
RT_assess$Assess[apply(RT_assess[,1:3],1,all)]<-"TruePos"
RT_assess$Assess[is.na(RT_assess$Assess)]<-"NA_Other"
}

if (exists("TS_assess")){
TS_assess$Assess[apply(!TS_assess[,1:6],1,any)]<-"FalsePos"
TS_assess$Assess[apply(TS_assess[,1:6],1,all)]<-"TruePos"
TS_assess$Assess[apply(TS_assess[,1:6],1,PartPos)]<-"PartPos"
TS_assess$Assess[is.na(TS_assess$Assess)]<-"NA_Other"
}

# Counts for True/Partial/Neg. By query.
if (exists("RT_assess")){
RT_counts<-as.data.frame(matrix(NA,nrow = length(unique(RT_table[,1])),ncol = 6))
colnames(RT_counts)<-c("SeqName","Total","TruePos","Partial","FalsePos","NA_Other")
RT_counts[,1]<-unique(RT_table[,1])
for (i in 1:nrow(RT_counts)){
    temp_index<-RT_table$SeqName==RT_counts$SeqName[i]
  RT_counts$Total[i]<-sum(temp_index)
  RT_counts$TruePos[i]<-sum(RT_assess$Assess[temp_index]=="TruePos")
  RT_counts$Partial[i]<-sum(RT_assess$Assess[temp_index]=="PartPos")
  RT_counts$FalsePos[i]<-sum(RT_assess$Assess[temp_index]=="FalsePos")
  RT_counts$NA_Other[i]<-sum(RT_assess$Assess[temp_index]=="NA_Other")
}
}

if (exists("TS_assess")){
TS_counts<-as.data.frame(matrix(NA,nrow = length(unique(TS_table[,1])),ncol = 6))
colnames(TS_counts)<-c("SeqName","Total","TruePos","Partial","FalsePos","NA_Other")
TS_counts[,1]<-unique(TS_table[,1])
for (i in 1:nrow(TS_counts)){
  temp_index<-TS_table$SeqName==TS_counts$SeqName[i]
  TS_counts$Total[i]<-sum(temp_index)
  TS_counts$TruePos[i]<-sum(TS_assess$Assess[temp_index]=="TruePos")
  TS_counts$Partial[i]<-sum(TS_assess$Assess[temp_index]=="PartPos")
  TS_counts$FalsePos[i]<-sum(TS_assess$Assess[temp_index]=="FalsePos")
  TS_counts$NA_Other[i]<-sum(TS_assess$Assess[temp_index]=="NA_Other")
}
}

# Define output based on what exists.
if (!is.null(RT_table) & !is.null(TS_table)){ #In case submits only TS or RT.
  counts_table<-rbind(RT_counts,TS_counts)
} else if (!is.null(RT_table) & is.null(TS_table)){# only RT
  counts_table<-RT_counts
} else if (is.null(RT_table) & !is.null(TS_table)){ # only TS
  counts_table<-TS_counts
}

# Table for True alignments
if (!is.null(RT_table) & !is.null(TS_table)){ # Just in case submits only TS or RT.
  PosTable<-rbind(RT_table[RT_assess$Assess=="TruePos",],TS_table[TS_assess$Assess=="TruePos",])
} else if (!is.null(RT_table) & is.null(TS_table)){# only RT
  PosTable<-RT_table[RT_assess$Assess=="TruePos",]
} else if (is.null(RT_table) & !is.null(TS_table)){ # only TS
  PosTable<-TS_table[TS_assess$Assess=="TruePos",]
}

write.table(RT_table,file = paste0(foldername,"/RT_table.tsv"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(TS_table,file = paste0(foldername,"/TS_table.tsv"),sep = "\t",row.names = F,col.names = T,quote = F)

# Pos Counts Sorted by Tissue/Sample/Metadata.
# Whole Metadata table merged with ALL chimeras as columns. Counts in matrix.
if (exists("MetaFile")){
MetaCounts<-matrix(0,nrow=nrow(Metatable),ncol=length(seqName))
MetaCounts<-cbind(Metatable,MetaCounts)
colnames(MetaCounts)<-c(colnames(Metatable),seqName)
rownames(MetaCounts)<-MetaCounts[,1]
for (i in 1:nrow(PosTable)){
  # Regex that will likely need to be adjusted for use by others. Sorry.
  SRR<-gsub("_pass","",gsub(".*(?:RT|TS)_(.*)\\.fastq.*",replacement="\\1",PosTable$Qname[i]))
  MetaCounts[as.character(SRR),PosTable$SeqName[i]]<-MetaCounts[as.character(SRR),PosTable$SeqName[i]]+1
}
}


# List of Partial and FalsePos Reads and NA reads
if (!is.null(RT_table) & !is.null(TS_table)){ # Just in case submits only TS or RT.
  FPTable<-rbind(RT_table[RT_assess$Assess=="FalsePos",],TS_table[TS_assess$Assess=="FalsePos",])
  PPTable<-rbind(RT_table[RT_assess$Assess=="PartPos",],TS_table[TS_assess$Assess=="PartPos",])
  NATable<-rbind(RT_table[RT_assess$Assess=="NA_Other",],TS_table[TS_assess$Assess=="NA_Other",])
} else if (!is.null(RT_table) & is.null(TS_table)){# only RT
  FPTable<-RT_table[RT_assess$Assess=="FalsePos",]
  PPTable<-RT_table[RT_assess$Assess=="PartPos",]
  NATable<-RT_table[RT_assess$Assess=="NA_Other",]
} else if (is.null(RT_table) & !is.null(TS_table)){ # only TS
  FPTable<-TS_table[TS_assess$Assess=="FalsePos",]
  PPTable<-TS_table[TS_assess$Assess=="PartPos",]
  NATable<-TS_table[TS_assess$Assess=="NA_Other",]
}

# Write counts_table
write.table(counts_table,file = paste0(foldername,"/countsTable.tsv"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(PosTable,file = paste0(foldername,"/PosTable.tsv"),sep = "\t",row.names = F,col.names = T,quote = F)
write.table(PosTable$Qname,file = paste0(foldername,"/Pos_sequences.fa"),sep = "\t",row.names = F,col.names = F,quote = F)

# Print other tables for QC if desired
#write.table(FPTable,file = paste0(foldername,"/FPTable.tsv"),sep = "\t",row.names = F,col.names = F,quote = F)
#write.table(PPTable,file = paste0(foldername,"/PPTable.tsv"),sep = "\t",row.names = F,col.names = F,quote = F)
#write.table(NATable,file = paste0(foldername,"/NATable.tsv"),sep = "\t",row.names = F,col.names = F,quote = F)

#write.table(FPTable$Qname,file = paste0(foldername,"/FP_sequences.fa"),sep = "\t",row.names = F,col.names = F,quote = F)
#write.table(PPTable$Qname,file = paste0(foldername,"/PP_sequences.fa"),sep = "\t",row.names = F,col.names = F,quote = F)
#write.table(NATable$Qname,file = paste0(foldername,"/NA_sequences.fa"),sep = "\t",row.names = F,col.names = F,quote = F)

#if (exists("RT_assess")){
#write.table(RT_assess,file = paste0(foldername,"/RT_assess.tsv"),sep = "\t",row.names = F,col.names = F,quote = F)
#}
# if (exists("TS_assess")){
#   write.table(TS_assess,file = paste0(foldername,"/TS_assess.tsv"),sep = "\t",row.names = F,col.names = F,quote = F)
# }

if (exists("MetaFile")){
  write.table(MetaCounts,file = paste0(foldername,"/MetaCounts.tsv"),sep = "\t",row.names = F,col.names = T,quote = F)
}
