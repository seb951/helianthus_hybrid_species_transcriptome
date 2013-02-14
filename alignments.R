#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from_cmd = as.numeric(args[1]) # Specify which sequences in "list_ind" file you want to align, directlty from the shell. Alternatively, you can do this from the "alignments" function itself.
to_cmd = as.numeric(args[2])

#############
### functions ###
#############
#alignments(from =  from_cmd, to =  from_cmd)) ### align  with bwa


################
###requirements #####
################

### YOU NEED SAMTOOLS AND BWA INSTALLED!!!!!! ###
setwd("/media/seb_1TB/transposon") ### setup the working directory in the right format
#system("bwa index reference/HA412_trinity_noAltSplice_400bpmin.fa") ### indexing the reference BWA
#system("samtools faidx reference//HA412_trinity_noAltSplice_400bpmin.fa")### indexing the reference SAMTOOLS

###########################
### defining the alignments function ###
###########################
alignments = function(list_ind= "reference/hybrid_species", from = 1, to = 1) {

individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0) # reference file
for(p in 1:nrow(individuals)) {
individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])];
if(length(grep("fq",individuals[p,1])) != 1) individuals[p,c(1,5)] = paste(individuals[p,c(1,5)],".fq", sep = "")
} # proper formatting

for(i in from:to){
### BWA ALIGNMENTS
bwa_aln1 = bwa_aln2 = bwa_sampe = "ls"

	if((individuals[i,3] == "ill") & (individuals[i,4] == "sanger")) # new illumina quality format specified in the "list_ind" file
{bwa_aln1 = paste("bwa aln -t 8 -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_trim1.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_trim1.fq",individuals[i,5], fixed = T), ".sai", sep = "");
bwa_aln2 = paste("bwa aln -t 8 -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_trim2.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_trim2.fq",individuals[i,5], fixed = T), ".sai", sep = "")}

	if((individuals[i,3] == "ill") & (individuals[i,4] == "Ill1.3")) # old illumina quality format specified in the "list_ind" file
{bwa_aln1 = paste("bwa aln -t 8 -I -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_trim1.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_trim1.fq",individuals[i,5], fixed = T), ".sai", sep = "");
bwa_aln2 = paste("bwa aln -t 8 -I -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_trim2.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_trim2.fq",individuals[i,5], fixed = T), ".sai", sep = "")}

	if(individuals[i,3] == "ill")
bwa_sampe = paste("bwa sampe reference/HA412_trinity_noAltSplice_400bpmin.fa alignments/",sub(".fq","_trim1.fq",individuals[i,5], fixed = T), ".sai"," alignments/",sub(".fq","_trim2.fq",individuals[i,5], fixed = T), ".sai ",sub(".fq","_trim1.fq",individuals[i,1], fixed = T)," ", sub(".fq","_trim2.fq",individuals[i,1], fixed = T)," >alignments/",individuals[i,5], ".sam", sep = "")

sam_view = paste("samtools view -bS -o alignments/",individuals[i,5], ".bam", " alignments/",individuals[i,5], ".sam", sep = "")

sam_sort = paste("samtools sort alignments/",individuals[i,5], ".bam ","alignments/",individuals[i,5], ".sorted",sep = "")

sam_index = paste("samtools index alignments/",individuals[i,5],".sorted.bam",sep = "") #index bam files

sam_rm = paste("rm alignments/",individuals[i,5], ".bam ", "alignments/",individuals[i,5], ".sam", sep = "")

system(bwa_aln1)
system(bwa_aln2) 
system(bwa_sampe)
system(sam_view)
system(sam_sort)
system(sam_index)
system(sam_rm)
}
system("rm alignments/*sai")
}

############
### running ###
############
alignments(from =  from_cmd, to =  to_cmd) ### align 

 

 
