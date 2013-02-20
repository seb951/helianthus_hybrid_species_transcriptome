#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from = as.numeric(args[1]) # Specify which sequences in "list_ind" file you want to count stats on. You should do all that you aligned otherwise, bug!
to = as.numeric(args[2])

setwd("/media/seb_1TB/transposon") ### setup the working directory in the right format

#####################
###stats alignment function###
####################
trimo_stats = function(list_ind= "reference/hybrid_species", from = 1, to = 1) {

individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0,0,0,0) # reference file
for(p in 1:nrow(individuals)) {
individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])];
if(length(grep("fq",individuals[p,1])) != 1) individuals[p,c(1,5)] = paste(individuals[p,c(1,5)],".fq", sep = "")
} # proper format

system("grep 'Pairs' trim.log >trim.log.perc_filter") # trim.log was created by running the trim command above with nohup and a log file #
percentage = read.delim("trim.log.perc_filter", header = F,sep = " ", stringsAsFactors = F)
colnames(individuals)[6:8] = c("perc_filter","mean_length","raw_seq_num")
individuals[,6] = percentage[,8]

	for(i in from:to) {
	
seq_1_out_p = gsub(".fq","_trim1.fq",individuals[i,5])
seq_1 = gsub(".fq","_1.fq",individuals[i,1],fixed = T)
system(paste("awk 'length > 36' data/",seq_1_out_p," | wc >wc.out",sep = "")) #find out mean sequence length of the trimmed sequences
system(paste("wc -l ",seq_1," >wc_",i,sep = "")) #how many raw sequences	

wc = read.table(paste("wc_",i,sep = ""), header = F, stringsAsFactors = F)[1,1]
wc.out = read.delim("wc.out", header = F,sep = " ", stringsAsFactors = F)

individuals[i,7] = (wc.out[1,3] - wc.out[1,1]) / wc.out[1,1]
individuals[i,8] = wc/4
							}
system("rm wc* trim.log.perc_filter")							

write.table(individuals[,c(6:8)],"results/trim_meanlength", row.names = F, col.names = T)

#return(individuals)
}


##########################################
### count the number of aligned sequences (samtools idxstats) ###
##########################################
idxstats = function(list_ind= "reference/hybrid_species",total_idxstats = NULL, from = 1, to = 1) {

total_seq = read.table("results/trim_meanlength",header = T, stringsAsFactors = F);total_seq = cbind(total_seq,0,0);colnames(total_seq)[4:5] = c("aligned_total","perc_aligned")
individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0) # reference file
for(p in 1:nrow(individuals))
	{
	individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])];
	if(length(grep("fq",individuals[p,1])) != 1) individuals[p,5] = paste(individuals[p,5],".fq",sep = "")
	} # proper format to further process with idxstats

	for(i in from:to) {

	#step 1: generate the alignment statistics#
	command_idxstats = paste("samtools idxstats alignments/",individuals[i,5],".sorted.bam"," >idxstats_results",  sep = "")
	system(command_idxstats)
	
	#step 2: read the covage and the alignment stats and then update the result_matrix file
	idxstats_results = read.delim("idxstats_results", stringsAsFactors = F, header = F) #don't read last line.
	
	total_seq[i,4] = sum(idxstats_results[,3]); total_seq[i,5] = total_seq[i,4] / (total_seq[i,3] * 2) 
	
	if(i == 1) {total_idxstats =  idxstats_results[1:(nrow(idxstats_results)-1),c(1:3)]; colnames(total_idxstats) = c("name","length",paste(individuals[i,5],"_number_reads",sep = "") )}
	if(i > 1) {total_idxstats = cbind(total_idxstats,idxstats_results[1:(nrow(idxstats_results)-1),3]); colnames(total_idxstats)[ncol(total_idxstats)] = paste(individuals[i,5],"_number_reads",sep = "")}

	}
	system("rm idxstats_results")
	write.table(total_idxstats,"results/coverage_per_gene_individuals.txt",row.names = F,col.names = T)
	write.table(total_seq,"results/number_total_aligned.txt",row.names = F)
}

############
### running ###
############
trimo_stats(list_ind= "reference/hybrid_species", from = from_cmd, to = to_cmd)
idxstats(from = from_cmd, to = to_cmd) ### count the number of aligned sequences



############
### extra stats you may want ###
############



if(file.exists("results/total_idxstats.txt")) {
total_idxstats = read.table("results/total_idxstats.txt",header = T, stringsAsFactors = F)

mean_med = NULL
for(i in 3:ncol(total_idxstats))
{
mean_med = rbind(mean_med,paste(signif(mean(total_idxstats[,i]),4)," (",median(total_idxstats[,i]),")",sep = ""))
}

write.table(mean_med,"results/mean_med", row.names = F, col.names = F, quote = F)
}




