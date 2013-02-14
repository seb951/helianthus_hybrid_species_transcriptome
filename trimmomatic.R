#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from_cmd = as.numeric(args[1]) # Specify which sequences in "list_ind" file you want to align, directlty from the shell. Alternatively, you can do this from the "alignments" function itself.
to_cmd = as.numeric(args[2])


setwd("/media/seb_1TB/transposon") ### setup the working directory in the right format

trimo = function(list_ind= "reference/hybrid_species", from = 1, to = 1) {

individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0) # reference file
for(p in 1:nrow(individuals)) {individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])]} # proper format

	for(i in from:to) {
		if(length(grep("fq",individuals[i,1])) != 1)  #fastq files are sometimes zipped!
			{		
			#1.what sequencing file are there
			list = paste("ls -1 ",individuals[i,1],"/*fastq.gz   >list_",i, sep = ""); system(list)
			list = read.table(paste("list_",i,sep = ""), header = F, stringsAsFactors = F)
			
			#2. unzip all the fastq.gz files
			list = cbind(list,0)
				for(z in 1:nrow(list))
				{
				file = strsplit(list[z,1], split = "/")[[1]]
				list[z,2] = gsub(".gz","",file[length(file)])
				unzip = paste("gunzip -c ",list[z,1] ," >", list[z,2], sep = "")
				system(unzip)
				}
				
			#3.concatenate into 2 files
			seq_1 = gsub("$","_1.fq",individuals[i,5])
			seq_2 = gsub("$","_2.fq",individuals[i,5])

			cat_1 = paste("cat ", paste(list[c(1:(nrow(list)/2)),2], collapse = " ")," >data/",seq_1, sep = "")
			cat_2 = paste("cat ", paste(list[c(((nrow(list)/2)+1):nrow(list)),2], collapse = " ")," >data/",seq_2, sep = "")
			system(cat_1); system(cat_2)
			
			individuals[i,c(1,5)] = paste(individuals[i,c(1,5)],".fq", sep = "")
			}
			system(paste("rm *fastq list*")) #remove the temp files created during the unzipping
 			
###trimommatic #####
seq_1_in = gsub(".fq","_1.fq",individuals[i,5])
seq_2_in = gsub(".fq","_2.fq",individuals[i,5])
seq_1_out_p = gsub(".fq","_trim1.fq",individuals[i,5])
seq_2_out_p = gsub(".fq","_trim2.fq",individuals[i,5])
seq_1_out_up = gsub(".fq","_unpaired1.fq",individuals[i,5])
seq_2_out_up = gsub(".fq","_unpaired2.fq",individuals[i,5])
log = gsub(".fq","_trimo.log",individuals[i,5])

trim = paste("java -classpath /home/seb/Trimmomatic-0.22/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 6 -phred33 -trimlog data/",log," /media/seb_1TB/transposon/data/", seq_1_in," /media/seb_1TB/transposon/data/",seq_2_in," /media/seb_1TB/transposon/data/",seq_1_out_p," /media/seb_1TB/transposon/data/",seq_1_out_up," /media/seb_1TB/transposon/data/",seq_2_out_p," /media/seb_1TB/transposon/data/", seq_2_out_up," ILLUMINACLIP:reference/IlluminaContaminants.fa:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:10:15 MINLEN:36",sep = "")

system(trim)
}
}

###################
### run TRIMMOMATIC ###
###########e########
trimo(list_ind= "reference/hybrid_species", from = from_cmd, to = to_cmd)


