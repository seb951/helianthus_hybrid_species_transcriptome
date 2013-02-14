#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
species = args[1] # can specify which species you want to align, directlty from the shell. Alternatively, you can do this from the "snp" function itself.
cpu = as.numeric(args[2])
cpu = args[3] #either snp OR snp_stats
################
###requirements #####
################

### YOU NEED SAMTOOLS and BCFTOOLS. You should also have vcfutils.pl in your "Rcode/" directory ###
#you should run snp() species per species, but not simultaneously. Otherwise there will be some problems with overwriting files. Uncomment the 
#you should run snp_stats() once you have ran snp for all 3 species.



###########################
### defining the SNP calling function ###
###########################
snp = function(list_ind= "reference/hybrid_species",user = "seb", species = "deserticola",genes = 51468,cpu = 16,wd = "/media/seb_1TB/transposon") {

setwd(wd) ### setup the working directory in the right format with the proper subdirectories.

###create "split" object for mpileup if it doesn't exist###
system("grep '>' reference/HA412_trinity_noAltSplice_400bpmin.fa | awk '{print $1}' >reference/temp")
temp = read.table("reference/temp", stringsAsFactors = F, header = F);temp[,1] = cbind(gsub(">","",temp[,1]),"A")

if(file.exists("reference/split") == F){
for(i in 1: nrow(temp)) {
cmd = paste("awk 'NR == ", i*2,"' reference/HA412_trinity_noAltSplice_400bpmin.fa | wc -m >reference/temp_wc", sep = ""); system(cmd)
temp[i,2] = paste("1-",read.table("reference/temp_wc",stringsAsFactors = F, header = F)-1,sep = "")
if(i %% 10000 == 0) print(paste(i,"of 51k, Time is:",Sys.time()))
}
split = as.data.frame(paste(temp[,1],temp[,2],sep = ":"));colnames(split) = "V1";rm(temp) #split file for mpileup
} else split = read.table("reference/split", stringsAsFactors = F)

###load individuals' names
individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0) # reference file
for(p in 1:nrow(individuals))
{individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])];
if(length(grep("fq",individuals[p,1])) == 1) individuals[p,5] = paste(individuals[p,5],".sorted.bam",sep = "");
if(length(grep("fq",individuals[p,1])) != 1) individuals[p,5] = paste(individuals[p,5],".fq.sorted.bam",sep = "")} # proper format to further process with idxstats

if(species == "deserticola") {individuals = individuals[regexpr("des",individuals[,2])>0,]}
if(species == "paradoxus") {individuals = individuals[regexpr("par",individuals[,2])>0,]}
if(species == "anomalus") {individuals = individuals[regexpr("ano",individuals[,2])>0,]}

all_final = paste("alignments/", individuals[,5], collapse = " ", sep = "")

############
### mpileup ###
############
#To run mpileup gene per gene, this way you can parallelize.#
for(j in 1:genes)#2000 genes only
{
system(paste("ps -u",user,"| grep 'samtools' | wc -l >prog"))
mpileup = paste("samtools mpileup -C50 -I -ugf reference/HA412_trinity_noAltSplice_400bpmin.fa ", all_final, " -r ",split[j,1]," | bcftools view -cvg - > mpileup/variants",j,".raw.vcf", sep = "")
Sys.sleep(ifelse(read.table("prog") <= cpu,0, (read.table("prog") - cpu)^4    ))
mpileup_exe = paste("mpileup/cmd/mpileup_",j,sep = "")
write.table(mpileup,mpileup_exe, row.names = F, col.names = F, quote = F)
system(paste("chmod +x",mpileup_exe))
system(paste("nohup ./", mpileup_exe, ">log&",sep = ""))
if(j %% 100 == 0) print(paste("Calling SNPs of gene number",j,"(out of  ",genes, ") The time is:",Sys.time()))
}

#cat the variants_j.raw.vcf files. 
Sys.sleep(60) #tidy up the processes running, once you are done
system("echo -n '' > mpileup/0_51k_variants.raw.vcf") #final file
system("cat mpileup/variants*  >>mpileup/0_51k_variants.raw.vcf") # cat verything
system("grep '#' -v  mpileup/0_51k_variants.raw.vcf >mpileup/0_51k_variants.clean.vcf") #tidy up
system("rm mpileup/variants*")
system("rm mpileup/cmd/*")


###Pre-cleaning step###
	system("wc -l mpileup/0_51k_variants.clean.vcf >mpileup/wordcount1")
	wordcount = read.table("mpileup/wordcount1")
	system("awk 'NR == 27' mpileup/0_51k_variants.raw.vcf  >mpileup/header")
	header = read.delim("mpileup/header", header = F, stringsAsFactors = F)
	header = gsub("alignments/","",header[10:length(header)], fixed = T)
	header = gsub(".fq.sorted.bam","",header, fixed = T)
	header = gsub(".sorted.bam","",header, fixed = T)
	header = gsub("_R1.fq.bz2","",header, fixed = T)
	header = gsub("/SciBorg/array0/renaut/repeatability/","",header, fixed = T)
	header = gsub("/SciBorg/array0/renaut/speciation_islands_individuals/","", header, fixed = T)
	write.table(t(as.matrix(c("reference", "position","reference_annuus_allele",header))),"results/snp_table",row.names = F, col.names = F, quote = F, sep = " ")

	temp = NULL
	mis = NULL
	con = file("mpileup/0_51k_variants.clean.vcf")
	open(con)
	for(i in 1:as.numeric(wordcount[1]))
	{
		x = readLines(con,1) #you will now read the mega big 0-51k_variants.clean.vcf file line by line in R. 
		xx = strsplit(x, split = "\t")[[1]]
		xxx = xx[(length(xx)-(length(header) - 1)):length(xx)]
		
		for(q in 1:length(xxx))	# this is to replace low quality calls (maximum Phred-scaled genotype likelihoods below 20). ie. essentially, you need at least 2 reads to call a SNP!
		{
		temp = strsplit(xxx[q], split = ":|,")[[1]]
		#if(max(as.numeric(temp[length(temp)])) < 20) temp[1] = "XX" else if(min(as.numeric(temp[2:4])) == temp[2]) temp[1] = "RR" else if(min(as.numeric(temp[2:4])) == temp[3]) temp[1] = "AR" else if(min(as.numeric(temp[2:4])) == temp[4]) temp[1] = "AA" #old way of doing things, may not be appropriate...
		if(max(as.numeric(temp[length(temp)])) < 20) xxx[q] = "XX" else {xxx[q] = temp[1]; xxx[q] = gsub("/","",xxx[q],fixed = T); xxx[q] = gsub("0","R",xxx[q],fixed = T); xxx[q] = gsub("1","A",xxx[q],fixed = T)} #new way of doing things...
		}
		xxx = gsub("R",xx[4],xxx)
		xxx = gsub("A",xx[5],xxx)
		xxx = c(xx[c(1,2,4)],xxx)
		if(nchar(xx[5]) > 1) xxx[4:length(xxx)] = rep("XX",length(xxx)-3)
		
		if((length(xxx[xxx == "XX"]) / length(xxx[4:length(xxx)])) < 0.4)	cat(t(as.matrix(xxx)),file = "results/snp_table",append = T, fill = F, "\n") else mis = c(mis, i) #only cat loci with less than 40% missing data. 
	#	if(regexpr(xxx[3],paste(xxx[4:21],collapse = "")) < 0) temp = rbind(temp,xxx)
		if(i %% 1000 == 0) print(paste(i,"of", wordcount[1], Sys.time()))
		}
close(con)

if(species == "deserticola") system("cat results/snp_table | sed 's/[ \t]*$//' >results/snp_table_deserticola")
if(species == "paradoxus")system("cat results/snp_table | sed 's/[ \t]*$//' >results/snp_table_paradoxus")
if(species == "anomalus") system("cat results/snp_table | sed 's/[ \t]*$//' >results/snp_table_anomalus")

system("rm results/snp_table reference/temp*") ### rm snp_table and other temp files 
system("rm mpileup/0_51* mpileup/header mpileup/wordcount1 log prog ") ### rm undesirable outputs
}



#########################
### more stats, who many bases with data###
#########################

snp_stats = function(wd = "/media/seb_1TB/transposon") {

setwd(wd)
species = c("deserticola","anomalus","paradoxus")

results = matrix(c(1:12),nrow =3, ncol = 4)

for(s in 1:3)
{
if(s == 1) {system("samtools mpileup -C50 -Iuf reference/HA412_trinity_noAltSplice_400bpmin.fa alignments/king141B.fq.sorted.bam | bcftools view -cg - | Rcode/vcfutils.pl vcf2fq > cns.fq");snp =nrow(read.table("results/snp_table_deserticola", stringsAsFactors = F))}
if(s == 2) {system("samtools mpileup -C50 -Iuf reference/HA412_trinity_noAltSplice_400bpmin.fa alignments/Sample_Ano1506.fq.sorted.bam | bcftools view -cg - | Rcode/vcfutils.pl vcf2fq > cns.fq");snp = nrow(read.table("results/snp_table_anomalus", stringsAsFactors = F))}
if(s == 3) {system("samtools mpileup -C50 -Iuf reference/HA412_trinity_noAltSplice_400bpmin.fa alignments/Sample_king1443.fq.sorted.bam | bcftools view -cg - | Rcode/vcfutils.pl vcf2fq > cns.fq");snp = nrow(read.table("results/snp_table_paradoxus", stringsAsFactors = F))}

cns = read.table("cns.fq", stringsAsFactors = F)

system(" grep -o 'n' cns.fq | wc -w >temp_n"); n = read.table("temp_n", stringsAsFactors = F)#number of "n"
system(" grep '@comp' cns.fq | wc >temp_seq"); seq = read.table("temp_seq", stringsAsFactors = F)#number of sequences
system("wc cns.fq >temp_char"); char = read.table("temp_char", stringsAsFactors = F) #number of characters

seq_length =  (char[3] - (char[2] + seq[2] + seq[3])) / 2 #total number of bp ### # char - length of name - length of quality identifier (2) - end of line character for whole file

results[s,1] = as.numeric(seq[1]) # 
results[s,2] = as.numeric((seq_length - n)) 
results[s,3] = as.numeric(snp)
results[s,4] = as.numeric((snp / (seq_length - n) * 100))
}
write.table(results,"results/snp_stats",row.names = species, col.names = T)
system("rm cns.fq temp*")
}



#########################
### RUN SNP CALLING FUNCTION ###
#########################

if(fun == "snp") {snp(list_ind= "reference/hybrid_species",user = "seb",species = args[1],cpu = as.numeric(args[2]),wd = "/media/seb_1TB/transposon")}
if(fun == "snp_stats") snp_stats()







