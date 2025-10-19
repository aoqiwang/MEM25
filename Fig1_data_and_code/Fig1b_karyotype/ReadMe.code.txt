For detailed information on the MEM25 chromosome and its annotations, we used the R package Rideogram (version 0.2.2) for plotting and presentation. 
It relies on the genomic gff file and the karyotype file containing the starting and ending positions of the filaments. 
It can be verified in R with the following code:

###code
library("RIdeogram")
human_karyotype <- read.table("karyotype.txt",sep ="\t",header =T,stringsAsFactors=F)
gene_density <- GFFex(input = "MEM25.all.gff3",karyotype = "karyotype.txt",feature = "snp",window=20000)
Random_RNAs_500 <- read.table("color_conf.txt",sep ="\t",header =T,stringsAsFactors=F)karyotype=human_karyotype
ideogram(karyotype, overlaid = NULL, label = NULL, label_type = NULL, synteny = NULL, colorset1, colorset2, width, Lx, Ly, output = "chromosome.svg")
