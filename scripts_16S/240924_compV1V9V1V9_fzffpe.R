setwd("D:/Professionnel/Thesescience/Microbiote/Analyses/Pancreas/Projects/comparatifs/v1v9fz_v1v9_ffpe")

#basics : https://bookdown.org/rdpeng/rprogdatascience/
{
library(ggplot2)
library(gridExtra)
library(phyloseq) #https://joey711.github.io/phyloseq/   , #https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html
library(dada2) #https://benjjneb.github.io/dada2/tutorial.html, https://benjjneb.github.io/dada2/bigdata_paired.html
library(DECIPHER) #http://www2.decipher.codes/index.html
library(phangorn) #https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html
library(vegan)
library(readxl)
library(DESeq2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ade4)
library("phyloseqGraphTest")
library("igraph")
library("ggnetwork")
library(stringr)
library(qutils)
library(edgeR)
library(corrplot)
library(decontam) #https://benjjneb.github.io/decontam/vignettes/decontam_intro.html; https://rpubs.com/microbiotic/833244
library(Hmisc)
library(ShortRead) 
library(survminer)
library(forestmodel)
library(survival)
library(writexl)
library(Maaslin2)
}

#process (https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html)
{
  pathF="./FastQ/Forward/4_paired"
  fnFs <- sort(list.files(pathF, pattern="R1.")) 
  fnFs <- file.path(pathF, fnFs) 
  
  pathR="./FastQ/Reverse/4_paired"
  fnRs <- sort(list.files(pathR, pattern="R2.")) 
  fnRs <- file.path(pathR, fnRs) 
  
  #plotQualityProfile(fnFs[c(3,14)])  #pas de read vide
  #plotQualityProfile(fnRs[c(3,14)])  #pas de read vide
  
  filtpathF <- file.path("./fastq/Forward/5_filtered_dada2")
  fastqFs <- sort(list.files(pathF, pattern="R1."))
  
  filtpathR <- file.path("./fastq/Reverse/5_filtered_dada2")
  fastqRs <- sort(list.files(pathR, pattern="R2."))
  
  out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                       rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                       compress=TRUE, multithread=F, verbose=T) 
  saveRDS(out,"R/rds/out.rds")
  out <- readRDS("R/rds/out.rds")

  filtFs <- list.files(filtpathF, pattern="R1.", full.names = TRUE)
  filtRs <- list.files(filtpathR, pattern="R2.", full.names = TRUE)
  sampleNames <- c(paste0("V663","_",1:23),paste0("V728","_",1:14)) 
  names(filtFs) <- sampleNames
  names(filtRs) <- sampleNames
  
  # Estimation of error rates 
  errF_V663=learnErrors(filtFs[1:23], multithread=T) 
  saveRDS(errF_V663,"R/rds/errF_V663.rds")
  errF_V663=readRDS("R/rds/errF_V663.rds")
  plotErrors(errF_V663, nominalQ=TRUE)
  
  errR_V663=learnErrors(filtRs[1:23], multithread=T) 
  saveRDS(errR_V663,"R/rds/errR_V663.rds")
  errR_V663=readRDS("R/rds/errR_V663.rds")
  plotErrors(errR_V663, nominalQ=TRUE)
  
  errF_V728=learnErrors(filtFs[24:37], multithread=T) 
  saveRDS(errF_V728,"R/rds/errF_V728.rds")
  errF_V728=readRDS("R/rds/errF_V728.rds")
  plotErrors(errF_V728, nominalQ=TRUE)
  
  errR_V728=learnErrors(filtRs[24:37], multithread=T) 
  saveRDS(errR_V728,"R/rds/errR_V728.rds")
  errR_V728=readRDS("R/rds/errR_V728.rds")
  plotErrors(errR_V728, nominalQ=TRUE)
  
  # Dereplication (remove paired sequences)
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  names(derepFs) <- sampleNames
  
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  names(derepRs) <- sampleNames
  
  # Infered composition of the sample (removing all sequencing errors to reveal members of community)
  dada_V663_F <- dada(derepFs[1:23], err=errF_V663, multithread=TRUE)
  saveRDS(dada_V663_F,"R/rds/dada_V663_F.rds")
  dada_V663_F=readRDS("R/rds/dada_V663_F.rds")
  
  dada_V663_R <- dada(derepRs[1:23], err=errR_V663, multithread=TRUE)
  saveRDS(dada_V663_R,"R/rds/dada_V663_R.rds")
  dada_V663_R=readRDS("R/rds/dada_V663_R.rds")
  
  dada_V728_F <- dada(derepFs[24:37], err=errF_V728, multithread=TRUE)
  saveRDS(dada_V728_F,"R/rds/dada_V728_F.rds")
  dada_V728_F=readRDS("R/rds/dada_V728_F.rds")
  
  dada_V728_R <- dada(derepRs[24:37], err=errR_V728, multithread=TRUE)
  saveRDS(dada_V728_R,"R/rds/dada_V728_R.rds")
  dada_V728_R=readRDS("R/rds/dada_V728_R.rds")
  
  # Remove chimeras
  mergers_V663 <- mergePairs(dada_V663_F, derepFs[1:23], dada_V663_R, derepRs[1:23])
  seqtab_V663 <- makeSequenceTable(mergers_V663)
  
  mergers_V728 <- mergePairs(dada_V728_F, derepFs[24:37], dada_V728_R, derepRs[24:37])
  seqtab_V728 <- makeSequenceTable(mergers_V728)
  
  seqtabAll <- mergeSequenceTables(seqtab_V663, seqtab_V728)
  
  seqtab <- removeBimeraDenovo(seqtabAll, method="consensus", multithread=TRUE, verbose=T)
  
  saveRDS(seqtab, "R/rds/seqtabNoC.rds")
  seqtab=readRDS("R/rds/seqtabNoC.rds")
  hist((nchar(getSequences(seqtab))))

  # Assign Taxonomy 
  fastaRef <- "R/rds/rdp_train_set_18.fa.gz" 
  taxTab <- assignTaxonomy(seqtab, refFasta = fastaRef, multithread=T, tryRC=T)
  saveRDS(taxTab, "R/rds/taxTab.rds")
  taxTab=readRDS("R/rds/taxTab.rds")
  1-length(which(is.na(taxTab[,"Genus"])))/nrow(taxTab) # % ASV assigned to the genus level
  
  fastaRefS <- "R/rds/rdp_species_assignment_18.fa" 
  taxTab_species <- addSpecies(taxTab, fastaRefS, verbose=TRUE, tryRC=T)
  saveRDS(taxTab_species, "R/rds/taxTab_species.rds")
  taxTab_species=readRDS("R/rds/taxTab_species.rds")
  1-length(which(is.na(taxTab_species[,"Species"])))/nrow(taxTab_species)
  ps_s <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(taxTab_species)) 
 
  #sanity check
  QC=do.call(rbind,lapply(1:length(sampleNames),function(x){
    fq=readFastq(fnFs[x])
    reads=sread(fq)
    quals=quality(fq)
    qscores=as(quals,"matrix")
    avgqscore=rowMeans(qscores, na.rm=T)
    
    data.frame(sample=sampleNames[x],
               length_reads=mean(reads@ranges@width),
               mean_score= mean(avgqscore)
    )
  }))
  
  blast=do.call(rbind,lapply(1:length(sampleNames),function(x){
    path <- "./FastQ/Forward/tag"
    tables <- sort(list.files(path, pattern=".txt")) 
    tables <- file.path(path, tables) 
    table <- read.delim(tables[x],header=T, skip=1)
    
    data.frame(sample=tables[x],
              total_reads=table$X.Reads_processed[1],
               reads_rdp=sum(table$X.One_hit_one_genome[2],table$X.Multiple_hits_one_genome[2]),
               reads_rdp_perc=sum(table$X.One_hit_one_genome.1[2],table$X.Multiple_hits_one_genome.1[2]),
               reads_human_perc=sum(table$X.One_hit_one_genome.1[1],table$X.Multiple_hits_one_genome.1[1]))
  }))
  rownames(blast)=sampleNames

  getN <- function(x) sum(getUniques(x))
  track_V663 <- cbind(out[1:23,], sapply(dada_V663_F, getN), sapply(dada_V663_R, getN), sapply(mergers_V663, getN), rowSums(seqtab[1:23,]))
  track_V728 <- cbind(out[24:37,], sapply(dada_V728_F, getN), sapply(dada_V728_R, getN), sapply(mergers_V728, getN), rowSums(seqtab[24:37,]))
  track <- rbind(track_V663, track_V728)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track)=sampleNames

  sanity=cbind(QC,track,blast)
  sanity$filt_perc=(sanity$filtered/sanity$input)*100
  sanity$ASV_perc=(sanity$denoisedF/sanity$filtered)*100
  sanity$merged_perc=(sanity$merged/sanity$denoisedF)*100
  sanity$chim_perc=(1-(sanity$nonchim/sanity$denoisedF))*100 

  saveRDS(sanity, "R/rds/sanity.rds")
  sanity=readRDS("R/rds/sanity.rds")
  
  # Combine clinical data
  clin=as.data.frame(read_excel("liste_cas.xlsx"))
  clin[,2:7] <- lapply(clin[,2:7] , factor)

  clin=cbind(sanity,clin)
  clin=clin[,-10]
  clin=clin[which(clin$Common=="T"),]
  
  a=ggplot(clin, aes(x=Type, y=length_reads, color=Type)) + 
    geom_boxplot()+ 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
    theme(legend.position = "none", axis.title.x=element_blank())+
    ggtitle("Mean Length Reads")
  
  b=ggplot(clin, aes(x=Type, y=filt_perc, color=Type)) + 
    geom_boxplot()+ 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
    theme(legend.position = "none", axis.title.x=element_blank())+
    ggtitle("% Unfiltered Reads")
  
  c=ggplot(clin, aes(x=Type, y=chim_perc, color=Type)) + 
    geom_boxplot()+ 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    theme(legend.position = "none", axis.title.x=element_blank())+
    ggtitle("% Chimeras")
  
  d=ggplot(clin, aes(x=Type, y=merged_perc, color=Type)) + 
    geom_boxplot()+ 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    theme(legend.position = "none", axis.title.x=element_blank())+
    ggtitle("% Merged")
  
  e=ggplot(clin, aes(x=Type, y=reads_rdp_perc, color=Type)) + 
    geom_boxplot()+ 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
    theme(legend.position = "none", axis.title.x=element_blank())+
    ggtitle("% Bacterial Reads")
  
  f=ggplot(clin, aes(x=Type, y=reads_human_perc, color=Type)) + 
    geom_boxplot()+ 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
    theme(legend.position = "none", axis.title.x=element_blank())+
    ggtitle("% Human Reads")
  
  ggarrange(a,c,e,f,
            ncol=2,nrow=2)
  
  # Using phyloseq 
  rownames(clin)=sampleNames
  ps_s <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), sample_data(clin),tax_table(taxTab)) 
  saveRDS(ps_s, "R/rds/ps_s.rds")
  ps_s=readRDS("R/rds/ps_s.rds")

  length(get_taxa_unique(ps_s, taxonomic.rank = "Genus")) 
  ps3_s = tax_glom(ps_s, "Genus", NArm = T)
  saveRDS(ps3_s, "R/rds/ps3_s.rds")

  ## control positif
  ps_ctrl_pos=subset_samples(ps3_s, Type_2 %in% c("cpost"))
  otu=as.data.frame(ps_ctrl_pos@otu_table)
  tax=as.data.frame(ps_ctrl_pos@tax_table)
  colnames(otu)=tax$Genus
  colnames(otu)[which(otu[1,]>0)]
  colnames(otu)[which(otu[2,]>0)]

  ps3ra = transform_sample_counts(ps_ctrl_pos, function(x){x / sum(x)})
  top5P.names = sort(tapply(taxa_sums(ps3ra), tax_table(ps3ra)[, "Genus"], sum), TRUE)[1:10]
  top5P = subset_taxa(ps3ra, Genus %in% names(top5P.names))
  plot_bar(top5P, x="Type",fill="Genus") +
    theme(legend.key.size = unit(0.3, 'cm'))
  
  ## controle negatif
  ps_ctrl_neg=subset_samples(ps3_s, Type_2 == "cneg")
  ps3ra = transform_sample_counts(ps_ctrl_neg, function(x){x / sum(x)})
  top5P.names = sort(tapply(taxa_sums(ps3ra), tax_table(ps3ra)[, "Genus"], sum), TRUE)[1:15]
  top5P = subset_taxa(ps3ra, Genus %in% names(top5P.names))
  plot_bar(top5P, x="Type",fill="Genus") + 
    facet_wrap(~Patient, scales="free_x", nrow=1)+
    theme(legend.key.size = unit(0.3, 'cm'))
  
  ## visualisation avant décontamination 
  ps3_s=subset_samples(ps3_s, ! Type %in% c("ctl_pos","ctl_pos_bis"))
  
  ps3_s@sam_data$Type[ps3_s@sam_data$Type=="ctl_neg_SC"]="ctl_neg"
  
  my_comparisons <- list( c("tumeur", "sain"),c("ctl_neg","tumeur"))
  plot_richness(ps3_s,x="Type", color="Type", measures=c("Shannon","Chao1"))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7)+
    stat_compare_means(comparisons = my_comparisons)+
    xlab("Tissue") +
    theme(legend.position = "none")
  
  my_comparisons <- list( c("tumeur", "sain"),c("ctl_neg","tumeur"))
  plot_richness(ps3_s,x="Type_bis", color="Type", measures=c("Shannon","Chao1"))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7)+
    stat_compare_means(comparisons = my_comparisons)+
    xlab("Tissue") 
  
  ## remove contaminants by prevalence  (attention, avoir vision critique pour ne pas enlever contaminants et bactéries présentes)
  ps_dec <- ps3_s
  sample_data(ps_dec)$is.neg <- sample_data(ps_dec)$Type_bis == "CTL"
  ps.pa <- transform_sample_counts(ps_dec, function(abund) 1*(abund>0))
  ps.pa.neg <- prune_samples(sample_data(ps.pa)$Type_bis == "CTL", ps.pa)
  ps.pa.pos <- prune_samples(sample_data(ps.pa)$Type_bis != "Negative control", ps.pa)
  
  ### method stat
  contamdf.prev <- isContaminant(ps_dec, method="prevalence", neg="is.neg", threshold=0.20)
  table(contamdf.prev$contaminant)
  df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=contamdf.prev$contaminant)
  tax=data.frame(ps_dec@tax_table)
  rownames(df.pa)=tax$Genus
  ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  contam=rownames(contamdf.prev)[which(contamdf.prev$contaminant==T)]  
  tax=tax[contam,]
  ps3_s <- prune_taxa(!contamdf.prev$contaminant, ps3_s) 
  
  ### method manuelle 
  df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg), 
                      contaminant=(taxa_sums(ps.pa.neg)>0)&(taxa_sums(ps.pa.pos)<4)) 
  ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  contam=rownames(df.pa[which(df.pa$contaminant==T),])
  tax=data.frame(ps_dec@tax_table)
  tax=tax[contam,]
  tax$Genus
  ps3_s <- prune_taxa(!(df.pa$contaminant==T), ps3_s)
  
  ### % contaminants connus
  df.contam=read_excel("D:/Professionnel/Thesescience/Microbiote/Analyses/Contaminants_MH.xlsx")
  otu=as.data.frame(ps3_s@otu_table)
  tax=as.data.frame(ps3_s@tax_table)
  colnames(otu)=tax$Genus
  
  fun=function(x){
    length( intersect (colnames(otu)[which(otu[x,]>0)],df.contam$Genre))/length(colnames(otu)[which(otu[x,]>0)])
  }
  perc_cont=unlist(lapply(1:37,fun))
  df_perc_cont=data.frame(perc_cont=perc_cont, 
                          type=ps3_s@sam_data$Type)
  df_perc_cont=df_perc_cont[which(df_perc_cont$type=="FZ"|df_perc_cont$type=="FFPE"),]
  ggplot(df_perc_cont, aes(x=type, y=perc_cont, color=type))+
    geom_boxplot()+
    geom_jitter(shape=16, position=position_jitter(0.2))+
    stat_compare_means()
  
  ps3_s_decont = subset_taxa(ps3_s, !Genus %in% df.contam$Genre)
  my_comparisons <- list( c("tumeur", "sain"),c("ctl_neg","tumeur"))
  plot_richness(ps3_s_decont,x="Type", color="Type", measures=c("Shannon","Chao1"))+
    geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7)+
    stat_compare_means(comparisons = my_comparisons)+
    xlab("Tissue") +
    theme(legend.position = "none")
  
  ## Compute prevalence 
  ps_samples=subset_samples(ps3_s, Type_2 == "tumor")
  prevdf = apply(X = otu_table(ps_samples),
                 MARGIN = ifelse(taxa_are_rows(ps_samples), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  prevdf = data.frame(Prevalence = prevdf, 
                      TotalAbundance = taxa_sums(ps_samples), 
                      tax_table(ps_samples)) 
  prev=plyr::ddply(prevdf, "Genus", function(df1){cbind(mean(df1$TotalAbundance),sum(df1$Prevalence))}) #Compute the total and average prevalences of the features in each phylum
  colnames(prev)=c("Genus","MeanTotalAbundance","SamplePrevalence")
  filt=prev[which(prev$`2`<1),1]  ##ne pas filtrer si bactéries super basse
  
  ps3_s = subset_taxa(ps3_s, !Genus %in% filt) 
  saveRDS(ps3_s, "R/rds/ps3_s.rds")
  ps3_s=readRDS("R/rds/ps3_s.rds")
  
  # Multiple alignement using DECIPHER
  seqs=colnames(ps3_s@otu_table)
  names(seqs) <- seqs
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE, processors=20)
  saveRDS(alignment, "R/rds/alignment.rds")
  alignment=readRDS("R/rds/alignment.rds")
  
  # Construction of phylogenetic tree 
  phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phangAlign)
  saveRDS(dm, "R/rds/dm.rds")
  dm=readRDS("R/rds/dm.rds")
  
  treeNJ <- NJ(dm) 
  saveRDS(treeNJ, "R/rds/treeNJ.rds")
  treeNJ=readRDS("R/rds/treeNJ.rds")
  
  fit = pml(treeNJ, data=phangAlign) # Fit a  GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree
  fit= update(fit, k=4, inv=0.2)
  
  fit <- optim.pml(fit,optInv=TRUE, optGamma=TRUE,model="GTR") 
  saveRDS(fit, "R/rds/fit.rds")
  fit=readRDS("R/rds/fit.rds")
  ps3_s=merge_phyloseq(ps3_s,phy_tree(fit$tree))
  saveRDS(ps3_s, "R/rds/ps3_s_tree.rds")
}

#analysis
ps3=readRDS("R/rds/ps3_s_tree.rds")
#rarecurve(otu_table(ps3), step=50, cex=0.5, label=F, xlim=c(0,3000)) 
ps3@sam_data$charge_bact=rowSums(ps3@otu_table)
ps3@sam_data$charge_bact_perc=rowSums(ps3@otu_table)/ps3@sam_data$total_reads
ps3@sam_data$nb_taxa=rowSums(ps3@otu_table != 0)

clin_ps=data.frame(ps3@sam_data)

## Relative abundance
ps_pair <- subset_samples(ps3, Common == "T")
ps3ra = transform_sample_counts(ps_pair, function(x){x / sum(x)})
top5P.names = sort(tapply(taxa_sums(ps3ra), tax_table(ps3ra)[, "Phylum"], sum), TRUE)[1:15]
top5P = subset_taxa(ps3ra, Phylum %in% names(top5P.names))
plot_bar(top5P, x="Type",fill="Phylum") + 
  facet_wrap(~Patient, scales="free_x", nrow=1)+
  theme(legend.key.size = unit(0.3, 'cm'))

top5P.names = sort(tapply(taxa_sums(ps3ra), tax_table(ps3ra)[, "Genus"], sum), TRUE)[1:10]
top5P = subset_taxa(ps3ra, Genus %in% names(top5P.names))
plot_bar(top5P, x="Type",fill="Genus") + 
  facet_wrap(~Patient, scales="free_x", nrow=1)+
  theme(legend.key.size = unit(0.3, 'cm'))
  
ps4= merge_samples(top5P, "Type_2")
ps4ra <- transform_sample_counts(ps4, function(x) x / sum(x))
plot_bar(ps4ra, fill="Phylum")+ 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")+
  theme(legend.key.size = unit(0.3, 'cm'))

###heatmap genres les plus fréquents
ps3_sub <- prune_taxa(names(sort(taxa_sums(ps_samples),TRUE)[1:20]), ps_samples)
plot_heatmap(ps3_sub, sample.label="Type",taxa.label="Genus",low="#66CCFF", high="#000033", na.value="white")

ps_pair <- subset_samples(ps3, Common == "T")
ps3ra = transform_sample_counts(ps_pair, function(x){x / sum(x)})
top5P.names = sort(tapply(taxa_sums(ps3ra), tax_table(ps3ra)[, "Genus"], sum), TRUE)[1:15]
top5P = subset_taxa(ps3ra, Genus %in% names(top5P.names))
plot_bar(top5P, x="Type",fill="Genus") + 
  facet_wrap(~Patient, scales="free_x", nrow=1)+
  theme(legend.key.size = unit(0.3, 'cm'))

ps4= merge_samples(top5P, "Biopsie")
ps4ra <- transform_sample_counts(ps4, function(x) x / sum(x))
plot_bar(ps4ra, fill="Genus")+ geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

# alpha div comparison
ps_pair <- subset_samples(ps3, Common == "T")
alpha_div=estimate_richness(ps_pair)
alpha_div$Kit=ps_pair@sam_data$Type_3
alpha_div$Patient=ps_pair@sam_data$Patient
alpha_div$Common=ps_pair@sam_data$Common
alpha_div$Type=ps_pair@sam_data$Type

a=ggplot(alpha_div, aes(x=Type, y=Chao1, color=Kit)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  theme(legend.position = "none", axis.title.x=element_blank())+
  stat_compare_means()+
  ggtitle("alpha-diversity")
b=ggplot(alpha_div, aes(x=Type, y=Shannon, color=Kit)) + 
  geom_boxplot()+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  theme(legend.position = "none", axis.title.x=element_blank())+
  stat_compare_means()+
  ggtitle("alpha-diversity")

ggarrange(a,b)

alpha_div=alpha_div[which(alpha_div$Common=="T"),]
alpha_FFPE=alpha_div[which(alpha_div$Kit=="FFPE"),]
alpha_FFPE=alpha_FFPE[order(alpha_FFPE$Patient),]
alpha_FZ=alpha_div[which(alpha_div$Kit=="FZ"),]
alpha_FZ=alpha_FZ[order(alpha_FZ$Patient),]

alpha_ggscat_chao1=data.frame(FZ=alpha_FZ$Chao1,FFPE=alpha_FFPE$Chao1)
alpha_ggscat_shannon=data.frame(FZ=alpha_FZ$Shannon,FFPE=alpha_FFPE$Shannon)

a=ggscatter(alpha_ggscat_chao1, x = "FZ", y = "FFPE",
          add = "reg.line",                               
          conf.int = TRUE)+
  stat_cor(method = "spearman")+
  ggtitle("Chao1 index")
b=ggscatter(alpha_ggscat_shannon, x = "FZ", y = "FFPE",
            add = "reg.line",                               
            conf.int = TRUE)+
  stat_cor(method = "spearman") +
  ggtitle("Shannon index")

ggarrange(a,b)

#comparer otu entre echantillons
##quali
ps_pair <- subset_samples(ps3, Common == "T")
ps3ra = transform_sample_counts(ps_pair, function(x){x / sum(x)})
top5P.names = sort(tapply(taxa_sums(ps3ra), tax_table(ps3ra)[, "Genus"], sum), TRUE)[1:10]
top5P = subset_taxa(ps3ra, Genus %in% names(top5P.names))

otu=as.data.frame(ps_s@otu_table)
tax=as.data.frame(ps_s@tax_table)
colnames(otu)=tax$Species
rownames(otu)[which(otu$somerae>0)]
rownames(otu)
otu$Kit=top5P@sam_data$Type_3
otu$Patient=top5P@sam_data$Patient
otu$Common=top5P@sam_data$Common

otu_V3V4=otu[which(otu$Kit=="V3V4"),]
otu_V3V4=otu_V3V4[order(otu_V3V4$Patient),]
otu_V1V9=otu[which(otu$Kit=="V1V9"),]
otu_V1V9=otu_V1V9[order(otu_V1V9$Patient),]

fun=function(x){
  length(intersect(
    colnames(otu_V3V4)[which(otu_V3V4[1,]>0)],
    colnames(otu_V1V9)[which(otu_V1V9[x,]>0)]
  )
  )
}
test=unlist(lapply(1:12,fun))
names(test)=otu_V1V9$Patient
test

##quanti 
###bacteria
ps_pair <- subset_samples(ps3, Common == "T")
ps3ra = transform_sample_counts(ps_pair, function(x){x / sum(x)})
top5P.names = sort(tapply(taxa_sums(ps3ra), tax_table(ps3ra)[, "Genus"], sum), TRUE)[1:110]
top5P = subset_taxa(ps3ra, Genus %in% names(top5P.names))

otu=as.data.frame(top5P@otu_table)
tax=as.data.frame(top5P@tax_table)
colnames(otu)=tax$Genus
otu$Kit=top5P@sam_data$Type_3
otu$Patient=top5P@sam_data$Patient
otu$Common=top5P@sam_data$Common

otu_V3V4=otu[which(otu$Kit=="V3V4"),]
otu_V3V4=otu_V3V4[order(otu_V3V4$Patient),]
otu_V1V9=otu[which(otu$Kit=="V1V9"),]
otu_V1V9=otu_V1V9[order(otu_V1V9$Patient),]

fun=function(x){
  cor.test(as.numeric(otu_V3V4[,x]),as.numeric(otu_V1V9[,x]), 
         method="spearman")$p.value
}
test_pval=unlist(lapply(1:110,fun))
hist(test_pval)
length(which(test_pval<0.05))
tax$Genus[which(test_pval<0.05)]

###samples
fun=function(x){
  cor.test(as.numeric(otu_V3V4[x,]),as.numeric(otu_V1V9[x,]), 
           method="spearman")$p.value
}
test_pval=unlist(lapply(1:12,fun))
length(which(test_pval<0.0001))

otu_ord=t(rbind(otu_V3V4,otu_V1V9))
otu_ord=otu_ord[-c(111:113),]
otu_ord <- as.data.frame(apply(otu_ord, 2, as.numeric)) 
a=ggscatter(otu_ord,x = colnames(otu_ord)[1], y = colnames(otu_ord)[13],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
b=ggscatter(otu_ord,x = colnames(otu_ord)[2], y = colnames(otu_ord)[14],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
c=ggscatter(otu_ord,x = colnames(otu_ord)[3], y = colnames(otu_ord)[15],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
d=ggscatter(otu_ord,x = colnames(otu_ord)[4], y = colnames(otu_ord)[16],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
e=ggscatter(otu_ord,x = colnames(otu_ord)[5], y = colnames(otu_ord)[17],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
f=ggscatter(otu_ord,x = colnames(otu_ord)[6], y = colnames(otu_ord)[18],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
g=ggscatter(otu_ord,x = colnames(otu_ord)[7], y = colnames(otu_ord)[19],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
h=ggscatter(otu_ord,x = colnames(otu_ord)[8], y = colnames(otu_ord)[20],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
i=ggscatter(otu_ord,x = colnames(otu_ord)[9], y = colnames(otu_ord)[21],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
j=ggscatter(otu_ord,x = colnames(otu_ord)[10], y = colnames(otu_ord)[22],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
k=ggscatter(otu_ord,x = colnames(otu_ord)[11], y = colnames(otu_ord)[23],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 
l=ggscatter(otu_ord,x = colnames(otu_ord)[12], y = colnames(otu_ord)[24],
            add = "reg.line",                               
            conf.int = TRUE,                                 
            add.params = list(color = "blue", fill = "lightgray"))+
  stat_cor(method = "spearman") 

ggarrange(a,b,c,d,e,f,g,h,i,j,k,l,
          nrow=4,
          ncol=3)

## PCA https://joey711.github.io/phyloseq/plot_ordination-examples.html
### Bray-Curtis
ps_pair <- subset_samples(ps3, Common == "T")
out.pcoa.log <- ordinate(ps_pair,  method = "MDS", distance = "bray")
evals_bc <- out.pcoa.log$values[,1]
plot_ordination(ps_pair, out.pcoa.log, color="Patient", shape="Type", axes=c(1:2)) +
  coord_fixed(sqrt(evals_bc[2] / evals_bc[1]))+
  geom_point(size=3)
plot_ordination(ps_pair, out.pcoa.log, color="Patient", shape="Type", axes=c(3:4)) +
  coord_fixed(sqrt(evals_bc[4] / evals_bc[3]))+
  geom_point(size=3)

### Double principal coordinates analysis (DPCoA)
out.dpcoa.log <- ordinate(ps_pair, method = "DPCoA")
evals_dpcoa <- out.dpcoa.log$eig
plot_ordination(ps_pair, out.dpcoa.log, color="Patient", shape="Run")+
  coord_fixed(sqrt(evals_dpcoa[2] / evals_dpcoa[1]))+
  geom_point(size=3, alpha=0.75) 

### Unifrac weighted
out.wuf.log <- ordinate(ps_pair, method = "PCoA", distance ="wunifrac")
evals_wu <- out.wuf.log$values$Eigenvalues
plot_ordination(ps3, out.wuf.log,color="Patient", shape="Run",axes=c(1:2))+
  coord_fixed(sqrt(evals_wu[2] / evals_wu[1]))+
  geom_point(size=3, alpha=0.75) 
plot_ordination(ps3, out.wuf.log,color="Patient", shape="Run",axes=c(3:4))+
  coord_fixed(sqrt(evals_wu[4] / evals_wu[3]))+
  geom_point(size=3, alpha=0.75) 

out.wuf.log <- ordinate(ps3_ds, method = "PCoA", distance ="wunifrac")
pca_res=plot_ordination(ps3_ds, out.wuf.log,axes=c(1:5),justDF=T)

##comparer 
ps3=readRDS("R/rds/ps3_s_tree.rds")
ps_pair <- subset_samples(ps3, Common == "T")
otu=as.data.frame(ps_pair@otu_table)
tax=as.data.frame(ps_pair@tax_table)
colnames(otu)=tax$Genus

## union
genus_fz=names(which(colSums(otu[1,])>0))
genus_ffpe=names(which(colSums(otu[12,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
a=ggscatter(otu_comm, x = colnames(otu_comm)[1], y = colnames(otu_comm)[12],
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"),
          xlim=c(0,0.05), ylim=c(0,0.05)
          )+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[2,])>0))
genus_ffpe=names(which(colSums(otu[13,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
b=ggscatter(otu_comm, x = colnames(otu_comm)[2], y = colnames(otu_comm)[13],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[3,])>0))
genus_ffpe=names(which(colSums(otu[14,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
c=ggscatter(otu_comm, x = colnames(otu_comm)[3], y = colnames(otu_comm)[14],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[4,])>0))
genus_ffpe=names(which(colSums(otu[15,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
d=ggscatter(otu_comm, x = colnames(otu_comm)[4], y = colnames(otu_comm)[15],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[5,])>0))
genus_ffpe=names(which(colSums(otu[16,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
e=ggscatter(otu_comm, x = colnames(otu_comm)[5], y = colnames(otu_comm)[16],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[6,])>0))
genus_ffpe=names(which(colSums(otu[17,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
f=ggscatter(otu_comm, x = colnames(otu_comm)[6], y = colnames(otu_comm)[17],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)


genus_fz=names(which(colSums(otu[7,])>0))
genus_ffpe=names(which(colSums(otu[18,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
g=ggscatter(otu_comm, x = colnames(otu_comm)[7], y = colnames(otu_comm)[18],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[8,])>0))
genus_ffpe=names(which(colSums(otu[19,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
h=ggscatter(otu_comm, x = colnames(otu_comm)[8], y = colnames(otu_comm)[19],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[9,])>0))
genus_ffpe=names(which(colSums(otu[20,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
i=ggscatter(otu_comm, x = colnames(otu_comm)[9], y = colnames(otu_comm)[20],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[10,])>0))
genus_ffpe=names(which(colSums(otu[21,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
j=ggscatter(otu_comm, x = colnames(otu_comm)[10], y = colnames(otu_comm)[21],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[11,])>0))
genus_ffpe=names(which(colSums(otu[22,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
k=ggscatter(otu_comm, x = colnames(otu_comm)[11], y = colnames(otu_comm)[22],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.05), ylim=c(0,0.05)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

ggarrange(a,b,c,d,e,f,g,h,i,j,k,nrow=4, ncol=3)


## intersect
genus_fz=names(which(colSums(otu[1,])>0))
genus_ffpe=names(which(colSums(otu[12,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
a=ggscatter(otu_comm, x = colnames(otu_comm)[1], y = colnames(otu_comm)[12],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[2,])>0))
genus_ffpe=names(which(colSums(otu[13,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
b=ggscatter(otu_comm, x = colnames(otu_comm)[2], y = colnames(otu_comm)[13],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[3,])>0))
genus_ffpe=names(which(colSums(otu[14,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
c=ggscatter(otu_comm, x = colnames(otu_comm)[3], y = colnames(otu_comm)[14],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[4,])>0))
genus_ffpe=names(which(colSums(otu[15,])>0))
comm_bact=union(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
d=ggscatter(otu_comm, x = colnames(otu_comm)[4], y = colnames(otu_comm)[15],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[5,])>0))
genus_ffpe=names(which(colSums(otu[16,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
e=ggscatter(otu_comm, x = colnames(otu_comm)[5], y = colnames(otu_comm)[16],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[6,])>0))
genus_ffpe=names(which(colSums(otu[17,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
f=ggscatter(otu_comm, x = colnames(otu_comm)[6], y = colnames(otu_comm)[17],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)


genus_fz=names(which(colSums(otu[7,])>0))
genus_ffpe=names(which(colSums(otu[18,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
g=ggscatter(otu_comm, x = colnames(otu_comm)[7], y = colnames(otu_comm)[18],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[8,])>0))
genus_ffpe=names(which(colSums(otu[19,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
h=ggscatter(otu_comm, x = colnames(otu_comm)[8], y = colnames(otu_comm)[19],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[9,])>0))
genus_ffpe=names(which(colSums(otu[20,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
i=ggscatter(otu_comm, x = colnames(otu_comm)[9], y = colnames(otu_comm)[20],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[10,])>0))
genus_ffpe=names(which(colSums(otu[21,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
j=ggscatter(otu_comm, x = colnames(otu_comm)[10], y = colnames(otu_comm)[21],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

genus_fz=names(which(colSums(otu[11,])>0))
genus_ffpe=names(which(colSums(otu[22,])>0))
comm_bact=intersect(genus_fz,genus_ffpe)                  
ps_comm = subset_taxa(ps_pair, Genus %in% comm_bact)
ps_commra = transform_sample_counts(ps_comm, function(x){x / sum(x)})
otu_comm=as.data.frame(t(ps_commra@otu_table))
k=ggscatter(otu_comm, x = colnames(otu_comm)[11], y = colnames(otu_comm)[22],
            add = "reg.line",                                 # Add regression line
            conf.int = TRUE,                                  # Add confidence interval
            add.params = list(color = "blue",
                              fill = "lightgray"),
            xlim=c(0,0.25), ylim=c(0,0.25)
)+
  stat_cor(method = "spearman",label.x = 0.03, label.y = 0.05)

ggarrange(a,b,c,d,e,f,g,h,i,j,k,nrow=4, ncol=3)
# tester package adonis et package mantel https://www.youtube.com/watch?v=oLf0EpMJ4yA