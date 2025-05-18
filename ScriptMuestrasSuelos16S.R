# Suelos 16S

library(dada2)
library (ggplot2)
library(plotly)
library(DECIPHER)
library(ShortRead)
library(Biostrings)
library(readr)
library("phyloseq")
library("DESeq2")
library("tibble")
library("dplyr")
library("tidyr")
library(openssl) 
library("vegan")
library("microbiome")
library("hrbrthemes") 
library("RColorBrewer")
library("data.table")
library(openxlsx)
library(ShortRead)
library(Biostrings)
library(readxl)


setwd("C:/Users/UsuarioPC/Desktop/Masterbioinf/TFM")

list.files(path = "./Muestras/Nat-Bloqu-05/", pattern = ".*suelos.*16S.*_001\\.fastq\\.gz$", 
           ignore.case = TRUE)

path<- './Muestras/Nat-Bloqu-05/'

fnFs <- sort(list.files(path, pattern="^_.*suelos.*16S.*_R1.*_001\\.fastq\\.gz$", full.names = TRUE)) 
#fnFs <- gsub("^\\./Muestras/Nat-Bloqu-05/", "", fnFs)
fnRs <- sort(list.files(path="./Muestras/Nat-Bloqu-05/", pattern=".*suelos.*16S.*_R2.*_001\\.fastq\\.gz$", full.names = TRUE))
#fnRs <- gsub("^\\./Muestras/Nat-Bloqu-05/", "", fnRs)

# Extraemos solo el nombre de las muestras
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2) 


# Cutadapt

cutadapt<-"./Muestras/Nat-Bloqu-05/cutadapt.exe"

FWD<- "MGGATTAGATACCCKGGT"
REV <- "ACGTCRTCCCCDCCTTCCT"
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

cut_dir <- file.path(".", "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

fnFs.cut <- file.path(cut_dir, basename(fnFs))
fnRs.cut <- file.path(cut_dir, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

names(fnFs.cut) <- sample.names
names(fnRs.cut) <- sample.names

minlen <- 150

cut_logs <- path.expand(file.path(cut_dir, paste0(sample.names, ".log")))

cutadapt_args <- c("-g", FWD, "-a", REV.RC, 
                   "-G", REV, "-A", FWD.RC,
                   "--match-read-wildcards",
                   "-n", 3, "--discard-untrimmed", "--minimum-length", minlen)

for (i in seq_along(fnFs)) {
  system2(cutadapt, 
          args = c(cutadapt_args,
                   "-o", fnFs.cut[i], "-p", fnRs.cut[i], 
                   fnFs[i], fnRs[i]),
          stdout = cut_logs[i])  
}

# Comprobamos
head(list.files(cut_dir))

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))



# Calidad Fordwards

forwplot<-ggplotly(plotQualityProfile(fnFs.cut[1:length(fnFs.cut)], aggregate=TRUE) + 
                     geom_hline(yintercept=c(15,25,35), 
                                color=c("red","blue","green"), 
                                size=0.5),
                   width =600) 
forwplot

# Calidad Reverse

revqplot<-ggplotly(plotQualityProfile(fnRs.cut[1:length(fnRs.cut)], aggregate=TRUE) + 
                     geom_hline(yintercept=c(15,25,35), 
                                color=c("red","blue","green"),
                                size=0.5),
                   width =600)
revqplot

# Filter and trimm 
filtFs <- file.path(".", "cutfiltered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(".", "cutfiltered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names #asignamos los nombres
out <- filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs, 
                     maxN=0, maxEE=c(2,5), trimLeft=c(17,17), trimRight=c(25,25), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
out<- cbind(out, perc.cons=round(out[, "reads.out"]/out[, "reads.in"]*100, digits=2))



# Aprender tasas de error 

errF <- learnErrors(filtFs,multithread=TRUE, nbases=130000000)
errR <- learnErrors(filtRs, multithread=TRUE, nbases=130000000)

     


#Ploteamos 

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)




# Inferencia de las muestras 

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

        


# Hacemos el merge
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
       

# Tabla de secuencias 
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 59 3295

# DISTRIBUTION

table(nchar(getSequences(seqtab))) 

# Eliminamos quimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab) # con lo que nos quedamos 0.96

# Resumen del filtrado

getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1:2], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track<- as.data.frame(track)
track$conservedperc<-round(with(track, nonchim/input*100),2)
print(track)

##############

# Asignamos taxonomía 

taxaSUasvs <- assignTaxonomy(seqtab.nochim, "./silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE, tryRC=TRUE)
taxaSUasvs <- assignTaxonomy(seqtab.nochim, "./silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE, tryRC=TRUE)


# Cargamos metadatos y lo modificamos

metadatossuelo <- read.xlsx("./Muestras/Nat-Bloqu-05/tax_otus_metadatos_suelo_12_09_24B.xlsx", sheet = 1)

# Por orden 
metadatossuelo <- metadatossuelo %>% 
  arrange(Column1)

samples.out <- rownames(seqtab.nochim) # Asignamos nombres
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

# Creamos variables para representar posteriormente

estres <- metadatossuelo$ELC
ciudad <- metadatossuelo$Poblacion
samdf <- data.frame(subject=subject, estres=estres, poblacion=ciudad)
samdf$estrespoblacion <- paste(estres,ciudad)
samdf$estres[samdf$estres == 1] <- "Estrés Alto"
samdf$estres[samdf$estres==2] <- "Estrés Medio"

samdf$estres[samdf$estres==3] <- "Estrés Bajo"
rownames(samdf) <- samples.out

#Creamos el objeto phyloseq

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                 sample_data(samdf), 
                 tax_table(taxaSUasvs))

# Eliminamos Mitocondria y cloroplastos

physeqbacsuelo16sASVs<-subset_taxa(ps, Order!="Mitochondria" & Order!="Chloroplast") 


# Evaluamos curvas de rarefracción

rarefactioncurve<-rarecurve(as.matrix(seqtab.nochim), step=100,las = 1, cex=0.75, tidy = TRUE) 
plotrarefaction <- ggplot(rarefactioncurve, aes(Sample, Species)) + 
  geom_line(aes(color=Site)) +
  xlab("Coverage") + ylab("ASV number") +
  scale_x_continuous(breaks = seq(0, max(rarefactioncurve$Sample), by = 70)) +  
  scale_y_continuous(breaks = seq(0, max(rarefactioncurve$Species), by = 2)) +  
  theme(legend.position = "none")
plotrarefaction <- ggplot(rarefactioncurve, aes(Sample, Species)) + 
  geom_line(aes(color=Site)) +
  xlab("Coverage") + ylab("ASV number") +
  theme(legend.position = "none") 
ggplotly(plotrarefaction)
rarefactioncurve$Sample <- substr(rarefactioncurve$Sample, 1, 10)


# Representación composición taxonómica

physeqbacsuelo16sASVs %>% transform(transform = "compositional") # transformamos

# Representamos 

physeq.agg<- aggregate_rare(physeqbacsuelo16sASVs %>% 
                              transform(transform="compositional"),
                            level="Genus", detection =0.05, prevalence=0.1)
taxa_names(physeq.agg)

# Modificamos los nombres para la representación 

taxa_names(physeq.agg)[taxa_names(physeq.agg) == "Other"] <- "Géneros <10% de abundancia"
taxa_names(physeq.agg)[taxa_names(physeq.agg) == "Unknown"] <- "NA"

getPalette <- colorRampPalette(brewer.pal(8, "Set2")) 
PhylaPalette <- getPalette(length(taxa(physeq.agg)))

# Representamos

taxcompplot<- plot_composition(physeq.agg, average_by="poblacion",
                               x_label="poblacion", group_by="estres")+
  scale_y_percent() +
  scale_fill_manual(values = PhylaPalette)
taxcompplot +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

taxcompplot$data$Tax <- factor(taxcompplot$data$Tax, levels = c(levels(taxcompplot$data$Tax)
                                                                [levels(taxcompplot$data$Tax)!="Unknown"], "Unknown"))
taxcompplot +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

##############################################
# Alfa diversidad 

# Definimos los índices

divIdx = c("Chao1", "Shannon", "Simpson")

# Representamos 

alphaplot_box <- plot_richness(physeqbacsuelo16sASVs, x = "estres", measures = divIdx, color = "estres", nrow = 1) + 
  geom_boxplot(aes(color = estres, fill = estres),
               alpha = 0.3,  
               size = 0.1,   
               show.legend = FALSE) + 
  geom_jitter(aes(color = estres), size = 0.6, width = 0.1, alpha = 0.7) +  
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Group", y = "Alpha Diversity Measure") + 
  stat_summary(aes(color = estres), fun = "median", geom = "point", shape = 18, size = 3)
alphaplot_box


# Análisis ANOVA

alpha_df <- estimate_richness(physeqbacsuelo16sASVs, measures = divIdx)
alpha_df$estres <- sample_data(physeqbacsuelo16sASVs)$estres

# Función para el ANOVA

anova_results <- lapply(divIdx, function(index) {
  formula <- as.formula(paste(index, "~ estres"))
  model <- aov(formula, data = alpha_df)
  summary(model)
})

names(anova_results) <- divIdx

# Ver resultados

anova_results



# Beta diversidad

# Método ordenación NMDS distancias Bray Curtis

ordinationNMDSbray <- ordinate(physeqbacsuelo16sASVs, method = "NMDS", distance = "bray")

# Representamos 

plotordinationbrayNMDS <- plot_ordination(physeqbacsuelo16sASVs, ordinationNMDSbray, color = "estres", shape = "estres") + 
  geom_point(size = 2) +  
  stat_ellipse(aes(group = estres), type = "t") + 
  theme_bw()   
plotordinationbrayNMDS




# Análisis de enriquecimiento

filtersamp<-genefilter_sample(physeqbacsuelo16sASVs, filterfun_sample(function(x) x > 1),
                              A=0.01*nsamples(physeqbacsuelo16sASVs))
physeqmayor <- prune_taxa(filtersamp, physeqbacsuelo16sASVs)
physeqmayor <- prune_samples(sample_sums(physeqmayor) > 0, physeqmayor)
phyestres <- phyloseq_to_deseq2(physeqmayor, ~ estres)
phyestres <- estimateSizeFactors(phyestres, type = "poscounts")

#Creamos el objeto DESeq

phyestres<- DESeq(phyestres, test="Wald", fitType="parametric")

# Alto vs bajo

altobajo <- results(phyestres, cooksCutoff = FALSE, contrast = c("estres", "Estrés Alto", "Estrés Bajo"))

# Significancia 

alpha <- 0.05
altobajo <- altobajo[which(altobajo$padj < alpha &(altobajo$log2FoldChange >=1 |
                                                      altobajo$log2FoldChange<=-1)), ]

altobajo = cbind(as(altobajo, "data.frame"), as(tax_table(physeqbacsuelo16sASVs)
                                                  [rownames(altobajo), ], "matrix"))

# Representamos a nivel de género 

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(altobajo$log2FoldChange, altobajo$Genus, function(x) max(x))
x = sort(x, TRUE)
altobajo$Genus = factor(as.character(altobajo$Genus), levels=names(x))
x = tapply(altobajo$log2FoldChange, altobajo$Family, function(x) max(x))
x = sort(x, TRUE)
altobajo$Family = factor(as.character(altobajo$Family), levels=names(x))
enrichplot <- ggplot(altobajo, aes(x = log2FoldChange, y = Genus, fill = Family)) +
  geom_col(width = 0.4, color = "black", size = 0.2) +  
  geom_vline(xintercept = 0, linetype = "solid", color = "gray40") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "#ff0000", size = 1, show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#00ff00", size = 1, show.legend = FALSE) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(),  
    panel.grid.major.x = element_line(color = "gray90", linetype = "dotted"),  
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.text.y = element_text(face = "bold.italic", size = 12),
    legend.text = element_text(face = "italic")
  ) +
  labs(
    title = "Estrés Alto vs Estrés Bajo Suelo 16S",
    x = "log2 Fold Change",
    y = "Género"
  )

enrichplot

##########################################################################################################