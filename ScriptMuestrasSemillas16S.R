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

setwd("C:/Users/UsuarioPC/Desktop/Masterbioinf/TFM/")

list.files(path = "./Muestras/Nat-Bloqu-05/", pattern = ".*semillas.*16S.*_001\\.fastq\\.gz$", 
           ignore.case = TRUE)

path<- './Muestras/Nat-Bloqu-05/'

fnFsSe <- sort(list.files(path, pattern="^_.*semillas.*16S.*_R1.*_001\\.fastq\\.gz$", full.names = TRUE)) 
#fnFs <- gsub("^\\./Muestras/Nat-Bloqu-05/", "", fnFs)
fnRsSe <- sort(list.files(path="./Muestras/Nat-Bloqu-05/", pattern=".*semillas.*16S.*_R2.*_001\\.fastq\\.gz$", full.names = TRUE))
#fnRs <- gsub("^\\./Muestras/Nat-Bloqu-05/", "", fnRs)


# Extraemos solo el nombre de las muestras de semillas

sample.namesSe <- sapply(strsplit(basename(fnFsSe), "_"), `[`, 2) 

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
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFsSe[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRsSe[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFsSe[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRsSe[[1]]))

cut_dirSe <- file.path(".", "cutadapt")
if (!dir.exists(cut_dirSe)) dir.create(cut_dirSe)

fnFsSe.cut <- file.path(cut_dirSe, basename(fnFsSe))
fnRsSe.cut <- file.path(cut_dirSe, basename(fnRsSe))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

names(fnFsSe.cut) <- sample.namesSe
names(fnRsSe.cut) <- sample.namesSe

minlen <- 150

cut_logs <- path.expand(file.path(cut_dirSe, paste0(sample.namesSe, ".log")))

cutadapt_args <- c("-g", FWD, "-a", REV.RC, 
                   "-G", REV, "-A", FWD.RC,
                   "--match-read-wildcards",
                   "-n", 3, "--discard-untrimmed", "--minimum-length", minlen)

for (i in seq_along(fnFsSe)) {
  system2(cutadapt, 
          args = c(cutadapt_args,
                   "-o", fnFsSe.cut[i], "-p", fnRsSe.cut[i], 
                   fnFsSe[i], fnRsSe[i]),
          stdout = cut_logs[i])  
}

# Comprobamos
head(list.files(cut_dirSe))

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFsSe.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRsSe.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFsSe.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRsSe.cut[[1]]))




# Calidad Fordwards

forwplotSe<-ggplotly(plotQualityProfile(fnFsSe.cut[1:length(fnFsSe.cut)], aggregate=TRUE) + 
                     geom_hline(yintercept=c(15,25,35), 
                                color=c("red","blue","green"), 
                                size=0.5),
                   width =600) 
forwplotSe

# Calidad Reverse
revqplotSe<-ggplotly(plotQualityProfile(fnRsSe.cut[1:length(fnRsSe.cut)], aggregate=TRUE) + 
                     geom_hline(yintercept=c(15,25,35), 
                                color=c("red","blue","green"),
                                size=0.5),
                   width =600)
revqplotSe

# Filter and trimm 

filtFsSe <- file.path(path, "filtered", paste0(sample.namesSe, "_F_filt.fastq.gz"))
filtRsSe <- file.path(path, "filtered", paste0(sample.namesSe, "_R_filt.fastq.gz"))
save(filtFsSe, file = "filtFsSe.RData")
names(filtFsSe) <- sample.namesSe
names(filtRsSe) <- sample.namesSe
outSe <- filterAndTrim(fnFsSe.cut, filtFsSe, fnRsSe.cut, filtRsSe, 
                     maxN=0, maxEE=c(2,5), trimLeft=c(17,17), trimRight=c(25,25), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)


# Aprender tasas de error 

errFSe <- learnErrors(filtFsSe, multithread=TRUE)
errRSe <- learnErrors(filtRsSe, multithread=TRUE)





# Inferencia de las muestras 

dadaFsSe <- dada(filtFsSe, err=errFSe, multithread=TRUE)
dadaRsSe <- dada(filtRsSe, err=errRSe, multithread=TRUE)


# Hacemos el merge

mergersSe <- mergePairs(dadaFsSe, filtFsSe, dadaRsSe, filtRsSe, verbose=TRUE)

# Tabla de secuencias 

seqtabSe <- makeSequenceTable(mergersSe)
dim(seqtabSe) # 66 12724

# DISTRIBUTION

table(nchar(getSequences(seqtabSe))) # DISTRIBUTION

# Eliminamos quimeras

seqtab.nochimSe <- removeBimeraDenovo(seqtabSe, method="consensus", multithread=TRUE, verbose=TRUE)

sum(seqtab.nochimSe)/sum(seqtabSe) # con lo que nos quedamos 0.66

# Resumen filtrado

getN <- function(x) sum(getUniques(x))
trackSe <- cbind(outSe, sapply(dadaFsSe, getN), sapply(dadaRsSe, getN), sapply(mergersSe, getN), rowSums(seqtab.nochimSe))
colnames(trackSe) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(trackSe) <- sample.namesSe
head(trackSe)
trackSe

# Asignamos taxonomía

taxaSeasvs <- assignTaxonomy(seqtab.nochimSe, "./silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE, tryRC=TRUE)
taxaSeasvs <- assignTaxonomy(seqtab.nochimSe, "./silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE, tryRC=TRUE)


# Cargamos metadatos y modificamos

metadatossemillas <- read.xlsx("./Muestras/Nat-Bloqu-05/tax_otus_metadatos_semillas_12_09_24.xlsx", sheet = 1)

metadatossemillas <- metadatossemillas %>% 
  arrange(Samples)
samples.outse <- rownames(seqtab.nochimSe)
subjectSe <- sapply(strsplit(samples.outse, "D"), `[`, 1)
estresSe <- metadatossemillas$ELC
ciudadSe <- metadatossemillas$Poblacion
samdfse <- data.frame(subjectSe=subjectSe, estresSe=estresSe, poblacionSe=ciudadSe)
samdfse$estrespoblacionse <- paste(estresSe,ciudadSe)
samdfse$estresSe[samdfse$estresSe == 1] <- "Estrés Alto"
samdfse$estresSe[samdfse$estresSe==2] <- "Estrés Medio"

samdfse$estresSe[samdfse$estresSe==3] <- "Estrés Bajo"
samdfse$estresSe[samdfse$estresSe=="Cultivar "] <- "C"
samdfse$estresSe[samdfse$estresSe=="Cultivar"] <- "C"
samdfse$estresSe[samdfse$estresSe=="Cultivar no fungicida"] <- "C"


rownames(samdfse) <- samples.outse

# Para que coincidan sample.names y no nos de error a la hora de hacer el objeto phyloseq

rownames(samdfse) <- sprintf("%03d", as.integer(rownames(samdfse))) 


# Creamos objeto phyloseq

psSe <- phyloseq(otu_table(seqtab.nochimSe, taxa_are_rows=FALSE), 
               sample_data(samdfse), 
               tax_table(taxaSeasvs))

# Eliminamos mitocondria y cloroplastos

physeqbacsemillas16s<-subset_taxa(psSe, Order!="Mitochondria" & Order!="Chloroplast") 




# Curvas de rarefracción

rarefactioncurve<-rarecurve(as.matrix(seqtab.nochimSe), step=100,las = 1, cex=0.75, tidy = TRUE) 
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


physeq.aggSe<- aggregate_rare(physeqbacsemillas16s %>% 
                              transform(transform="compositional"),
                            level="Genus", detection =0.05, prevalence=0.1)

physeq.filt <- subset_samples(physeq.aggSe)

# Cambiamos nombres para la representación

taxa_names(physeq.filt)[taxa_names(physeq.filt) == "Other"] <- "Géneros <10% de abundancia"
taxa_names(physeq.filt)[taxa_names(physeq.filt) == "Unknown"] <- "NA"

# Colores

getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
PhylaPalette <- getPalette(length(taxa(physeq.filt)))

# Representamos

taxcompplot<- plot_composition(physeq.filt, average_by="poblacionSe",
                               x_label="poblacionSe", group_by="estresSe")+
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

alphaplot_box <- plot_richness(physeqbacsemillas16s, x = "estresSe", measures = divIdx, color = "estresSe", nrow = 1) + 
  geom_boxplot(aes(color = estresSe, fill = estresSe), 
               alpha = 0.3,  
               size = 0.3,   
               show.legend = FALSE) +  
  geom_jitter(aes(color = estresSe), size = 0.6, width = 0.1, alpha = 0.7) +  
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Group", y = "Alpha Diversity Measure") + 
  stat_summary(aes(color = estresSe), fun = "median", geom = "point", shape = 18, size = 3)  
alphaplot_box


# Análisis ANOVA

# Filtramos

physeq.filt <- subset_samples(physeqbacsemillas16s, 
                              !poblacionSe %in% c("OSK", "ZAPN"))

alpha_df <- estimate_richness(physeq.filt, measures = divIdx)
alpha_df$estres <- sample_data(physeq.filt)$estresSe

# Función para el ANOVA

anova_results <- lapply(divIdx, function(index) {
  formula <- as.formula(paste(index, "~ estres"))
  model <- aov(formula, data = alpha_df)
  summary(model)
})
names(anova_results) <- divIdx
anova_results

# Beta diversidad

# Método ordenación NMDS distancias Bray Curtis

ordinationNMDSbray<-ordinate(physeq.filt, method="NMDS", distance="bray")


plotordinationbray<-plot_ordination(physeq.filt, ordinationNMDSbray, color="estresSe",
                                    shape = "estresSe") +
  geom_point(size=2) + stat_ellipse(aes(group = estresSe), type = "t") +
  theme_bw()  
plotordinationbray

#####################################################################################

# Análisis de enriquecimiento

filtersampSe<-genefilter_sample(physeq.filt, filterfun_sample(function(x) x > 1),
                              A=0.01*nsamples(physeq.filt))
physeqmayorSe <- prune_taxa(filtersampSe, physeq.filt)
physeqmayorSe <- prune_samples(sample_sums(physeqmayorSe) > 0, physeqmayorSe)
phyestresSe <- phyloseq_to_deseq2(physeqmayorSe, ~ estresSe)
phyestresSe <- estimateSizeFactors(phyestresSe, type = "poscounts")

# Creamos el objeto DESeq

phyestresSe <- DESeq(phyestresSe, test = "Wald", fitType = "parametric")

# Alto vs bajo

altobajoSe <- results(phyestresSe, cooksCutoff = FALSE, contrast = c("estresSe", "Estrés Alto", "Estrés Bajo"))

# Significancia 

alpha <- 0.05
altobajoSe <- altobajoSe[which(altobajoSe$padj < alpha &(altobajoSe$log2FoldChange >=1 |
                                                      altobajoSe$log2FoldChange<=-1)), ]

altobajoSe = cbind(as(altobajoSe, "data.frame"), as(tax_table(physeqbacsemillas16s)
                                                  [rownames(altobajoSe), ], "matrix"))

# Filtramos los que el género es NA para no representarlos

altobajoSe <- altobajoSe[!is.na(altobajoSe$Genus), ]

# Filtro para los cutibacterium con specie NA

cutibacterium_na <- altobajoSe %>%
  filter(Genus == "Cutibacterium", is.na(Species))



# Representamos

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(altobajoSe$log2FoldChange, altobajoSe$Genus, function(x) max(x))
x = sort(x, TRUE)
altobajoSe$Genus = factor(as.character(altobajoSe$Genus), levels=names(x))
altobajoSe$Species = factor(as.character(altobajoSe$Species), levels=names(x))
x = tapply(altobajoSe$log2FoldChange, altobajoSe$Family, function(x) max(x))
x = sort(x, TRUE)
altobajoSe$Family = factor(as.character(altobajoSe$Family), levels=names(x))
enrichplot <- ggplot(altobajoSe, aes(x = log2FoldChange, y = Genus, fill = Family)) +
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
    title = "Estrés Alto vs Estrés Bajo Semillas 16S",
    x = "log2 Fold Change",
    y = "Género"
  )

enrichplot

#############################################################################

