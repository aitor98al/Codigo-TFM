## ITS SUELOS 
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

list.files(path = "./Muestras/Nat-Bloqu05_merge/", pattern = ".*suelos.*ITS.*\\.fastq\\.gz$", 
           ignore.case = TRUE)
path<- './Muestras/Nat-Bloqu05_merge/'

fnFITSsu <- sort(list.files(path, pattern="^_.*suelos.*ITS.*_R1.*\\.fastq\\.gz$", full.names = TRUE)) 
fnRITSsu <- sort(list.files(path, pattern="^_.*suelos.*ITS.*_R2.*\\.fastq\\.gz$", full.names = TRUE))

# Nos quedamos solo con los nombres de las muestras
sample.namesITSsu <- sapply(strsplit(basename(fnFITSsu), "_"), `[`, 2)  

# Cutadapt

FWD <- "AACTTTYRRCAAYGGATCWCT"  
REV <- "AGCCTCCGCTTATTGATATGCTTAART"

allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFITSsu.filtN <- file.path(path, "filtN", basename(fnFITSsu)) 
fnRITSsu.filtN <- file.path(path, "filtN", basename(fnRITSsu))
filterAndTrim(fnFITSsu, fnFITSsu.filtN, fnRITSsu, fnRITSsu.filtN, maxN = 0, multithread = FALSE)

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFITSsu.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                          primerHits, fn = fnRITSsu.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                                       fn = fnFITSsu.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRITSsu.filtN[[1]]))
# Cutadapt 

cutadapt<-"./Muestras/Nat-Bloqu05_merge/cutadapt.exe"
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFITSsu.cut <- file.path(path.cut, basename(fnFITSsu))
fnRITSsu.cut <- file.path(path.cut, basename(fnRITSsu))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFITSsu)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, 
                             "-o", fnFITSsu.cut[i], "-p", fnRITSsu.cut[i], 
                             fnFITSsu.filtN[i], fnRITSsu.filtN[i])) 
}


# Comprobamos

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFITSsu.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                        primerHits, fn = fnRITSsu.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                                 fn = fnFITSsu.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRITSsu.cut[[1]]))
cutFs <- sort(list.files(path.cut, pattern = "_R1_all.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_all.fastq.gz", full.names = TRUE))

# Para el nombre de las muestras

get.sample.name <- function(fname) {
  sample_name <- sub("^.+cutadapt/_0*([0-9]+).*", "\\1", fname)
  return(sample_name)
}
sample.names <- unname(sapply(cutFs, get.sample.name))


# Calidad forwards

forwplotITSsu<-ggplotly(plotQualityProfile(cutFs[1:length(cutFs)], aggregate=TRUE) + 
                       geom_hline(yintercept=c(15,25,35), 
                                  color=c("red","blue","green"), 
                                  size=0.5),
                     width =600)
forwplotITSsu

# Calidad Reverse

revqplotITSsu<-ggplotly(plotQualityProfile(cutRs[1:length(cutRs)], aggregate=TRUE) + 
                       geom_hline(yintercept=c(15,25,35), 
                                  color=c("red","blue","green"),
                                  size=0.5),
                      600)
revqplotITSsu

# Filter and trimm

filtFITSsu <- file.path(path.cut, "filteredITSsuelo", basename(cutFs))
filtRITSsu <- file.path(path.cut, "filteredITSsuelo", basename(cutRs))
names(filtFITSsu) <- sample.namesITSsu
names(filtRITSsu) <- sample.namesITSsu
out <- filterAndTrim(cutFs, filtFITSsu, cutRs, filtRITSsu, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE) 

# Aprendemos error

errFITSsu <- learnErrors(filtFITSsu, multithread=TRUE)
errRITSsu <- learnErrors(filtRITSsu, multithread=TRUE)

# Inferencia de las muestras

dadaFITSsu <- dada(filtFITSsu, err = errFITSsu, multithread = TRUE)
dadaRITSsu <- dada(filtRITSsu, err = errRITSsu, multithread = TRUE)

# Merge de las secuencias

mergersITSsu <- mergePairs(dadaFITSsu, filtFITSsu, dadaRITSsu, filtRITSsu, verbose=TRUE)

seqtabITSsu <- makeSequenceTable(mergersITSsu)
dim(seqtabITSsu)

# Eliminamos quimeras

seqtab.nochimITSsu <- removeBimeraDenovo(seqtabITSsu, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochimITSsu)/sum(seqtabITSsu) # nos quedamos con el 98%


# Resumen filtrado

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFITSsu, getN), sapply(dadaRITSsu, getN), sapply(mergersITSsu, getN),
               rowSums(seqtab.nochimITSsu))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

## Asignamos Taxonomía

taxaITSsuasvs <- assignTaxonomy(seqtab.nochimITSsu, "./Muestras/Nat-Bloqu05_merge/sh_general_release_dynamic_04.04.2024.fasta" , multithread = TRUE, tryRC = TRUE)

# Cargamos metadatos y modificamos

metadatosITSsu <- read.xlsx("./Muestras/Nat-Bloqu-05/tax_otus_metadatos_suelo_12_09_24B.xlsx", sheet = 1)
metadatosITSsu <- metadatosITSsu %>% 
  arrange(Column1)
rownames(seqtab.nochimITSsu) <- sample.namesITSsu
samples.outITSsu <- rownames(seqtab.nochimITSsu)
subjectITSsu <- sample.namesITSsu
estresITSsu <- metadatosITSsu$ELC
samdfITSsu <- data.frame(subjectITSsu=subjectITSsu, estresITSsu=estresITSsu)
ciudadITSSu <- metadatosITSsu$Poblacion
samdfITSsu <- data.frame(subjectITSsu=subjectITSsu, estresITSsu=estresITSsu, poblacionITSsu=ciudadITSSu)
samdfITSsu$estrespoblacionITSsu <- paste(estresITSsu,ciudadITSSu)
samdfITSsu$estresITSsu[samdfITSsu$estresITSsu == 1] <- "Estrés Alto"
samdfITSsu$estresITSsu[samdfITSsu$estresITSsu==2] <- "Estrés Medio"

samdfITSsu$estresITSsu[samdfITSsu$estresITSsu==3] <- "Estrés Bajo"
rownames(samdfITSsu) <- samples.outITSsu

# Creamos phyloseq

psITSsu <- phyloseq(otu_table(seqtab.nochimITSsu, taxa_are_rows=FALSE), 
                 sample_data(samdfITSsu), 
                 tax_table(taxaITSsuasvs))

# Eliminamos posibles contaminaciones

physeqbacITSsu<-subset_taxa(psITSsu, Order!="Mitochondria" & Order!="Chloroplast") # que no sea mitocondria ni cloroplastos

# Representación taxonomia

physeq.aggITSSu <- aggregate_rare(physeqbacITSsu %>% 
                                    transform(transform="compositional"),
                                  level="Genus", detection=0.05, prevalence=0.1)

# Eliminamos a nivel de género el "_"
tax_table(physeq.aggITSSu)[, "Genus"] <- gsub("^g[_]*", "", tax_table(physeq.aggITSSu)[, "Genus"])
taxa_names(physeq.aggITSSu) <- tax_table(physeq.aggITSSu)[, "Genus"]


physeq.filt <- subset_samples(physeq.aggITSSu)

# Renombramos para la representación

taxa_names(physeq.filt)[taxa_names(physeq.filt) == "Other"] <- "Géneros <10% de abundancia"
taxa_names(physeq.filt)[taxa_names(physeq.filt) == "Unknown"] <- "NA"

getPalette <- colorRampPalette(brewer.pal(8, "Set2")) 
PhylaPalette <- getPalette(length(taxa(physeq.filt)))

# Representamos 

taxcompplot <- plot_composition(physeq.filt, average_by="poblacionITSsu",
                                x_label="poblacionITSsu", group_by="estresITSsu") +
  scale_y_percent() +
  scale_fill_manual(values = PhylaPalette)

taxcompplot +
  theme(axis.text.x = element_text(size = 10),  
        axis.text.y = element_text(size = 10)) 




# Curvas de rarefracción 

rarefactioncurveITSsu<-rarecurve(as.matrix(seqtab.nochimITSsu), step=100,las = 1, cex=0.75, tidy = TRUE)  
plotrarefaction <- ggplot(rarefactioncurveITSsu, aes(Sample, Species)) + 
  geom_line(aes(color=Site)) +
  xlab("Coverage") + ylab("ASV number") +
  scale_x_continuous(breaks = seq(0, max(rarefactioncurveITSsu$Sample), by = 70)) +  
  scale_y_continuous(breaks = seq(0, max(rarefactioncurveITSsu$Species), by = 2)) +  
  theme(legend.position = "none")
plotrarefactionITSsu <- ggplot(rarefactioncurveITSsu, aes(Sample, Species)) + 
  geom_line(aes(color=Site)) +
  xlab("Coverage") + ylab("ASV number") +
  theme(legend.position = "none") 
ggplotly(plotrarefactionITSsu)


# Alfa diversidad

divIdx = c("Chao1", "Shannon", "Simpson")

alphaplot_boxITS <- plot_richness(physeqbacITSsu, x = "estresITSsu", measures = divIdx, color = "estresITSsu", nrow = 1) + 
  geom_boxplot(aes(color = estresITSsu, fill = estresITSsu), 
               alpha = 0.3,  
               size = 0.3,   
               show.legend = FALSE) +  
  geom_jitter(aes(color = estresITSsu), size = 0.6, width = 0.1, alpha = 0.7) +  
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Group", y = "Alpha Diversity Measure") + 
  stat_summary(aes(color = estresITSsu), fun = "median", geom = "point", shape = 18, size = 3)
alphaplot_boxITS


# Análisis ANOVA

alpha_df <- estimate_richness(physeqbacITSsu, measures = divIdx)
alpha_df$estres <- sample_data(physeqbacITSsu)$estresITSsu

anova_results <- lapply(divIdx, function(index) {
  formula <- as.formula(paste(index, "~ estres"))
  model <- aov(formula, data = alpha_df)
  summary(model)
})
names(anova_results) <- divIdx
anova_results

# Análisis TUKEY

model <- aov(Chao1 ~ estres, data = alpha_df)
TukeyHSD(model)

# Beta diversidad

ordinationNMDSbrayITS<- ordinate(physeqbacITSsu, method = "NMDS", distance = "bray")

# Ploteamos

plotordinationbrayNMDSITS <- plot_ordination(physeqbacITSsu, ordinationNMDSbrayITS, color = "estresITSsu", shape = "estresITSsu") + 
  geom_point(size = 2) +  
  stat_ellipse(aes(group = estresITSsu), type = "t") +  
  theme_bw()   
plotordinationbrayNMDSITS


# Análisis enriquecimiento

filtersamp<-genefilter_sample(physeqbacITSsu, filterfun_sample(function(x) x > 2),
                              A=0.3*nsamples(physeqbacITSsu)) # 0.3 por que salen muchas 
physeqmayor <- prune_taxa(filtersamp, physeqbacITSsu)
phyestres <- phyloseq_to_deseq2(physeqmayor, ~ estresITSsu)
phyestres<- DESeq(phyestres, test="Wald", fitType="parametric")

# Quitamos prefijos 

rownames(phyestres) <- gsub("^[a-z]__", "", rownames(phyestres))  

# Alto vs bajo

altobajo <- results(phyestres, cooksCutoff = FALSE, contrast = c("estresITSsu", "Estrés Alto", "Estrés Bajo"))
alpha <- 0.05
altobajo <- altobajo[which(altobajo$padj < alpha &(altobajo$log2FoldChange >=1 |
                                                      altobajo$log2FoldChange<=-1)), ]

altobajo = cbind(as(altobajo, "data.frame"), as(tax_table(physeqbacITSsu)
                                                  [rownames(altobajo), ], "matrix"))
# Eliminamos prefijos

altobajo$Genus <- gsub("^[a-z]__+", "", altobajo$Genus)
altobajo$Family <- gsub("^[a-z]__+", "", altobajo$Family)

# Para filtrar penicillium y buscar la especie

altobajo$Genus <- as.character(altobajo$Genus)
penicillium_na_idx <- which(altobajo$Genus == "Penicillium" & is.na(altobajo$Species))
altobajo$Genus[penicillium_na_idx[1]] <- "Penicillium restrictum"
altobajo$Genus[penicillium_na_idx[2]] <- "Penicillium onobense"

# Representamos

n_families <- length(unique(altobajo$Family))
palette_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_families)
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
  scale_fill_manual(values = palette_colors) +
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
    title = "Estrés Alto vs Estrés Bajo Suelo ITS",
    x = "log2 Fold Change",
    y = "Género"
  )

enrichplot

###########################################################################################
