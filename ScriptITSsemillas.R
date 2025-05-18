# ITS semillas
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

list.files(path = "./Muestras/Nat-Bloqu05_merge/", pattern = ".*semillas.*ITS.*\\.fastq\\.gz$", 
           ignore.case = TRUE)
path<- './Muestras/Nat-Bloqu05_merge/'

fnFITSse <- sort(list.files(path, pattern="^_.*semillas.*ITS.*_R1.*\\.fastq\\.gz$", full.names = TRUE)) 
fnRITSse <- sort(list.files(path, pattern="^_.*semillas.*ITS.*_R2.*\\.fastq\\.gz$", full.names = TRUE)) 

# Para quedarnos con el nombre de las muestras

sample.namesITSse <- sapply(strsplit(basename(fnFITSse), "_"), `[`, 2)  

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

fnFITSse.filtN <- file.path(path, "filtN", basename(fnFITSse)) 
fnRITSse.filtN <- file.path(path, "filtN", basename(fnRITSse))
filterAndTrim(fnFITSse, fnFITSse.filtN, fnRITSse, fnRITSse.filtN, maxN = 0, multithread = FALSE)

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFITSse.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                              primerHits, fn = fnRITSse.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,                                                                                                                                                                               fn = fnFITSse.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRITSse.filtN[[1]]))
cutadapt<-"./Muestras/Nat-Bloqu05_merge/cutadapt.exe"

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFITSse.cut <- file.path(path.cut, basename(fnFITSse))
fnRITSse.cut <- file.path(path.cut, basename(fnRITSse))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFITSse)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFITSse.cut[i], "-p", fnRITSse.cut[i], 
                             fnFITSse.filtN[i], fnRITSse.filtN[i])) 
}

#comprobamos que hemos eliminado todo 

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFITSse.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                            primerHits, fn = fnRITSse.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                                           fn = fnFITSse.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRITSse.cut[[1]]))

cutFse <- sort(list.files(path.cut, pattern = "semillas_ITS_R1_all.fastq.gz", full.names = TRUE))
cutRse <- sort(list.files(path.cut, pattern = "semillas_ITS_R2_all.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) {
  sample_nameITSse <- sub("^.+cutadapt/_0*([0-9]+).*", "\\1", fname)
  return(sample_nameITSse)
}

sample.namesITSse <- unname(sapply(cutFse, get.sample.name))


# calidad forward

forwplotITSse<-ggplotly(plotQualityProfile(cutFse[1:length(cutFse)], aggregate=TRUE) + 
                          geom_hline(yintercept=c(15,25,35), 
                                     color=c("red","blue","green"), 
                                     size=0.5),
                        width =600)
forwplotITSse

# Calidad reverse 

revqplotITSse<-ggplotly(plotQualityProfile(cutRse[1:length(cutRse)], aggregate=TRUE) + 
                          geom_hline(yintercept=c(15,25,35), 
                                     color=c("red","blue","green"),
                                     size=0.5),
                        width =600)
revqplotITSse

# Filter and trimm

filtFITSse <- file.path(path.cut, "filteredITSse", basename(cutFse))
filtRITSse <- file.path(path.cut, "filteredITSse", basename(cutRse))
names(filtFITSse) <- sample.namesITSse
names(filtRITSse) <- sample.namesITSse
outsemillas <- filterAndTrim(cutFse, filtFITSse, cutRse, filtRITSse, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  

# Aprendemos errores

errFITSse <- learnErrors(filtFITSse, multithread=TRUE)
errRITSse <- learnErrors(filtRITSse, multithread=TRUE)

# Inferencia de las muestras 

dadaFITSse <- dada(filtFITSse, err = errFITSse, multithread = TRUE)
dadaRITSse <- dada(filtRITSse, err = errRITSse, multithread = TRUE)

# Merge de las secuencias 

mergersITSse <- mergePairs(dadaFITSse, filtFITSse, dadaRITSse, filtRITSse, verbose=TRUE)

seqtabITSse <- makeSequenceTable(mergersITSse)
dim(seqtabITSse)


# Eliminamos Quimeras

seqtab.nochimITSse <- removeBimeraDenovo(seqtabITSse, method="consensus", multithread=TRUE, verbose=TRUE)

sum(seqtab.nochimITSse)/sum(seqtabITSse) # Nos quedamos con el 94%

# Resumen del filtrado 

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFITSse, getN), sapply(dadaRITSse, getN), sapply(mergersITSse, getN),
               rowSums(seqtab.nochimITSse))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.namesITSse

# Asignamos taxonomía

taxaITSseasvs <- assignTaxonomy(seqtab.nochimITSse, "./Muestras/Nat-Bloqu05_merge/sh_general_release_dynamic_04.04.2024.fasta" , multithread = TRUE, tryRC = TRUE)

# Cargamos metadatos y modificamos

metadatossemillas <- read.xlsx("./Muestras/Nat-Bloqu-05/tax_otus_metadatos_semillas_12_09_24.xlsx", sheet = 1)

metadatossemillas <- metadatossemillas %>% 
  arrange(Samples)
samples.outITSse <- rownames(seqtab.nochimITSse)
subjectITSSe <- sapply(strsplit(samples.outITSse, "D"), `[`, 1)
estresITSSe <- metadatossemillas$ELC
ciudadITSSe <- metadatossemillas$Poblacion
samdfITSse <- data.frame(subjectITSSe=subjectITSSe, estresITSSe=estresITSSe, poblacionITSSe=ciudadITSSe)
samdfITSse$estrespoblacionITSse <- paste(estresITSSe,ciudadITSSe)
samdfITSse$estresITSSe[samdfITSse$estresITSSe == 1] <- "Estrés Alto"
samdfITSse$estresITSSe[samdfITSse$estresITSSe==2] <- "Estrés Medio"

samdfITSse$estresITSSe[samdfITSse$estresITSSe==3] <- "Estrés Bajo"
samdfITSse$estresITSSe[samdfITSse$estresITSSe=="Cultivar "] <- "CNF"
samdfITSse$estresITSSe[samdfITSse$estresITSSe=="Cultivar"] <- "C"
samdfITSse$estresITSSe[samdfITSse$estresITSSe=="Cultivar no fungicida"] <- "CNF"


rownames(samdfITSse) <- samples.outITSse

# Para igualar los nombres y no de error al crear el objeto phyloseq

rownames(samdfITSse) <- sprintf("%03d", as.integer(rownames(samdfITSse)))
otu <- otu_table(seqtab.nochimITSse, taxa_are_rows = FALSE)
sample_names(otu) <- rownames(seqtab.nochimITSse)
rownames(samdfITSse) <- rownames(seqtab.nochimITSse)
sampdata <- sample_data(samdfITSse)
tax <- tax_table(taxaITSseasvs)
sample_names(otu) <- sample_names(sampdata)
sample_names(tax) <- sample_names(otu)

# Creamos el objeto phyloseq

psITSSe <- phyloseq(otu, sampdata, tax)


# Objeto phyloseq limpio

physeqbacsemillasITS<-subset_taxa(psITSSe, Order!="Mitochondria" & Order!="Chloroplast") 


#################################################


# Representación

# Convertimos a data frame para modificar 
tax_df <- as.data.frame(tax_table(physeqbacsemillasITS))

# Eliminamos los taxones de Incertae por posible contaminación 

tax_to_keep <- !apply(tax_df, 1, function(x) any(grepl("Incertae", x, ignore.case = TRUE)))

# Filtramos el objeto phyloseq 

physeqbacsemillasITS <- prune_taxa(tax_to_keep, physeqbacsemillasITS)

# Representamos 

physeq.aggITSSe <- aggregate_rare(physeqbacsemillasITS %>% 
                                    transform(transform="compositional"),
                                  level="Genus", detection=0.05, prevalence=0.1)

# Eliminamos lo que va antes del "_"

taxa_names(physeq.aggITSSe) <- gsub(".*_", "", taxa_names(physeq.aggITSSe))
physeq.filt <- subset_samples(physeq.aggITSSe)


# Renombramos para la representación

taxa_names(physeq.filt)[taxa_names(physeq.filt) == "Other"] <- "Géneros <10% de abundancia"
taxa_names(physeq.filt)[taxa_names(physeq.filt) == "Unknown"] <- "NA"

getPalette <- colorRampPalette(brewer.pal(8, "Set2")) 
PhylaPalette <- getPalette(length(taxa(physeq.filt)))

# Graficamos

taxcompplot <- plot_composition(physeq.filt, average_by="poblacionITSSe",
                                x_label="poblacionITSSe", group_by="estresITSSe") +
  scale_y_percent() +
  scale_fill_manual(values = PhylaPalette)

taxcompplot$data <- taxcompplot$data %>%
  complete(xlabel, Tax, fill = list(Abundance = 0))

# Normalizamos para el total

taxcompplot$data <- taxcompplot$data %>%
  group_by(xlabel) %>%
  mutate(Abundance = Abundance / sum(Abundance)) 

taxcompplot$data$Tax <- factor(taxcompplot$data$Tax, levels = c(levels(taxcompplot$data$Tax)
                                                                [levels(taxcompplot$data$Tax)!="Unknown"], "Unknown"))


taxcompplot +
  theme(axis.text.x = element_text(size = 10),  
        axis.text.y = element_text(size = 10))




# Curvas de rarefracción

rarefactioncurveITSse<-rarecurve(as.matrix(seqtab.nochimITSse), step=100,las = 1, cex=0.75, tidy = TRUE) 
plotrarefaction <- ggplot(rarefactioncurveITSse, aes(Sample, Species)) + 
  geom_line(aes(color=Site)) +
  xlab("Coverage") + ylab("ASV number") +
  scale_x_continuous(breaks = seq(0, max(rarefactioncurveITSse$Sample), by = 70)) +  
  scale_y_continuous(breaks = seq(0, max(rarefactioncurveITSse$Species), by = 2)) +  
  theme(legend.position = "none")
plotrarefactionITSse <- ggplot(rarefactioncurveITSse, aes(Sample, Species)) + 
  geom_line(aes(color=Site)) +
  xlab("Coverage") + ylab("ASV number") +
  theme(legend.position = "none") 
ggplotly(plotrarefactionITSse)
rarefactioncurveITSse$Sample <- substr(rarefactioncurveITSse$Sample, 1, 10)

# Alfa diversidad 

divIdx = c("Chao1", "Shannon", "Simpson")

alphaplot_boxITSse <- plot_richness(physeqbacsemillasITS, x = "estresITSSe", measures = divIdx, color = "estresITSSe", nrow = 1) + 
  geom_boxplot(aes(color = estresITSSe, fill = estresITSSe), 
               alpha = 0.3,  
               size = 0.3,   
               show.legend = FALSE) +  
  geom_jitter(aes(color = estresITSSe), size = 0.6, width = 0.1, alpha = 0.7) +  
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Group", y = "Alpha Diversity Measure") + 
  stat_summary(aes(color = estresITSSe), fun = "median", geom = "point", shape = 18, size = 3)  
alphaplot_boxITSse


# Análisis ANOVA

physeq.filt <- subset_samples(physeqbacsemillasITS, 
                              !poblacionITSSe %in% c("OSK", "ZAPN"))
alpha_df <- estimate_richness(physeq.filt, measures = divIdx)
alpha_df$estres <- sample_data(physeq.filt)$estresITSSe

# Ejecutamos ANOVA

anova_results <- lapply(divIdx, function(index) {
  formula <- as.formula(paste(index, "~ estres"))
  model <- aov(formula, data = alpha_df)
  summary(model)
})
names(anova_results) <- divIdx
anova_results

########################################################

# Beta diversidad

# Filtramos 

physeq.filt <- subset_samples(physeq.aggITSSe, !poblacionITSSe %in% c("OSK", "ZAPN"))

ordinationNMDSbray<-ordinate(physeq.filt, method="NMDS", distance="bray")
plotordinationbray<-plot_ordination(physeq.filt, ordinationNMDSbray, color="estresITSSe",
                                    shape = "estresITSSe") +
  geom_point(size=2) + stat_ellipse(aes(group = estresITSSe), type = "t") +
  theme_bw() + labs(title = "Non rarefied Bray Curtis NMDS")
plotordinationbray

# No tenemos muestras suficientes muestras por lo que no podemos hacerlo 


# Análisis de enriquecimiento

filtersampITSSe<-genefilter_sample(physeqbacsemillasITS, filterfun_sample(function(x) x > 1),
                                A=0.01*nsamples(physeqbacsemillasITS))
physeqmayorITSSe <- prune_taxa(filtersampITSSe, physeqbacsemillasITS)

physeqmayorITSSe <- subset_samples(physeqmayorITSSe, 
                                   !estresITSSe %in% c("C"))

table(sample_data(physeqmayorITSSe)$estresITSSe)

phyestresITSSe <- phyloseq_to_deseq2(physeqmayorITSSe, ~ estresITSSe)
phyestresITSSe <- estimateSizeFactors(phyestresITSSe, type = "poscounts")

# No podemos crear el objeto phyloseq por falta de datos 

phyestresITSSe<- DESeq(phyestresITSSe, test="Wald", fitType="parametric")


