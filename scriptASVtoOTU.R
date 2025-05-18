# Para las correlaciones pasamos de ASVs a OTUS
 
library(tibble)
library(dplyr)
library(DECIPHER)
library(Biostrings)
library(tibble)
library(iNEXT)
library(openxlsx)
library(phyloseq)
library(microbiome)
library(readxl)
library(tidyverse)
library(magrittr)
library(ggnested)
library(knitr)
library(gridExtra)
library(scales)
library(ggplot2)
library(vegan)
library(tidyr)
library(dplyr)
library(ggalt)
library(wesanderson)
library(RColorBrewer)
library(vegan)
library(corrplot)
library(ggpubr)
library(car)
library(cowplot)
library(lattice)
library(ggpubr)
library(ADER)
library(patchwork)

# Numero de nucleos a utilizar

nproc <- 6


# SUELOS 16S

# Tabla de conteos de ASVs de suelos 16S

seqtab <- seqtab.nochimfiltered

# Obtenemos secuencias para el alineamiento 

asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

## Alineamos las secuencias

aln <- DECIPHER::AlignSeqs(dna, processors = nproc)

# Matriz de distancias 

d <- DECIPHER::DistanceMatrix(aln, processors = nproc)

# Agrupamos ASVs por 97% de similitud

clusters <- DECIPHER::TreeLine(
  myDistMatrix=d,
  method = "complete",
  cutoff = 0.03, 
  type = "clusters",
  processors = nproc,
  showPlot = TRUE,
  verbose = TRUE)

# Vemos el número de clusters

summary(clusters$cluster)
length(unique(clusters$cluster))

# Creamos la nueva tabla de OTUs

merged_seqtab <- seqtab %>% 
  t %>%
  rowsum(clusters$cluster) %>%
  t

# Los nombramos como OTU1, OTU2...

colnames(merged_seqtab) <- paste0("OTU", colnames(merged_seqtab))


# Asociamos ASVs con OTUS 

asv_to_otu <- data.frame(
  ASV = asv_sequences,
  OTU = clusters$cluster
)

# Elegimos secuencia más representativa 

find_centroid <- function(otu_group) {
  otu_asvs <- otu_group$ASV
  otu_indices <- which(asv_to_otu$ASV %in% otu_asvs)
  dist_subset <- d[otu_indices, otu_indices]
  avg_dist <- colMeans(as.matrix(dist_subset))
  centroid_index <- which.min(avg_dist)
  return(otu_asvs[centroid_index])
}

# Aplicamos la función 

representative_otus <- asv_to_otu %>%
  group_by(OTU) %>%
  summarize(representative_sequence = find_centroid(cur_data()), .groups = 'drop')


# Creamos un FASTA para las secuencias más representativas 

otu_dna <- DNAStringSet(representative_otus$representative_sequence)
names(otu_dna) <- paste0("OTU", representative_otus$OTU)

Biostrings::writeXStringSet(otu_dna, filepath = "./OTUs_sequences.fasta")

# Asignamos de nuevo taxonomía 

taxaSuOTUs <- assignTaxonomy("OTUs_sequences.fasta", "./silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)

# Crearemos el objeto phyloseq de nuevo 

#Cargamos los archivos

seqtable_nochim_otu_suelos <- merged_seqtab
taxonomy_otu_suelos <- taxaSuOTUs

# Modificamos

rownames(taxonomy_otu_suelos) <- colnames(seqtable_nochim_otu_suelos)
rownames(taxonomy_otu_suelos)
class(taxonomy_otu_suelos)
class(seqtable_nochim_otu_suelos)
t_seqtable_nochim_otu_suelos <- t(seqtable_nochim_otu_suelos) 
m_taxonomy_otu_suelos <- as.matrix(taxonomy_otu_suelos) 
class(m_taxonomy_otu_suelos)
class(t_seqtable_nochim_otu_suelos)


# Cargamos metadatos y modificamos

metadata_suelos <- sample_data(read.xlsx("./Muestras/Nat-Bloqu-05/tax_otus_metadatos_suelo_12_09_24B.xlsx",
                                         rowNames = TRUE, sheet = "Metadatos"))
OTU_SUE <- otu_table(t_seqtable_nochim_otu_suelos, taxa_are_rows = TRUE)
TAX_OTU_SUE <- tax_table(m_taxonomy_otu_suelos)

# Creamos el objeto phyloseq 

phylobject_otu_sue<- phyloseq(OTU_SUE, TAX_OTU_SUE, metadata_suelos)

# Limpiamos 

phylobject_otu_sue<-subset_taxa(phylobject_otu_sue, Order!="Mitochondria" & Order!="Chloroplast") 




# Semillas 16S

# Modificamos nombres 

rownames(seqtab.nochimSe) <- gsub("^0", "", rownames(seqtab.nochimSe))

# Cargamos tabla de asvs

seqtabSe <- seqtab.nochimSe

# Obtenemos secuencias para el alineamiento 

asv_sequencesSe <- colnames(seqtabSe)
sample_namesSe <- rownames(seqtabSe)
dnaSe <- Biostrings::DNAStringSet(asv_sequencesSe)


# Alineamos

alnSe <- DECIPHER::AlignSeqs(dnaSe, processors = nproc)

# Matriz de distancias 

dSe <- DECIPHER::DistanceMatrix(alnSe, processors = nproc)


# Agrupamos ASVs por 97% de similitud

clustersSe <- DECIPHER::TreeLine(
  myDistMatrix=dSe,
  method = "complete",
  cutoff = 0.03, 
  type = "clusters",
  processors = nproc,
  showPlot = TRUE,
  verbose = TRUE)


# Verificamos numero de clusters

summary(clustersSe$cluster)
length(unique(clustersSe$cluster))

# Creamos nueva tabla de OTUS

merged_seqtabSe <- seqtabSe %>% 
  t %>%
  rowsum(clustersSe$cluster) %>%
  t

# Renombramos los OTUS

colnames(merged_seqtabSe) <- paste0("OTU", colnames(merged_seqtabSe))



# Asociamos OTUS y ASVs
asv_to_otuSe <- data.frame(
  ASV = asv_sequencesSe,
  OTU = clustersSe$cluster
)


# Función para encontrar la secuencia más representativa

find_centroid <- function(otu_group) {
  otu_asvs <- otu_group$ASV
  otu_indices <- which(asv_to_otuSe$ASV %in% otu_asvs)
  dist_subset <- dSe[otu_indices, otu_indices]
  avg_dist <- colMeans(as.matrix(dist_subset))
  centroid_index <- which.min(avg_dist)
  return(otu_asvs[centroid_index])
}

# Aplicamos la función 

representative_otusSe <- asv_to_otuSe %>%
  group_by(OTU) %>%
  summarize(representative_sequence = find_centroid(cur_data()), .groups = 'drop')
otu_dnaSe <- DNAStringSet(representative_otusSe$representative_sequence)
names(otu_dnaSe) <- paste0("OTU", representative_otusSe$OTU)

# Creamos FASTA 

Biostrings::writeXStringSet(otu_dnaSe, filepath = "./OTUs_sequencessemillas.fasta")

# Asignamos taxonomía 

taxaSeOTUs <- assignTaxonomy("OTUs_sequencessemillas.fasta", "./silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)



# Cargamos los archivos
rownames(merged_seqtabSe) <- sub("^0+", "", rownames(merged_seqtabSe))

seqtable_nochim_otu_semillas <- merged_seqtabSe
taxonomy_otu_semillas <- taxaSeOTUs

# Modificamos 

rownames(taxonomy_otu_semillas) <- colnames(seqtable_nochim_otu_semillas)
rownames(taxonomy_otu_semillas)
class(taxonomy_otu_semillas)
class(seqtable_nochim_otu_semillas)
t_seqtable_nochim_otu_semillas <- t(seqtable_nochim_otu_semillas) 
m_taxonomy_otu_semillas <- as.matrix(taxonomy_otu_semillas) 
class(m_taxonomy_otu_semillas)
class(t_seqtable_nochim_otu_semillas)

# Cargamos metadatos

metadatossemillasOTUS <- sample_data(read.xlsx("./Muestras/Nat-Bloqu-05/tax_otus_metadatos_semillas_12_09_24.xlsx",
                                               rowNames = FALSE, sheet = "metadatos"))


OTU_SEM <- otu_table(t_seqtable_nochim_otu_semillas, taxa_are_rows = TRUE)
TAX_OTU_SEM <- tax_table(m_taxonomy_otu_semillas)

# Creamos phyloseq

phylobject_otu_sem <- phyloseq(OTU_SEM, TAX_OTU_SEM, metadatossemillasOTUS)

# Limpiamos
phylobject_otu_sem<-subset_taxa(phylobject_otu_sem, Order!="Mitochondria" & Order!="Chloroplast") 

#########


#SUELOS ITS

# Modificamos Nombres

rownames(seqtabITSsu) <- gsub("^_(\\d+)_.*", "\\1", rownames(seqtabITSsu))
rownames(seqtabITSsu) <- gsub("^0", "", rownames(seqtabITSsu))


# Cargamos la tabla de ASVs

seqtabITSsu <- seqtab.nochimITSsu

# Obtenemos secuencias para el alineamiento 

asv_sequencesITSsu <- colnames(seqtabITSsu)
sample_namesITSsu <- rownames(seqtabITSsu)
dnaITSsu <- Biostrings::DNAStringSet(asv_sequencesITSsu)


# Alineamos 

alnITSsu <- DECIPHER::AlignSeqs(dnaITSsu, processors = nproc)

# Matriz de distancias

dITSsu <- DECIPHER::DistanceMatrix(alnITSsu, processors = nproc)

# Buscamos clusters OTUS con 97%
clustersITSsu <- DECIPHER::TreeLine(
  myDistMatrix=dITSsu,
  method = "complete",
  cutoff = 0.03, 
  type = "clusters",
  processors = nproc,
  showPlot = TRUE,
  verbose = TRUE)


# Verificamos el numero de clusters

summary(clustersITSsu$cluster)
length(unique(clustersITSsu$cluster))

# Creamos nueva tabla de OTUS

merged_seqtabITSsu <- seqtabITSsu %>% 
  t %>%
  rowsum(clustersITSsu$cluster) %>%
  t


# Renombramos OTUS

colnames(merged_seqtabITSsu) <- paste0("OTU", colnames(merged_seqtabITSsu))



# Asociamos OTUS y ASVs

asv_to_otuITSsu <- data.frame(
  ASV = asv_sequencesITSsu,
  OTU = clustersITSsu$cluster
)

# Función para encontrar la secuencia más representativa

find_centroid <- function(otu_group) {
  otu_asvs <- otu_group$ASV
  otu_indices <- which(asv_to_otuITSsu$ASV %in% otu_asvs)
  dist_subset <- dITSsu[otu_indices, otu_indices]
  avg_dist <- colMeans(as.matrix(dist_subset))
  centroid_index <- which.min(avg_dist)
  return(otu_asvs[centroid_index])
}

# Aplicamos la función 

representative_otusITSsu <- asv_to_otuITSsu %>%
  group_by(OTU) %>%
  do({
    otu_group <- .
    centroid <- find_centroid(otu_group)
    tibble(representative_sequence = centroid)
  }) %>%
  ungroup()
otu_dnaITSsu <- DNAStringSet(representative_otusITSsu$representative_sequence)
names(otu_dnaITSsu) <- paste0("OTU", representative_otusITSsu$OTU)

#Creamos FASTA 

Biostrings::writeXStringSet(otu_dnaITSsu, filepath = "./OTUs_sequencesITSsu.fasta")

# Asignamos taxonomía

taxaITSsuOTUS <- assignTaxonomy("OTUs_sequencesITSsu.fasta", "./Muestras/Nat-Bloqu05_merge/sh_general_release_dynamic_04.04.2024.fasta", multithread=TRUE)


seqtable_nochim_otu_ITSsuelo <- merged_seqtabITSsu
taxonomy_otu_ITSsuelo <- taxaITSsuOTUS



# Modificamos 

rownames(taxonomy_otu_ITSsuelo) <- colnames(seqtable_nochim_otu_ITSsuelo)
rownames(taxonomy_otu_ITSsuelo)
class(taxonomy_otu_ITSsuelo)
class(seqtable_nochim_otu_ITSsuelo)
t_seqtable_nochim_otu_ITSsuelo <- t(seqtable_nochim_otu_ITSsuelo) 
m_taxonomy_otu_ITSsuelo <- as.matrix(taxonomy_otu_ITSsuelo) 
class(m_taxonomy_otu_ITSsuelo)
class(t_seqtable_nochim_otu_ITSsuelo)



OTU_ITSSU <- otu_table(t_seqtable_nochim_otu_ITSsuelo, taxa_are_rows = TRUE)
TAX_OTU_ITSSU <- tax_table(m_taxonomy_otu_ITSsuelo)

# Cargamos metadatos

metadata_suelos <- sample_data(read.xlsx("./Muestras/Nat-Bloqu-05/tax_otus_metadatos_suelo_12_09_24B.xlsx",
                                         rowNames = TRUE, sheet = "Metadatos"))
sample_names(OTU_ITSSU) <- gsub(".*_(\\d+)_suelos_ITS.*", "\\1", sample_names(OTU_ITSSU))

sample_names(OTU_ITSSU) <- sub("^0+", "", sample_names(OTU_ITSSU))


# Creamos el objeto phyloseq

phylobject_otu_ITSsu<- phyloseq(OTU_ITSSU, TAX_OTU_ITSSU, metadata_suelos)
phylobject_otu_ITSsu<-subset_taxa(phylobject_otu_ITSsu, Order!="Mitochondria" & Order!="Chloroplast") 

# Eliminamos incertae

phylobject_otu_ITSsu <- subset_taxa(
  phylobject_otu_ITSsu,
  !grepl("incertae", Order, ignore.case = TRUE) &
    Order != "Mitochondria" &
    Order != "Chloroplast"
)

#####################################################################################################