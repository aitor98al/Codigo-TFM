#Script correlaciones

library(phyloseq)
library(dplyr)
library(tidyr)



# SUELOS 16s

# Cargamos los datos y modificamos

otu_table_suelos <- otu_table(phylobject_otu_sue)

# Convertimos en data frame

otu_df_suelos <- as.data.frame(otu_table_suelos)
otu_taxonomy_suelos <- tax_table(phylobject_otu_sue)
otu_taxonomy_suelos <- as.data.frame(otu_taxonomy_suelos)
otu_df_suelos <- as.data.frame(t(otu_table_suelos))
otu_idsuelo <- colnames(otu_df_suelos)
genero <- otu_taxonomy_suelos$Genus[match(otu_idsuelo, rownames(otu_taxonomy_suelos))]
familia <- otu_taxonomy_suelos$Family[match(otu_idsuelo, rownames(otu_taxonomy_suelos))]

colnames(otu_df_suelos) <- genero

# Eliminamos NAs

colnames(otu_df_suelos) <- ifelse(is.na(genero), otu_idsuelo, genero)


# Semillas 16S

# Cargamos los datos y modificamos

otu_table_semillas <- otu_table(phylobject_otu_sem)

# Convertimos en data frame

otu_df_semillas <- as.data.frame(otu_table_semillas)
otu_taxonomy_semillas <- tax_table(phylobject_otu_sem)
otu_taxonomy_semillas <- as.data.frame(otu_taxonomy_semillas)
otu_df_semillas <- as.data.frame(t(otu_table_semillas))
otu_idsemillas <- colnames(otu_df_semillas)
generosemillas <- otu_taxonomy_semillas$Genus[match(otu_idsemillas, rownames(otu_taxonomy_semillas))]
familiasemillas <- otu_taxonomy_semillas$Family[match(otu_idsemillas, rownames(otu_taxonomy_semillas))]
colnames(otu_df_semillas) <- generosemillas

# Eliminamos NAs

colnames(otu_df_semillas) <- ifelse(is.na(generosemillas), otu_idsemillas, generosemillas)


# Análisis por estreses

#Estrés alto

# Suelos 16s

otu_df_suelosalto<- otu_df_suelos[rownames(otu_df_suelos) %in% c("83","84","85","91","92","93","97","98","99","100","101","102","103","104","107", "108", "109","110","111","112","116","120","121","122","123","124","125"), ]

# Semillas 16S

otu_df_semillasalto<-otu_df_semillas[rownames(otu_df_semillas) %in% c ("1","2","3","4","5","6","10","11","12","13","14","15","16","17","18","25","26","27", "37","38","39","55","56","57","58","59","60")]

nombres_suelosalto <- colnames(otu_df_suelosalto)
nombres_semillasalto <- colnames(otu_df_semillasalto)


# Correlación Suelos 16S estrés alto 

# Método spearman

cor_matrixalto <- cor(otu_df_suelosalto, method = "spearman")  

#Función para calcular matriz de p-valores

cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}

# Matriz de p-valores aplicando la función

p_matrixsuelosalto <- cor_pmat(otu_df_suelosalto, method = "spearman")

# Convertimos a vector para aplicar corrección
p_vec_suelosalto <- as.vector(p_matrixsuelosalto)

# Aplicamos corrección
p_adj_vec_suelosalto <- p.adjust(p_vec_suelosalto, method = "BH")
p_matrixsuelosalto <- matrix(p_adj_vec_suelosalto, nrow = nrow(p_matrixsuelosalto), ncol = ncol(p_matrixsuelosalto))


# Convertimos las matrices a formato largo

cor_df <- reshape2::melt(cor_matrixalto)
pval_df <- reshape2::melt(p_matrixsuelosalto)

# Combinamos en un data frame

cor_pval_df <- cbind(cor_df, pval = pval_df$value)

# Filtramos 

cor_pval_fuerte_sig <- subset(cor_pval_df, value > 0.8 & value <0.99 & pval > 0 & pval < 0.05 & Var1 != Var2)
cor_pval_fuerte_sig <- cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ]
cor_pval_fuerte_sig_filtradosuelos16s <- subset(cor_pval_fuerte_sig, 
                                       !grepl("^OTU", Var1) & !grepl("^OTU", Var2))

# Creamos un nuevo grafo

gsuelos <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosuelos16s, directed = FALSE)

# Calculamos variables

transitivity(gsuelos, type = "global")
transitivity(gsuelos, type = "local")
btw <- betweenness(gsuelos, directed = FALSE)
V(gsuelos)$betweenness <- btw
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
edge_density(gsuelos)

# Clusterizacion

transitivity(gsuelos, type = "global")

# Grado
V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)

# Puntos de articulación

articulation_points(gsuelos)

# Ponemos capa para una mejor visualización

lay <- layout_with_fr(gsuelos, niter = 1000)

# Paleta de colores

n_clusters <- length(unique(V(gsuelos)$cluster))
colores <- brewer.pal(min(n_clusters, 8), "Set2")

# Ploteamos

plot(gsuelos,
     layout = lay,
     vertex.color = colores[V(gsuelos)$cluster],
     vertex.size = log(V(gsuelos)$degree + 1) * 5,
     vertex.label.cex = 0.6,
     vertex.label.color = "black", 
     edge.width = E(gsuelos)$value * 2,
     main = "Red de correlación suelos 16s estrés alto")

# Creamos tabla de métricas

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$betweenness <- betweenness(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)

# Clusterizacion
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)

# Puntos de articulación

art_points <- articulation_points(gsuelos)
V(gsuelos)$articulacion <- V(gsuelos)$name %in% V(gsuelos)[art_points]$name

# Creamos un data frame y guardamos

tabla_metricassuelosestresalto <- data.frame(
  Nodo = V(gsuelos)$name,
  Grado = V(gsuelos)$degree,
  Betweenness = V(gsuelos)$betweenness,
  Closeness = V(gsuelos)$closeness,
  Cluster = V(gsuelos)$cluster,
  Articulacion = V(gsuelos)$articulacion
)

write.table(tabla_metricassuelosestresalto, 
            file = "tabla_metricassuelosestresalto.tsv",  
            sep = "\t",  
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE) 
##########################################################################################

# Estrés medio

# Suelos 16s

otu_df_suelosmedio<- otu_df_suelos[rownames(otu_df_suelos) %in% c("70","71","72","73","74","75","76","77","78","79","88","89","90","117","118", "119"), ]

# Semillas 16S

otu_df_semillasmedio<-otu_df_semillas[rownames(otu_df_semillas) %in% c ("22","23","24","46","47","48","49","50","51","52","53","54","61","62","63")]

nombres_suelosmedio <- colnames(otu_df_suelosmedio)
nombres_semillasmedio <- colnames(otu_df_semillasmedio)


# Correlación suelos medio

cor_matrixmedio <- cor(otu_df_suelosmedio, method = "spearman") 
cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}

# Matriz de p-valores

p_matrixsuelosmedio <- cor_pmat(otu_df_suelosmedio, method = "spearman")

# Convertimos a vector para aplicar corrección
p_vec_suelosmedio <- as.vector(p_matrixsuelosmedio)

# Aplicamos corrección
p_adj_vec_suelosmedio <- p.adjust(p_vec_suelosmedio, method = "BH")
p_matrixsuelosmedio <- matrix(p_adj_vec_suelosmedio, nrow = nrow(p_matrixsuelosmedio), ncol = ncol(p_matrixsuelosmedio))

# Convertimos las matrices a formato largo

cor_df <- reshape2::melt(cor_matrixmedio)
pval_df <- reshape2::melt(p_matrixsuelosmedio)

# Combinamos 

cor_pval_df <- cbind(cor_df, pval = pval_df$value)

# Filtramos

cor_pval_fuerte_sig <- subset(cor_pval_df, value > 0.8 & value <0.99 & pval > 0 & pval < 0.05 & Var1 != Var2)
cor_pval_fuerte_sig <- cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ]
cor_pval_fuerte_sig_filtradosuelosmedio <- subset(cor_pval_fuerte_sig, 
                                                !grepl("^OTU", Var1) & !grepl("^OTU", Var2))



# Creamos Grafo

gsuelos <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosuelosmedio, directed = FALSE)

# Calculamos variables

transitivity(gsuelos, type = "global")
transitivity(gsuelos, type = "local")
btw <- betweenness(gsuelos, directed = FALSE)
V(gsuelos)$betweenness <- btw
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
edge_density(gsuelos)

# Clusterización

transitivity(gsuelos, type = "global")

# Grado

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)

# Puntos de articulación

articulation_points(gsuelos)

# Capa para mejor visualización

lay <- layout_with_fr(gsuelos, niter = 1000)

# Paleta 
n_clusters <- length(unique(V(gsuelos)$cluster))
colores <- brewer.pal(min(n_clusters, 8), "Set2")

# Graficamos

plot(gsuelos,
     layout = lay,
     vertex.color = colores[V(gsuelos)$cluster],
     vertex.size = log(V(gsuelos)$degree + 1) * 5,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",  
     edge.width = E(gsuelos)$value * 2,
     main = "Red de correlación suelos 16s estrés medio")


# Tabla de métricas

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$betweenness <- betweenness(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)

# Clusterización

clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)

# Puntos de articulación

art_points <- articulation_points(gsuelos)
V(gsuelos)$articulacion <- V(gsuelos)$name %in% V(gsuelos)[art_points]$name

# Tabla de métricas y guardamos

tabla_metricassuelosestresmedio <- data.frame(
  Nodo = V(gsuelos)$name,
  Grado = V(gsuelos)$degree,
  Betweenness = V(gsuelos)$betweenness,
  Closeness = V(gsuelos)$closeness,
  Cluster = V(gsuelos)$cluster,
  Articulacion = V(gsuelos)$articulacion
)

write.table(tabla_metricassuelosestresmedio, 
            file = "tabla_metricassuelosestresmedio.tsv",  
            sep = "\t",  # Delimitador de tabulaciones
            row.names = FALSE,  
            col.names = TRUE, 
            quote = FALSE) 



###################################################################################

# Estrés bajo

# SUELOS 16S

otu_df_suelosbajo<- otu_df_suelos[rownames(otu_df_suelos) %in% c("67","68","80","81","82","86","87","94","95","96","105","106","113","114","115"), ]

#Semillas 16S

otu_df_semillasbajo<-otu_df_semillas[rownames(otu_df_semillas) %in% c ("7","8","19","20","21","28","29","30","34","35","36","40","41","42")]

nombres_suelosbajo <- colnames(otu_df_suelosbajo)
nombres_semillasbajo <- colnames(otu_df_semillasbajo)


# Correlación suelos medio

cor_matrixbajo <- cor(otu_df_suelosbajo, method = "spearman")  

# Función para calcular matriz de p-valores

cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}

# Calculamos matriz de p-valores

p_matrixsuelosbajo<- cor_pmat(otu_df_suelosbajo, method = "spearman")

# Convertimos a vector para aplicar corrección
p_vec_suelosbajo<- as.vector(p_matrixsuelosbajo)

# Aplicamos corrección
p_adj_vec_suelosbajo <- p.adjust(p_vec_suelosbajo, method = "BH")
p_matrixsuelosbajo <- matrix(p_adj_vec_suelosbajo, nrow = nrow(p_matrixsuelosbajo), ncol = ncol(p_matrixsuelosbajo))

# Convertimos las matrices a formato largo

cor_df <- reshape2::melt(cor_matrixbajo)
pval_df <- reshape2::melt(p_matrixsuelosbajo)

# Combinamos en un data frame

cor_pval_df <- cbind(cor_df, pval = pval_df$value)

# Filtramos

cor_pval_fuerte_sig <- subset(cor_pval_df, value > 0.8 & value <0.99 & pval > 0 & pval < 0.05 & Var1 != Var2)

cor_pval_fuerte_sig <- cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ]

cor_pval_fuerte_sig_filtradosuelosbajos <- subset(cor_pval_fuerte_sig, 
                                                  !grepl("^OTU", Var1) & !grepl("^OTU", Var2))



# Creamos grafo

gsuelos <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosuelosbajos, directed = FALSE)

# Visualizamos 
plot(gsuelos, vertex.label.cex = 0.7, edge.width = E(gsuelos)$value * 5, main = "Correlaciones Fuertes y Significativas")

# Calculamos métricas

transitivity(gsuelos, type = "global")
transitivity(gsuelos, type = "local")
btw <- betweenness(gsuelos, directed = FALSE)
V(gsuelos)$betweenness <- btw
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
edge_density(gsuelos)

# Clusterizacion

transitivity(gsuelos, type = "global")

# Grado

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)

#Puntos de articulacion 

articulation_points(gsuelos)


# Layout
lay <- layout_with_fr(gsuelos, niter = 1000)

# Paleta para los clusters
n_clusters <- length(unique(V(gsuelos)$cluster))
colores <- brewer.pal(min(n_clusters, 12), "Paired")

# Graficamos grafo

plot(gsuelos,
     layout = lay,
     vertex.color = colores[V(gsuelos)$cluster],
     vertex.size = log(V(gsuelos)$degree + 1) * 5,
     vertex.label.cex = 0.6,
     vertex.label.color = "black", 
     edge.width = E(gsuelos)$value * 2,
     main = "Red de correlación suelos 16s estrés bajo")

# Tabla de métricas

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$betweenness <- betweenness(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)

# Clusterización

clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)

# Puntos de articulación

art_points <- articulation_points(gsuelos)
V(gsuelos)$articulacion <- V(gsuelos)$name %in% V(gsuelos)[art_points]$name

# Creamos data frame y guardamos

tabla_metricassuelosestresbajo <- data.frame(
  Nodo = V(gsuelos)$name,
  Grado = V(gsuelos)$degree,
  Betweenness = V(gsuelos)$betweenness,
  Closeness = V(gsuelos)$closeness,
  Cluster = V(gsuelos)$cluster,
  Articulacion = V(gsuelos)$articulacion
)

write.table(tabla_metricassuelosestresbajo, 
            file = "tabla_metricassuelosestresbajo.tsv",  
            sep = "\t",  
            row.names = FALSE,  
            col.names = TRUE, 
            quote = FALSE) 


#################################################################################

# Semillas 16S estres alto 

# Correlación

cor_matrixsemillasalto <- cor(otu_df_semillasalto, method = "spearman") 

# Función para calcular matriz de p-valores

cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}

# Calculamos p valores

p_matrixsemillasalto <- cor_pmat(otu_df_semillasalto, method = "spearman")

# Convertimos a vector para aplicar corrección

p_vec_semillasalto <- as.vector(p_matrixsemillasalto)

# Aplicamos corrección

p_adj_vec_semillasalto <- p.adjust(p_vec_semillasalto, method = "BH")
p_matrixsemillasalto <- matrix(p_adj_vec_semillasalto, nrow = nrow(p_matrixsemillasalto), ncol = ncol(p_matrixsemillasalto))



# Convertimos las matrices a formato largo

cor_dfsemillasalto <- reshape2::melt(cor_matrixsemillasalto)
pval_dfsemillasalto <- reshape2::melt(p_matrixsemillasalto)


# Combinamos en un solo data frame
cor_pval_dfsemillasalto <- cbind(cor_dfsemillasalto, pval = pval_dfsemillasalto$value)

# Filtramos

cor_pval_fuerte_sigsemillasalto <- subset(cor_pval_dfsemillasalto, value > 0.7 & value <0.99  & pval > 0 & pval < 0.05 & Var1 != Var2)

cor_pval_fuerte_sigsemillasalto <- cor_pval_fuerte_sigsemillasalto[order(-cor_pval_fuerte_sigsemillasalto$value), ]

cor_pval_fuerte_sig_filtradosemillasalto <- subset(cor_pval_fuerte_sigsemillasalto, 
                                                   !grepl("^OTU", Var1) & !grepl("^OTU", Var2))

# Graficamos

gsemillasalto <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosemillasalto, directed = FALSE)

# Visualizamos

plot(gsemillasalto, vertex.label.cex = 0.7, edge.width = E(gsemillasalto)$value * 5, main = "Correlaciones Fuertes y Significativas")

# Métricas

transitivity(gsemillasalto, type = "global")
transitivity(gsemillasalto, type = "local")
btw <- betweenness(gsemillasalto, directed = FALSE)
V(gsemillasalto)$betweenness <- btw
clusters_btwalto <- cluster_edge_betweenness(gsemillasalto)
V(gsemillasalto)$cluster <- membership(clusters_btwalto)
edge_density(gsemillasalto)
V(gsemillasalto)$degree <- degree(gsemillasalto)
V(gsemillasalto)$closeness <- closeness(gsemillasalto)
articulation_points(gsemillasalto)


# Layout

lay <- layout_with_fr(gsemillasalto, niter = 1000)

# Paleta para los clusters

n_clusterssemillasalto <- length(unique(V(gsemillasalto)$cluster))
coloressemillasalto <- colorRampPalette(brewer.pal(12, "Set3"))(15)

# Graficamos 

plot(gsemillasalto,
     layout = lay,
     vertex.color = coloressemillasalto[V(gsemillasalto)$cluster],
     vertex.size = log(V(gsemillasalto)$degree + 1) * 5,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = E(gsemillasalto)$value * 2,  # antes era *5
     main = "Correlación Semillas 16s Estrés Alto")


# Tabla de métricas 

V(gsemillasalto)$degree <- degree(gsemillasalto)
V(gsemillasalto)$betweenness <- betweenness(gsemillasalto)
V(gsemillasalto)$closeness <- closeness(gsemillasalto)

# Clusterización

clusters_btwalto <- cluster_edge_betweenness(gsemillasalto)
V(gsemillasalto)$cluster <- membership(clusters_btwalto)

# Puntos de articulación

art_points <- articulation_points(gsemillasalto)
V(gsemillasalto)$articulacion <- V(gsemillasalto)$name %in% V(gsemillasalto)[art_points]$name

# Creamos data frame y guardamos

tabla_metricassemillasalto <- data.frame(
  Nodo = V(gsemillasalto)$name,
  Grado = V(gsemillasalto)$degree,
  Betweenness = V(gsemillasalto)$betweenness,
  Closeness = V(gsemillasalto)$closeness,
  Cluster = V(gsemillasalto)$cluster,
  Articulacion = V(gsemillasalto)$articulacion
)

write.table(tabla_metricassemillasalto, 
            file = "tabla_metricassemillasestresalto.tsv",  
            sep = "\t",  
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE) 

####################

# Semillas 16S estrés medio

# Matriz de correlación

cor_matrixsemillasmedio <- cor(otu_df_semillasmedio, method = "spearman") 

# Función para calcular matriz de p-valores

cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}

# Calculamos p valores
p_matrixsemillasmedio <- cor_pmat(otu_df_semillasmedio, method = "spearman")

# Convertimos a vector para aplicar corrección

p_vec_semillasmedio <- as.vector(p_matrixsemillasmedio)

# Aplicamos corrección

p_adj_vec_semillasmedio <- p.adjust(p_vec_semillasmedio, method = "BH")
p_matrixsemillasmedio <- matrix(p_adj_vec_semillasmedio, nrow = nrow(p_matrixsemillasmedio), ncol = ncol(p_matrixsemillasmedio))



# Convertimos las matrices a formato largo

cor_dfsemillasmedio <- melt(cor_matrixsemillasmedio)
pval_dfsemillasmedio <- melt(p_matrixsemillasmedio)

# Combinamos en un data frame

cor_pval_dfsemillasmedio <- cbind(cor_dfsemillasmedio, pval = pval_dfsemillasmedio$value)

# Filtramos

cor_pval_fuerte_sigsemillasmedio <- subset(cor_pval_dfsemillasmedio, value > 0.7 & value <0.99  & pval > 0 & pval < 0.05 & Var1 != Var2)

cor_pval_fuerte_sigsemillasmedio <- cor_pval_fuerte_sigsemillasmedio[order(-cor_pval_fuerte_sigsemillasmedio$value), ]

cor_pval_fuerte_sig_filtradosemillasmedio <- subset(cor_pval_fuerte_sigsemillasmedio, 
                                                   !grepl("^OTU", Var1) & !grepl("^OTU", Var2))

# Creamos grafo

gsemillasmedio <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosemillasmedio, directed = FALSE)

# Visualizamos

plot(gsemillasmedio, vertex.label.cex = 0.7, edge.width = E(g)$value * 5, main = "Correlaciones Fuertes y Significativas")

# Calculamos métricas

transitivity(gsemillasmedio, type = "global")
transitivity(gsemillasmedio, type = "local")
btw <- betweenness(gsemillasmedio, directed = FALSE)
V(gsemillasmedio)$betweenness <- btw
clusters_btwalto <- cluster_edge_betweenness(gsemillasmedio)
V(gsemillasmedio)$cluster <- membership(clusters_btwalto)
edge_density(gsemillasmedio)

#Clusterizacion

transitivity(gsemillas, type = "global")

#Grado

V(gsemillasmedio)$degree <- degree(gsemillasmedio)
V(gsemillasmedio)$closeness <- closeness(gsemillasmedio)

# Puntos de articulación 

articulation_points(gsemillasmedio)


# Layout 

lay <- layout_with_fr(gsemillasmedio, niter = 1000)

# Paleta para los clusters

n_clusterssemillasmedio <- length(unique(V(gsemillasmedio)$cluster))
coloressemillasmedio <- colorRampPalette(brewer.pal(12, "Set3"))(15)

# Grafo

plot(gsemillasmedio,
     layout = lay,
     vertex.color = coloressemillasmedio[V(gsemillasmedio)$cluster],
     vertex.size = log(V(gsemillasmedio)$degree + 1) * 5,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = E(gsemillasmedio)$value * 2,  # antes era *5
     main = "Correlación Semillas 16s Estrés Medio")

# Creamos data frame y guardamos

tabla_metricassemillasmedio <- data.frame(
  Nodo = V(gsemillasmedio)$name,
  Grado = V(gsemillasmedio)$degree,
  Betweenness = V(gsemillasmedio)$betweenness,
  Closeness = V(gsemillasmedio)$closeness,
  Cluster = V(gsemillasmedio)$cluster,
  Articulacion = V(gsemillasmedio)$articulacion
)

write.table(tabla_metricassemillasmedio, 
            file = "tabla_metricassemillasestresmedio.tsv",  
            sep = "\t",  
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE) 


#############################

# Semillas estrés bajo

# Correlación

cor_matrixsemillasbajo <- cor(otu_df_semillasbajo, method = "spearman")  

# Función para calcular matriz de p-valores

cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}


# P valores
p_matrixsemillasbajo <- cor_pmat(otu_df_semillasbajo, method = "spearman")

# Convertimos a vector para aplicar corrección

p_vec_semillasbajo <- as.vector(p_matrixsemillasbajo)

# Aplicamos corrección

p_adj_vec_semillasbajo <- p.adjust(p_vec_semillasbajo, method = "BH")
p_matrixsemillasbajo <- matrix(p_adj_vec_semillasbajo, nrow = nrow(p_matrixsemillasbajo), ncol = ncol(p_matrixsemillasbajo))

# Convertimos las matrices a formato largo

cor_dfsemillasbajo <- melt(cor_matrixsemillasbajo)
pval_dfsemillasbajo <- melt(p_matrixsemillasbajo)

# Combinamos en un solo data frame

cor_pval_dfsemillasbajo <- cbind(cor_dfsemillasbajo, pval = pval_dfsemillasbajo$value)

# Filtramos 

cor_pval_fuerte_sigsemillasbajo <- subset(cor_pval_dfsemillasbajo, value > 0.7 & value <0.95  & pval > 0 & pval < 0.05 & Var1 != Var2)

cor_pval_fuerte_sigsemillasbajo <- cor_pval_fuerte_sigsemillasbajo[order(-cor_pval_fuerte_sigsemillasbajo$value), ]

cor_pval_fuerte_sig_filtradosemillasbajo <- subset(cor_pval_fuerte_sigsemillasbajo, 
                                                    !grepl("^OTU", Var1) & !grepl("^OTU", Var2))

# Creamos grafo 

gsemillasbajo <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosemillasbajo, directed = FALSE)

# Visualizamos
plot(gsemillasbajo, vertex.label.cex = 0.7, edge.width = E(g)$value * 5, main = "Correlaciones Fuertes y Significativas")

# Métricas

transitivity(gsemillasbajo, type = "global")
transitivity(gsemillasbajo, type = "local")
btw <- betweenness(gsemillasbajo, directed = FALSE)
V(gsemillasbajo)$betweenness <- btw
clusters_btwalto <- cluster_edge_betweenness(gsemillasbajo)
V(gsemillasbajo)$cluster <- membership(clusters_btwalto)
edge_density(gsemillasbajo)

#Clusterizacion

transitivity(gsemillas, type = "global")

# Grado

V(gsemillasbajo)$degree <- degree(gsemillasbajo)
V(gsemillasbajo)$closeness <- closeness(gsemillasbajo)

# Puntos de articulación

articulation_points(gsemillasbajo)


# Layout 

lay <- layout_with_fr(gsemillasbajo, niter = 1000)

# Paleta para los clusters

n_clusterssemillasbajo <- length(unique(V(gsemillasbajo)$cluster))
coloressemillasbajo <- colorRampPalette(brewer.pal(12, "Set3"))(15)

# Graficamos

plot(gsemillasbajo,
     layout = lay,
     vertex.color = coloressemillasbajo[V(gsemillasbajo)$cluster],
     vertex.size = log(V(gsemillasbajo)$degree + 1) * 5,
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = E(gsemillasbajo)$value * 2, 
     main = "Correlación Semillas 16s Estrés Bajo")




# Métricas

V(gsemillasbajo)$degree <- degree(gsemillasbajo)
V(gsemillasbajo)$betweenness <- betweenness(gsemillasbajo)
V(gsemillasbajo)$closeness <- closeness(gsemillasbajo)

# Clusterización
clusters_btwalto <- cluster_edge_betweenness(gsemillasbajo)
V(gsemillasbajo)$cluster <- membership(clusters_btwalto)

# Puntos de articulación

art_points <- articulation_points(gsemillasbajo)
V(gsemillasbajo)$articulacion <- V(gsemillasbajo)$name %in% V(gsemillasbajo)[art_points]$name

# Creamos data frame y guardamos

tabla_metricassemillasbajo<- data.frame(
  Nodo = V(gsemillasbajo)$name,
  Grado = V(gsemillasbajo)$degree,
  Betweenness = V(gsemillasbajo)$betweenness,
  Closeness = V(gsemillasbajo)$closeness,
  Cluster = V(gsemillasbajo)$cluster,
  Articulacion = V(gsemillasbajo)$articulacion
)

write.table(tabla_metricassemillasbajo, 
            file = "tabla_metricassemillasestresbajo.tsv",  
            sep = "\t", 
            row.names = FALSE,  
            col.names = TRUE, 
            quote = FALSE) 




#########################################################################################

# Suelos ITS

# Cargamos los datos

otu_table_suelosITS <- otu_table(phylobject_otu_ITSsu)

# Convertimos a data frame para modificar

otu_df_suelosITS <- as.data.frame(otu_table_suelosITS)
otu_taxonomy_suelosITS <- tax_table(phylobject_otu_ITSsu)
otu_taxonomy_suelosITS <- as.data.frame(otu_taxonomy_suelosITS)
otu_df_suelosITS <- as.data.frame(t(otu_table_suelosITS))
otu_idsueloITS <- colnames(otu_df_suelosITS)
generoITS <- otu_taxonomy_suelosITS$Genus[match(otu_idsueloITS, rownames(otu_taxonomy_suelosITS))]
colnames(otu_df_suelosITS) <- generoITS

# Eliminamos NAs

colnames(otu_df_suelosITS) <- ifelse(is.na(generoITS), otu_idsueloITS, generoITS)


#Suelos estrés alto ITS

otu_df_suelosaltoITS<- otu_df_suelosITS[rownames(otu_df_suelosITS) %in% c("83","84","85","91","92","93","97","98","99","100","101","102","103","104","107", "108", "109","110","111","112","116","120","121","122","123","124","125"), ]

nombres_suelosaltoITS <- colnames(otu_df_suelosaltoITS)


# Correlación

cor_matrixalto <- cor(otu_df_suelosaltoITS, method = "spearman")  

# Calculamos p valores

cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}

p_matrixsuelosalto <- cor_pmat(otu_df_suelosaltoITS, method = "spearman")

# Convertimos a vector para aplicar corrección

p_vec_suelosalto <- as.vector(p_matrixsuelosalto)

# Aplicamos corrección

p_adj_vec_suelosalto <- p.adjust(p_vec_suelosalto, method = "BH")
p_matrixsuelosalto <- matrix(p_adj_vec_suelosalto, nrow = nrow(p_matrixsuelosalto), ncol = ncol(p_matrixsuelosalto))

# Convertimos a formato largo 

cor_df <- reshape2::melt(cor_matrixalto)
pval_df <- reshape2::melt(p_matrixsuelosalto)

# Combinamos

cor_pval_df <- cbind(cor_df, pval = pval_df$value)

# Filtramos

cor_pval_fuerte_sig <- subset(cor_pval_df, value > 0.8 & value <0.89 & pval > 0 & pval < 0.05 & Var1 != Var2)

cor_pval_fuerte_sig <- cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ]



# Cogemos las 200 primeras

cor_pval_fuerte_sig_top200 <- head(cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ], 200)

cor_pval_fuerte_sig_filtradosuelosITS <- subset(cor_pval_fuerte_sig_top200,
                                                !grepl("^OTU", Var1) & 
                                                  !grepl("^OTU", Var2) &
                                                  !grepl("_Incertae_", Var1) & 
                                                  !grepl("_Incertae_", Var2))



# Creamos el grafo 

gsuelos <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosuelosITS, directed = FALSE)

# Visualizamos

plot(gsuelos, vertex.label.cex = 0.7, edge.width = E(gsuelos)$value * 5, main = "Correlaciones Fuertes y Significativas")
plot(gsuelos,
     vertex.label = gsub("^_?[a-z]__?", "", V(gsuelos)$name),
     vertex.label.cex = 0.7,
     edge.width = E(gsuelos)$value * 5,
     main = "Correlaciones Fuertes y Significativas")


# Metricas

transitivity(gsuelos, type = "global")
transitivity(gsuelos, type = "local")
btw <- betweenness(gsuelos, directed = FALSE)
V(gsuelos)$betweenness <- btw
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
edge_density(gsuelos)

# Clusterizacion

transitivity(gsuelos, type = "global")

# Grado

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)

# Puntos de articulacion

articulation_points(gsuelos)


# Layout
lay <- layout_with_fr(gsuelos, niter = 1000)

# Paleta para los clusters

n_clusters <- length(unique(V(gsuelos)$cluster))
colores <- colorRampPalette(brewer.pal(12, "Set3"))(19)

# Grafo

plot(gsuelos,
     layout = lay,
     vertex.color = colores[V(gsuelos)$cluster],
     vertex.size = log(V(gsuelos)$degree + 1) * 5,
     vertex.label = gsub("^_?[a-z]__+", "", V(gsuelos)$name),  # quitamos los "_"
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = E(gsuelos)$value * 2,
     main = "Red de correlación Suelos ITS estrés alto")

# Tabla de métricas

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$betweenness <- betweenness(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
art_points <- articulation_points(gsuelos)
V(gsuelos)$articulacion <- V(gsuelos)$name %in% V(gsuelos)[art_points]$name

# Creamos data frame y guardamos

tabla_metricassuelosestresaltoITS <- data.frame(
  Nodo = V(gsuelos)$name,
  Grado = V(gsuelos)$degree,
  Betweenness = V(gsuelos)$betweenness,
  Closeness = V(gsuelos)$closeness,
  Cluster = V(gsuelos)$cluster,
  Articulacion = V(gsuelos)$articulacion
)

write.table(tabla_metricassuelosestresaltoITS, 
            file = "tabla_metricassuelosestresaltoITS.tsv", 
            sep = "\t", 
            row.names = FALSE,  
            col.names = TRUE,  
            quote = FALSE) 



#########################################################

# Suelos estrés medio

otu_df_suelosmedioITS<- otu_df_suelosITS[rownames(otu_df_suelosITS) %in% c("70","71","72","73","74","75","76","77","78","79","88","89","90","117","118", "119"), ]

nombres_suelosmedioITS <- colnames(otu_df_suelosmedioITS)

cor_matrixmedio <- cor(otu_df_suelosmedioITS, method = "spearman")  

# Calculamos p valores

cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}

# p valores

p_matrixsuelosmedio <- cor_pmat(otu_df_suelosmedioITS, method = "spearman")

# Convertimos a vector para aplicar corrección

p_vec_suelosmedio <- as.vector(p_matrixsuelosmedio)

# Aplicamos corrección

p_adj_vec_suelosmedio <- p.adjust(p_vec_suelosmedio, method = "BH")
p_matrixsuelosmedio <- matrix(p_adj_vec_suelosmedio, nrow = nrow(p_matrixsuelosmedio), ncol = ncol(p_matrixsuelosmedio))

# Convertimos 

cor_df <- reshape2::melt(cor_matrixmedio)
pval_df <- reshape2::melt(p_matrixsuelosmedio)

# Combinamos

cor_pval_df <- cbind(cor_df, pval = pval_df$value)

# Filtramos

cor_pval_fuerte_sig <- subset(cor_pval_df, value > 0.8 & value <0.89 & pval > 0 & pval < 0.05 & Var1 != Var2)

cor_pval_fuerte_sig <- cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ]



# 200 primeras

cor_pval_fuerte_sig_top200 <- head(cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ], 200)

cor_pval_fuerte_sig_filtradosuelosmedioITS <- subset(cor_pval_fuerte_sig_top200,
                                                !grepl("^OTU", Var1) & 
                                                  !grepl("^OTU", Var2) &
                                                  !grepl("_Incertae_", Var1) & 
                                                  !grepl("_Incertae_", Var2))


# Graficamos

gsuelos <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosuelosmedioITS, directed = FALSE)

# Visualizamos 

plot(gsuelos, vertex.label.cex = 0.7, edge.width = E(gsuelos)$value * 5, main = "Correlaciones Fuertes y Significativas")
plot(gsuelos,
     vertex.label = gsub("^_?[a-z]__?", "", V(gsuelos)$name),
     vertex.label.cex = 0.7,
     edge.width = E(gsuelos)$value * 5,
     main = "Correlaciones Fuertes y Significativas")

# Metricas

transitivity(gsuelos, type = "global")
transitivity(gsuelos, type = "local")
btw <- betweenness(gsuelos, directed = FALSE)
V(gsuelos)$betweenness <- btw
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
edge_density(gsuelos)

#Clusterizacion

transitivity(gsuelos, type = "global")

# Grado

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)

#Puntos de articulacion 

articulation_points(gsuelos)


# Layout
lay <- layout_with_fr(gsuelos, niter = 1000)

# Paleta para los clusters

n_clusters <- length(unique(V(gsuelos)$cluster))
colores <- colorRampPalette(brewer.pal(12, "Set3"))(19)

# Grafico final

plot(gsuelos,
     layout = lay,
     vertex.color = colores[V(gsuelos)$cluster],
     vertex.size = log(V(gsuelos)$degree + 1) * 5,
     vertex.label = gsub("^_?[a-z]__+", "", V(gsuelos)$name),  # eliminamos los "_"
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = E(gsuelos)$value * 2,
     main = "Red de correlación Suelos ITS estrés alto")

# Metricas

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$betweenness <- betweenness(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
art_points <- articulation_points(gsuelos)
V(gsuelos)$articulacion <- V(gsuelos)$name %in% V(gsuelos)[art_points]$name

# Creamos data frame y guardamos
tabla_metricassuelosestresmedioITS <- data.frame(
  Nodo = V(gsuelos)$name,
  Grado = V(gsuelos)$degree,
  Betweenness = V(gsuelos)$betweenness,
  Closeness = V(gsuelos)$closeness,
  Cluster = V(gsuelos)$cluster,
  Articulacion = V(gsuelos)$articulacion
)

write.table(tabla_metricassuelosestresmedioITS, 
            file = "tabla_metricassuelosestresmedioITS.tsv",  
            sep = "\t",  
            row.names = FALSE,  
            col.names = TRUE,  
            quote = FALSE) 

############################################################################

# Estrés bajo ITS suelos

otu_df_suelosbajoITS<- otu_df_suelosITS[rownames(otu_df_suelosITS) %in% c("67","68","80","81","82","86","87","94","95","96","105","106","113","114","115"), ]

nombres_suelosbajoITS <- colnames(otu_df_suelosbajoITS)

# Correlación 

cor_matrixbajo <- cor(otu_df_suelosbajoITS, method = "spearman")  

# p- valores
cor_pmat <- function(x, method = "spearman") {
  n <- ncol(x)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- colnames(x)
  rownames(p.mat) <- colnames(x)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(x[, i], x[, j], method = method)
      p.mat[i, j] <- test$p.value
      p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}


# Calculamos 

p_matrixsuelosbajoITS <- cor_pmat(otu_df_suelosbajoITS, method = "spearman")

# Convertimos a vector para aplicar corrección

p_vec_suelosbajo <- as.vector(p_matrixsuelosbajo)

# Aplicamos corrección

p_adj_vec_suelosbajo <- p.adjust(p_vec_suelosbajo, method = "BH")
p_matrixsuelosbajo <- matrix(p_adj_vec_suelosbajo, nrow = nrow(p_matrixsuelosbajo), ncol = ncol(p_matrixsuelosbajo))

# Convertimos

cor_df <- reshape2::melt(cor_matrixbajo)
pval_df <- reshape2::melt(p_matrixsuelosbajoITS)

# Combinamos

cor_pval_df <- cbind(cor_df, pval = pval_df$value)

# Filtramos

cor_pval_fuerte_sig <- subset(cor_pval_df, value > 0.8 & value <0.89 & pval > 0 & pval < 0.05 & Var1 != Var2)

cor_pval_fuerte_sig <- cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ]



# 200 primeras
cor_pval_fuerte_sig_top200 <- head(cor_pval_fuerte_sig[order(-cor_pval_fuerte_sig$value), ], 200)

cor_pval_fuerte_sig_filtradosuelosmedioITS <- subset(cor_pval_fuerte_sig_top200,
                                                     !grepl("^OTU", Var1) & 
                                                       !grepl("^OTU", Var2) &
                                                       !grepl("_Incertae_", Var1) & 
                                                       !grepl("_Incertae_", Var2))

# Grafo

gsuelos <- graph_from_data_frame(cor_pval_fuerte_sig_filtradosuelosmedioITS, directed = FALSE)

# Visualizamos

plot(gsuelos, vertex.label.cex = 0.7, edge.width = E(gsuelos)$value * 5, main = "Correlaciones Fuertes y Significativas")
plot(gsuelos,
     vertex.label = gsub("^_?[a-z]__?", "", V(gsuelos)$name),
     vertex.label.cex = 0.7,
     edge.width = E(gsuelos)$value * 5,
     main = "Correlaciones Fuertes y Significativas")


# Metricas

transitivity(gsuelos, type = "global")
transitivity(gsuelos, type = "local")
btw <- betweenness(gsuelos, directed = FALSE)
V(gsuelos)$betweenness <- btw
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
edge_density(gsuelos)
transitivity(gsuelos, type = "global")
V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)
articulation_points(gsuelos)


# Layout

lay <- layout_with_fr(gsuelos, niter = 1000)

# Paleta para los clusters

n_clusters <- length(unique(V(gsuelos)$cluster))
colores <- colorRampPalette(brewer.pal(12, "Set3"))(19)

# Grafo final

plot(gsuelos,
     layout = lay,
     vertex.color = colores[V(gsuelos)$cluster],
     vertex.size = log(V(gsuelos)$degree + 1) * 5,
     vertex.label = gsub("^_?[a-z]__+", "", V(gsuelos)$name),  # eliminamos los "_"
     vertex.label.cex = 0.6,
     vertex.label.color = "black",
     edge.width = E(gsuelos)$value * 2,
     main = "Red de correlación Suelos ITS estrés bajo")

# Metricas

V(gsuelos)$degree <- degree(gsuelos)
V(gsuelos)$betweenness <- betweenness(gsuelos)
V(gsuelos)$closeness <- closeness(gsuelos)
clusters_btw <- cluster_edge_betweenness(gsuelos)
V(gsuelos)$cluster <- membership(clusters_btw)
art_points <- articulation_points(gsuelos)
V(gsuelos)$articulacion <- V(gsuelos)$name %in% V(gsuelos)[art_points]$name

# Creamos data frame y guardamos

tabla_metricassuelosestresbajoITS <- data.frame(
  Nodo = V(gsuelos)$name,
  Grado = V(gsuelos)$degree,
  Betweenness = V(gsuelos)$betweenness,
  Closeness = V(gsuelos)$closeness,
  Cluster = V(gsuelos)$cluster,
  Articulacion = V(gsuelos)$articulacion
)

write.table(tabla_metricassuelosestresbajoITS, 
            file = "tabla_metricassuelosestresbajoITS.tsv", 
            sep = "\t", 
            row.names = FALSE,  
            col.names = TRUE,  
            quote = FALSE)

