#PROYECTO INTEGRADOR BT1013

#A01752953 Frida Cano Falcón  
#A01745037 Gala Flores García  
#A01750150 Hortencia Alejandra Ramírez Vázquez  
#A01751655 Gabriela Cortés Olvera  
#A01751568 Harumi Cristal Manzano Yáñez  
#A01746530 Diana Paola López Espinosa   


#Lectura de conjunto de datos

library(GEOquery)
library(limma)
knitr::opts_chunk$set(echo = TRUE)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

# Lee el conjunto de datos
gset <- getGEO("GSE161097", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Obtiene los valores de expresión
ex <- exprs(gset)

probes <- getEAWP(gset)$probes
probes$`Gene ID`[1]

#View(gset)
as.factor(gset$source_name_ch1)

primary_colon_cancer_mucinous <- exprs(gset[,gset$source_name_ch1 =="primary colon cancer" & gset$histology =="mucinous"])
primary_colon_cancer_mucinous <- data.frame(primary_colon_cancer_mucinous)

primary_colon_cancer_nonmucinous <- exprs(gset[,gset$source_name_ch1 =="primary colon cancer" & gset$histology =="non mucinous"])
primary_colon_cancer_nonmucinous <- data.frame(primary_colon_cancer_nonmucinous)


peritoneal_metastases_mucinous <- exprs(gset[,gset$source_name_ch1 =="peritoneal metastases" & gset$histology =="mucinous"])
peritoneal_metastases_mucinous <- data.frame(peritoneal_metastases_mucinous)

peritoneal_metastases_nonmucinous <- exprs(gset[,gset$source_name_ch1 =="peritoneal metastases" & gset$histology =="non mucinous"])
peritoneal_metastases_nonmucinous <- data.frame(peritoneal_metastases_nonmucinous)


#Agrupación de los 4 dataframes

microarray=cbind(primary_colon_cancer_mucinous,primary_colon_cancer_nonmucinous,peritoneal_metastases_mucinous,peritoneal_metastases_nonmucinous)

# 3. Normalizacion
raw_means = apply(microarray,2,mean,trim=0.02) 

microarray_norm = sweep(microarray, 2, raw_means, "/") * 100

# 4. Medias
primary_colon_cancer_mucinous_mean = rowMeans(microarray_norm[,1:7])
primary_colon_cancer_nonmucinous_mean = rowMeans(microarray_norm[,8:15])
peritoneal_metastases_mucinous_mean = rowMeans(microarray_norm[,16:22])
peritoneal_metastases_nonmucinous_mean = rowMeans(microarray_norm[,23:30])

microarray_means = data.frame(primary_colon_cancer_mucinous_mean,primary_colon_cancer_nonmucinous_mean,peritoneal_metastases_mucinous_mean,peritoneal_metastases_nonmucinous_mean)

# 5. Proporciones (mucinous/nonmucinous)

primary_colon_cancer_ratios = microarray_means$primary_colon_cancer_mucinous_mean / microarray_means$primary_colon_cancer_nonmucinous_mean

peritoneal_metastases_ratios = microarray_means$peritoneal_metastases_mucinous_mean / microarray_means$peritoneal_metastases_nonmucinous_mean

microarray_ratios = data.frame(primary_colon_cancer_ratios, peritoneal_metastases_ratios)

row.names(microarray_ratios) = row.names(microarray)

# 6. Cambio a log2

microarray_norm = log2(microarray_norm)
microarray_means = log2(microarray_means)
microarray_ratios = log2(microarray_ratios)

# 7. t-test

get_pvalue <- function(values, idx1, idx2) {
  return(t.test(values[idx1], values[idx2])$p.value)
  
}

primary_colon_cancer_p = apply(microarray_norm, 1, get_pvalue, 1:7, 8:15)
peritoneal_metastases_p = apply(microarray_norm, 1, get_pvalue, 16:22, 23:30)

# Seleccion de genes con p values menor a 0.05

filtered_primary_colon_cancer_p = primary_colon_cancer_p[primary_colon_cancer_p < 0.05]
filtered_peritoneal_metastases_p = peritoneal_metastases_p[peritoneal_metastases_p < 0.05]

#********************************************************************************

# Seleccionar genes con menor p_value y con un tamaño del efecto "grande"

  #Selección de genes de tamaño del efecto grande

filtred_primary_colon_cancer_MA = primary_colon_cancer_p[primary_colon_cancer_p < 0.05 & microarray_ratios$primary_colon_cancer_ratios > 0.25]
filtred_primary_colon_cancer_ME = primary_colon_cancer_p[primary_colon_cancer_p < 0.05 & microarray_ratios$primary_colon_cancer_ratios < -0.25]

filtred_peritoneal_metastases_MA = primary_colon_cancer_p[peritoneal_metastases_p < 0.05 & microarray_ratios$peritoneal_metastases_ratios > 0.25]
filtred_peritoneal_metastases_ME = primary_colon_cancer_p[peritoneal_metastases_p < 0.05 & microarray_ratios$peritoneal_metastases_ratios < -0.25]

#Ordenar datos
datos_pvalues_primary_colon_MA <- sort(filtred_primary_colon_cancer_MA)
datos_pvalues_primary_colon_ME <- sort(filtred_primary_colon_cancer_ME)

datos_pvalues_peritoneal_metastases_MA <- sort(filtred_peritoneal_metastases_MA)
datos_pvalues_peritoneal_metastases_ME <- sort(filtred_peritoneal_metastases_ME)

#Selección de los primeros 10 genes con menor pvalue y tamaño de afecto grande

#Primary colon

#Mayor expresión en tejido mucinoso
  
datos_primary_colonTop10_MA = datos_pvalues_primary_colon_MA[1:10]

df_datos_primary_colonTop10_MA<-data.frame(datos_primary_colonTop10_MA)

datos_primary_colonTop10_names_MA=row.names(df_datos_primary_colonTop10_MA)
Top10_genes_primary_MA=probes[datos_primary_colonTop10_names_MA,3] #Top
df_Top10_genes_primary_MA<-data.frame(Top10_genes_primary_MA)

#Mayor expresión en tejido no mucinoso

datos_primary_colonTop10_ME = datos_pvalues_primary_colon_ME[1:3]

df_datos_primary_colonTop10_ME<-data.frame(datos_primary_colonTop10_ME)

datos_primary_colonTop10_names_ME=row.names(df_datos_primary_colonTop10_ME)
Top10_genes_primary_ME=probes[datos_primary_colonTop10_names_ME,3] #Top
df_Top10_genes_primary_ME<-data.frame(Top10_genes_primary_ME)

#Peritoneal Metastases

#Mayor expresión en tejido mucinoso

datos_peritoneal_Top10_MA = datos_pvalues_peritoneal_metastases_MA[1:10]

df_datos_peritoneal_Top10_MA<-data.frame(datos_peritoneal_Top10_MA)

datos_peritoneal_Top10_names_MA=row.names(df_datos_peritoneal_Top10_MA)
Top10_genes_peritoneal_MA=probes[datos_peritoneal_Top10_names_MA,3] #Top
df_Top10_genes_peritoneal_MA<-data.frame(Top10_genes_peritoneal_MA)

#Mayor expresión en tejido no mucinoso

datos_peritoneal_Top10_ME = datos_pvalues_peritoneal_metastases_ME[1:10]

df_datos_peritoneal_Top10_ME<-data.frame(datos_peritoneal_Top10_ME)

datos_peritoneal_Top10_names_ME=row.names(df_datos_peritoneal_Top10_ME)
Top10_genes_peritoneal_ME=probes[datos_peritoneal_Top10_names_ME,3] #Top
df_Top10_genes_peritoneal_ME<-data.frame(Top10_genes_peritoneal_ME)

  
#*******************************************************************************
# Justificación de la selección de los genes

#Para hacer la selección de los genes diferencialmente expresados se aplican dos 
#métodos para el filtrado: primero seleccionando aquellos con un p-value menor a 
#0.05 y con un tamaño del efecto mayor a 0.25 y menor a -0.25, obteniendo una lista 
#de genes con un nivel de expresión 2 veces mayor en tejidos con presencia de mucina
#y 2 veces mayor con ausencia de mucina para cada tejido analizado.

#Posteriormente los resultados obtenidos se ordenaron de menor a mayor en función del p-value.
#De cada una de las listas se extrajeron los primeros diez genes, obteniendo así dos 
#listas de diez genes con mayor expresión en un tejido mucinoso y con mayor expresión 
#en un tejido con ausencia de mucina, esto para los dos casos analizados: cáncer de colon 
#en etapa inicial y metástasis peritoneal. A estas listas obtenidas se les asignó el ID 
#Gene para lograr una mejor identificación de ellos.  

#Aplicando los filtros anteriormente mencionados en la sección de metodología sobre
#cada uno de nuestros grupos, se redujo significativamente el número de genes; para 
#el grupo cáncer de colon en etapa inicial que contiene originalmente 780 genes se redujo 
#a 13 con presencia de mucina y 3 genes sin mucina, mientras que en el tejido de metástasis 
#peritoneal  se redujo de 780 genes a 12 con mucina y 23 sin mucina.


#*******************************************************************************

#Selección de datos que están en ambas listas
filtered_probes = sort(intersect(names(filtered_primary_colon_cancer_p ), names(filtered_peritoneal_metastases_p)))
microarray_selection = microarray_means[filtered_probes,]
microarray_selection_ratios = microarray_ratios[filtered_probes,]


# Plots

# Mapa de calor y dendograma

medias <- rowMeans(microarray_selection)
devs = apply(microarray_selection, 1, sd)

centered_microarray_selection = sweep(microarray_selection, 1, medias)
centered_microarray_selection = sweep(centered_microarray_selection, 1, devs, "/")

names(centered_microarray_selection) = c("PM", "PNM", "MM", "MNM") #Cambia el nombre de los grupos en graficas
hclustering = hclust(dist(centered_microarray_selection)) #Hace el clusting /dist que tan diferentes son los genes

plot(hclustering, main = "Dendograma", ylab = "Altura")
names(microarray_selection) = c("PM", "PNM", "MM", "MNM")

heatmap(as.matrix(centered_microarray_selection), Colv = NA)

#Primary Colon cancer

# Gráfica de dispersión

plot(microarray_means$primary_colon_cancer_mucinous_mean, 
     microarray_means$primary_colon_cancer_nonmucinous_mean, 
     xlim = c(0,16), ylim = c(0,16),
     xaxt="n", yaxt="n",
     main = "Expresión en cáncer de colon en etapa inicial: mucinoso vs no mucinoso",
     xlab = "Mucinoso (valor de expresión log2)",
     ylab = "No mucinoso (valor de expresión log2)")
axis(1, at=seq(0,16,2))
axis(2, at=seq(0,16,2))

abline(lm(microarray_means$primary_colon_cancer_nonmucinous_mean ~ microarray_means$primary_colon_cancer_mucinous_mean),
       col = "red")

# Gráfica R-I
plot(microarray_means$primary_colon_cancer_mucinous_mean + microarray_means$primary_colon_cancer_nonmucinous_mean, #Aqui es por los logaritmos
     microarray_means$primary_colon_cancer_mucinous_mean - microarray_means$primary_colon_cancer_nonmucinous_mean,
     main = "Gráfica R-I para cáncer de colon en etapa inicial no mucinoso vs mucinoso",
     xlab = "log2(mucinoso * no mucinoso)",
     ylab = "log2(mucinoso / no mucinoso)")

# Gráfica de volcán
colores = rep(1, length(primary_colon_cancer_p))
colores[primary_colon_cancer_p < 0.05 & microarray_ratios$primary_colon_cancer_ratios < -0.25] = "red"
colores[primary_colon_cancer_p < 0.05 & microarray_ratios$primary_colon_cancer_ratios > 0.25] = "blue"

plot(microarray_ratios$primary_colon_cancer_ratios,primary_colon_cancer_p, col = colores,
     log = "y", ylim = rev(range(primary_colon_cancer_p)),  #escala
     main = "Expresión diferencial para cancer de colon en etapa inicial no mucinoso vs mucinoso",
     xlab = "log2: no mucinoso vs mucinoso",
     ylab = "p-value")

#Peritoneal metastases

# Gráfica de dispersión

plot(microarray_means$peritoneal_metastases_mucinous_mean, 
     microarray_means$peritoneal_metastases_nonmucinous_mean, 
     xlim = c(0,16), ylim = c(0,16),
     xaxt="n", yaxt="n",
     main = "Expresión en metástasis peritoneal: mucinoso vs non mucinoso",
     xlab = "Mucinoso (valor de expresión log2)",
     ylab = "No mucinoso (valor de expresión log2")
axis(1, at=seq(0,16,2))
axis(2, at=seq(0,16,2))

abline(lm(microarray_means$peritoneal_metastases_nonmucinous_mean ~ microarray_means$peritoneal_metastases_mucinous_mean),
       col = "red")

# Gráfica R-I
plot(microarray_means$peritoneal_metastases_mucinous_mean + microarray_means$peritoneal_metastases_nonmucinous_mean, #Aqui es por los logaritmos
     microarray_means$peritoneal_metastases_mucinous_mean - microarray_means$peritoneal_metastases_nonmucinous_mean,
     main = "Gráfica R-I para Metástasis Peritoneal no mucinoso vs mucinoso",
     xlab = "log2(mucinoso * no mucinoso)",
     ylab = "log2(mucinoso / no mucinoso)")

# Gráfica de volcán
colores = rep(1, length(peritoneal_metastases_p))
colores[peritoneal_metastases_p < 0.05 & microarray_ratios$peritoneal_metastases_ratios < -0.25] = "red"
colores[peritoneal_metastases_p < 0.05 & microarray_ratios$peritoneal_metastases_ratios > 0.25] = "blue"

plot(microarray_ratios$peritoneal_metastases_ratios,peritoneal_metastases_p, col = colores,
     log = "y", ylim = rev(range(peritoneal_metastases_p)),  #escala
     main = "Expresión diferencial para Metástasis Peritoneal no mucinoso vs mucinoso",
     xlab = "log2: no mucinoso vs mucinoso",
     ylab = "p-value")

