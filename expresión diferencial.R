BiocManager :: install ("rhdf5") #c:
install.packages("devtools") # c:

BiocManager::install("pachterlab/sleuth") #### esta es la buena


####### Expresión diferencial ######

#Carmen Monserrat Anistro Romero 

library("sleuth")

#Esta funci?n permite mapear, a partir de la base de datos de biomaRT
tx2gene <- function(){
  
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # se cambia la base de datos, para  humanos
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}
# ayudará a conectar una base de datos preexistente con la que tu generarás 
# con los datos que tu deseas analizar

t2g <- tx2gene() # la función se asigna a un objeto que se usará más adelante
t2g #hay que actualizar dplyr


#dirección de la carpet donde estan todos los archivos que deseo analizar
base_dir <-"C:/Users/manis/Desktop/abundancias/adundancias"


# voy a comparar las primeras 3 con las últimas 3, entonces se seleccionan
samples <- paste0("sample", c("1","2","3",
                                   "10","11","12"))
samples

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
kal_dirs # dentro del objeto tengo los que voy a comparar con las direcciones de los archivos


#Selected samples

s2c <- data.frame(path=kal_dirs, sample=samples, muestras <- c("control1","control2","control3", "exp1",
                                                              "exp2", "exp3"), stringsAsFactors=FALSE)
s2c
#hice una base de dtos donde tengo la muestra que voy a comparar con su dirección y se les asigno si son o no controles

so <- sleuth_prep(s2c, ~muestras, target_mapping = t2g, extra_bootstrap_summary = TRUE) 
# normalización, los resultados se leen en kallisto
so <- sleuth_fit(so)
#
so <- sleuth_wt(so, which_beta="muestrasexp") ## se hace la comparación                  
sleuth_live(so) # se necesita shiny para visualizar los momentos hasta el momento 

setwd("C:/Users/manis/Desktop/abundancias/adundancias") 
# regresamos a nuestros archivos

resultados<-read.table("test_table.csv",sep=",", header=TRUE)
# test_table sale en los resultados de setwd, see descarga y se sube

significativos<-which(resultados$qval<0.1) 
#se seleccionan los resultados con un q valor <.1

significativos<-resultados[significativos,] 
# al objeto significativos se le asignan los resultados previamente seleccionados

# lo siguiente se hace para saber cuales genes estan sobre o sub expresados, esto se puede saber
# si b es mayor o menor que ceero, por eso se seleccionn los objetos según esa condición y se asignan a
# objetos diferentes, para tener separdos los resultados

upregulated<-which(significativos$b>0)
upregulated<-significativos[upregulated,]
downregulated<-which(significativos$b<0)
downregulated<-significativos[downregulated,]


write.table(upregulated,file="C:/Users/manis/Desktop/abundancias/adundancias",sep="\t")
write.table(downregulated,file="C:/Users/manis/Desktop/abundancias/adundancias",sep="\t")

upregulated$ext_gene
#de los genes selecciona los sobreexpresados, y te d el listado de sus ID