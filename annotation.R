
##pakages:

library(dplyr)
library(usethis)
library(roxygen2)

if ((params$peakAnnotation$isotopeLikelihoodScoreThreshold > 1) || (params$peakAnnotation$isotopeLikelihoodScoreThreshold < 0))
{
  stop("Error: isotopeLikelihoodScoreThreshold in the ProcParams object must be between 0 and 1")
}

if((params$peakAnnotation$ppmMassTolerance>100))
{
  writeLines("Mass tolerances close or bigger than 100 ppms may produce wrong annotations.")
}

if (!require("rMSI2", character.only = TRUE, quietly=TRUE)) {
  devtools::install_github("https://github.com/prafols/rMSI2.git", dependencies = TRUE)
  library(rMSI2)
}

if (!require("SummarizedExperiment", quietly = TRUE)){
  install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install("SummarizedExperiment")
  library(SummarizedExperiment)
}
if (!require("HDF5Array", character.only = TRUE)) {
  BiocManager::install("HDF5Array")
  library(HDF5Array)
}

if (!require("dplyr", character.only = TRUE, quietly=TRUE)) {
  install.packages("dplyr")
  library(dplyr)
}

#setup

data_path<- "/home/nmolto/Desktop/rMSI2_Nadia/Dataset"
rMSIprocPeakMatrix<- rMSI2::LoadPeakMatrix(file.path(data_path,"tiroides_roger.pkmat"))

AdductAnnotation <- function(rMSIPeacMatrix, AdductList, tolerance) {}

params<-rMSI2::ProcessingParameters()

params$peakAnnotation$adductElementsTable <- as.data.frame()

MetaboCoreUtils::adducts(polarity="negative")

## rMSI2 aductes negatius

negative_adducts<- data.frame(
  name = c("Cl-", "Br-", "-H"),
  mass = c(-34.969402, -78.918885, -1.007276),
  priority = c(0, 0, 0))

negative_adducts<-rbind(addduct_tables_neg[[2]], addduct_tables_neg[[3]])
params$peakAnnotation$adductElementsTable <- negative_adducts
params$peakAnnotation$ppmMassTolerance<-14
aductes <- rMSI2::peakAnnotation(rMSIprocPeakMatrix, params=params)
params$peakAnnotation$isotopeLikelihoodScoreThreshold<- 0.6


## MSIAnnotator

#settings
params<-rMSI2::ProcessingParameters()
params$peakAnnotation$ppmMassTolerance<-15
NH4<- data.frame(name="+NH4", mass=18.033823, priority=0)
params$peakAnnotation$adductElementsTable<- rbind(params$peakAnnotation$adductElementsTable,NH4)

adducts_C_df <- params$peakAnnotation$adductElementsTable

aductes <- rMSI2::peakAnnotation(rMSIprocPeakMatrix, params=params)

#dataframe format from rMSI2 output

process_peak <- function(peak_id) {
  adducts_A <- aductes$A[aductes$A$Adduct1Index == peak_id | aductes$A$Adduct2Index == peak_id,]
  if (nrow(adducts_A) > 0) {
    adducts_A$Level <- "A"
  }
  adducts_B <- aductes$B[aductes$B$Adduct1Index == peak_id | aductes$B$Adduct2Index == peak_id,]
  if (nrow(adducts_B) > 0) {
    adducts_B$Level <- "B"
  }
  monoisotopic_results <- aductes$C[aductes$C$MonoisotopicIndex == peak_id,]
  if (nrow(monoisotopic_results) > 0) {
    monoisotopic_results$Level <- "C"
  }
  combined_results <- plyr::rbind.fill(adducts_A, adducts_B, monoisotopic_results)
  return(combined_results)
}

annotation_mat <- lapply(seq_along(rMSIprocPeakMatrix$mass), process_peak)
annotation_df <- do.call(rbind, annotation_mat)
annotation_df <- annotation_df[!duplicated(annotation_df), ]

#separate rows for each adduct and include the monoisotopic mass
expand_adducts <- function(df, adducts_C_df) {
  df_expanded <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    if (row$Level %in% c("A", "B")) {
      adduct_pairs <- unlist(strsplit(as.character(row$Adducts), " & "))
      adduct_masses <- c(row$Adduct1Mass, row$Adduct2Mass)
      expanded_rows <- lapply(seq_along(adduct_pairs), function(j) {
        new_row <- row
        new_row$PutativeAdduct <- adduct_pairs[j]
        new_row$MonoisotopicMass <- adduct_masses[j]
        return(new_row)
      })
      return(do.call(rbind, expanded_rows))
    } else if (row$Level == "C") {
      expanded_rows <- lapply(1:nrow(adducts_C_df), function(j) {
        new_row <- row
        new_row$PutativeAdduct <- paste0("[M", adducts_C_df$name[j], "]")
        return(new_row)
      })
      return(do.call(rbind, expanded_rows))
    }
  }))
  return(df_expanded)
}

annotation_df_expanded <- expand_adducts(annotation_df, adducts_C_df)
annotation_df_expanded <- annotation_df_expanded[!duplicated(annotation_df_expanded), ]

#Order columns
column_order <- c("PutativeAdduct", "Level", setdiff(names(annotation_df_expanded), c("PutativeAdduct", "Level")))
annotation_df_expanded <- annotation_df_expanded[, column_order]

#Every mass in C level will be assessed considering all the adducts included in the initial element table
#C level are M+0 excluded in A or B levels

##D level:
#Rules:
#Isotopes with a treshold upper than the selected excluded
#masses in A, B and C levels excluded

#ISOTOPE ANOTATION FOR FILTERING:

isotopes<-rMSI2:::isotopeAnnotation(rMSIprocPeakMatrix, params)

#Isotopes in dataframe format:

result2<-isotopes[names(result) != "isotopicPeaks"] 
result2<-result2[names(result2) != "monoisotopicPeaks"] 

list_df<- vector("list",length(result2)) 
names(list_df) <- names(result2)

for(m in seq_along(result2)) { #recorrer índice de num. de M
  Mx<- result2[[m]] # filtro los valores de M0 cada nivel anotado (M+1, M+1,...)
  param_names<-c("mz_MN","abs_ppm_error", "M0_index","MN_index", "total_score", "morphology_score","intensity_score","mass_error_score","isotopic_ratio","C_atoms")
  mz_values<-names(Mx) #Todos los M0 de cada nivel isotopico anotado
  mat<- matrix(nrow=length(param_names), ncol=length(mz_values))
  rownames(mat)<-param_names
  colnames(mat)<-mz_values
  
  for(mz in seq_along(Mx)) {#indice de cada feature en M
    #isotope_data<-numeric(length(param_names))
    isotope_data<-unname(unlist(Mx[mz]))
    #print(isotope_data)
    mat[,mz]<-isotope_data
  }
  df<-as.data.frame(t(mat))
  list_df[[m]]<-df
}
dfs<-list()
df_names <- c()
for (level in names(list_df)) {
  df_name <- paste0("level_", level)
  #assign(df_name, list_df[[level]], envir = .GlobalEnv) eliminat per modificar l'entorn de treball global
  dfs[[df_name]]<-list_df[[level]]
  df_names<- c(df_names,df_name) #llista de dataframes amb cada nivell isotòpic
  
}

for (df_name in df_names){
  dfs[[df_name]]$isotope_level <- df_name
}

df_unico<- do.call(rbind, dfs) #dataframe que inclou el nivell isotòpic
df_unico$mz_M0<- rMSIprocPeakMatrix$mass[df_unico$M0_index]
df_unico$mz_MN<- rMSIprocPeakMatrix$mass[df_unico$MN_index]
isotopes_df_format <- df_unico[, c(12,1,11,3,4,5,2,6,7,8,9,10)]

isotopes_mass_vector <-rMSIprocPeakMatrix$mass[isotopes$isotopicPeaks]
isotopes_df_format$scoretreshold<- ifelse(isotopes_df_format$mz_MN %in% mass_vector, "YES", "NO")


#FILTERING FOR D GROUP:

all_masses_rounded <-round(rMSIprocPeakMatrix$mass, digits=4)

isotopic_peaks_th<- round(rMSIprocPeakMatrix$mass[isotopes$isotopicPeaks], digits=4)


#Filtering masses annotated in the adduct annotation step and isotopes:
monoisotopic_masses <- round(annotation_df_expanded$MonoisotopicMass, digits=4)
all_masses_rounded <- round(all_masses, digits=4)
masses_without_adducts <- all_masses_rounded[!all_masses_rounded %in% monoisotopic_masses]
masses_without_adducts_isotopes<- masses_without_adducts[!masses_without_adducts %in% isotopic_peaks_th]



adduct_elements <- params$peakAnnotation$adductElementsTable
all_columns <- names(annotation_df_expanded)
new_annotations <- data.frame(matrix(ncol = length(all_columns), nrow = 0))
colnames(new_annotations) <- all_columns

for (mass in missing_masses) {
  for (i in 1:nrow(adduct_elements)) {
    new_row <- data.frame(
      Level = "D",
      MonoisotopicMass = mass,
      PutativeAdduct = paste0("[M", adduct_elements$name[i], "]"),
      stringsAsFactors = FALSE
    )
    new_annotations <- rbind(new_annotations, new_row)
  }
}

for (col in all_columns) {
  if (!col %in% names(new_annotations)) {
    new_annotations[[col]] <- NA
  }
}

new_annotations <- new_annotations[all_columns]

annotation_df_expanded <- rbind(annotation_df_expanded, new_annotations)







#Calcul de probabilitats 


library(dplyr)

#Filtrar A i B
annotation_A_B <- annotation_df_expanded %>%
  filter(Level %in% c("A", "B"))

# Calcular la densitat per aducte i massa neutra en A i B
density_models <- annotation_A_B %>%
  group_by(PutativeAdduct) %>%
  summarise(density = list(density(NeutralMass, na.rm = TRUE)), .groups = 'drop')

# Funció densitat 
get_density_value <- function(mass, density_model) {
  if (is.null(density_model) || length(density_model$x) < 2) {
    return(NA)  # NA si no hi ha suficients dades
  }
  approx(density_model$x, density_model$y, xout = mass, rule = 2)$y
}

# Filtrar C i D
annotation_C_D <- annotation_df_expanded %>%
  filter(Level %in% c("C", "D"))

# Càlcul pb per a C i D
annotation_C_D <- annotation_C_D %>%
  rowwise() %>%
  mutate(Probability = {
    adduct_name <- gsub("\\[M\\+", "", PutativeAdduct)  # Extraer el nombre del aducto
    adduct_name <- gsub("\\]", "", adduct_name)  # Eliminar el corchete de cierre
    
    # Filtrar massa del aducte
    adduct_mass <- params$peakAnnotation$adductElementsTable$mass[
      params$peakAnnotation$adductElementsTable$name == paste0("+", adduct_name)
    ]
    
    if (length(adduct_mass) == 0) {
      warning(paste("No se encontró masa para el aducto:", PutativeAdduct))
      return(NA)
    }
    
    # Càlcul massa d'aducte
    combined_mass <- MonoisotopicMass + adduct_mass
    
    # Valor densitat
    adduct_density <- density_models$density[match(PutativeAdduct, density_models$PutativeAdduct)]
    
    if (length(adduct_density) == 0) {
      warning(paste("No se encontró modelo de densidad para el aducto:", PutativeAdduct))
      return(NA)
    }
    
    # Calcular pb
    prob <- get_density_value(combined_mass, adduct_density[[1]])
    
    prob
  }) %>%
  ungroup()

# Combinar
annotation_df_expanded_v3 <- annotation_df_expanded %>%
  filter(!Level %in% c("C", "D")) %>%
  bind_rows(annotation_C_D)


##Anotar databases



db <- read.csv("/home/nmolto/Desktop/rMSI2_Nadia/PubChemLite_31Oct2020.csv")

match_anotation <- anotation_df_expanded_v3 %>% filter(MonoisotopicMAss %in% db)

###################ADUCTES

MetaboCoreUtils::adducts(polarity = c("positive","negative"))


#####Database

cmpdb <- CompDb(all_minus_CompDb.sqlite)
