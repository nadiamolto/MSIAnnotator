library(dplyr)
library(usethis)
library(roxygen2)
library(readr)
library(data.table)
library(purrr)
library(tidyr)
library(ggplot2)



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

data_path<- "/home/nmolto/Desktop/datasets/toufik_dataset"
rMSIprocPeakMatrix<- rMSI2::LoadPeakMatrix(file.path(data_path,"240430TM_DHB-ANI_LMWC_brain_GFM4_30um.pkmat"))

#AdductAnnotation <- function(rMSIPeacMatrix, AdductList, tolerance) {}

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
Na2_H<- data.frame(name="+2Na-H", mass=44.971164, priority=0)
K2_H<-data.frame(name="+2K-H", mass=76.919044, priority=0)
params$peakAnnotation$adductElementsTable<- rbind(params$peakAnnotation$adductElementsTable,NH4,Na2_H,K2_H)

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
annotation_df <- annotation_df[!duplicated(annotation_df), ] #As we are filtering by adduct1 and adduct2 index, rows are duplicated


#separate rows for each adduct and include the monoisotopic mass

# Function to expand adducts
expand_adducts <- function(df, adducts_C_df) {
  df_expanded <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    if (row$Level %in% c("A", "B")) {
      adduct_pairs <- unlist(strsplit(as.character(row$Adducts), " & "))
      adduct_masses <- c(row$Adduct1Mass, row$Adduct2Mass)
      adduct_indices <- c(row$Adduct1Index, row$Adduct2Index)
      expanded_rows <- lapply(seq_along(adduct_pairs), function(j) {
        new_row <- row
        new_row$PutativeAdduct <- adduct_pairs[j]
        new_row$AdductMass <- adduct_masses[j]
        new_row$AdductIndex <- adduct_indices[j]
        return(new_row)
      })
      return(do.call(rbind, expanded_rows))
    } else if (row$Level == "C") {
      expanded_rows <- lapply(1:nrow(adducts_C_df), function(j) {
        new_row <- row
        new_row$PutativeAdduct <- paste0("[M", adducts_C_df$name[j], "]")
        new_row$AdductMass <- NA
        new_row$AdductIndex <- NA
        return(new_row)
      })
      return(do.call(rbind, expanded_rows))
    }
  }))
  return(df_expanded)
}

# Expanding the dataframe
annotation_df_expanded <- expand_adducts(annotation_df, adducts_C_df)
annotation_df_expanded <- annotation_df_expanded[!duplicated(annotation_df_expanded), ]

# Calculate the neutral mass for level C
annotation_df_expanded <- annotation_df_expanded %>%
  mutate(AdductName = gsub("\\[M(\\+.*)\\]", "\\1", PutativeAdduct)) %>%
  left_join(adducts_C_df, by = c("AdductName" = "name")) %>%
  mutate(NeutralMass = ifelse(Level == "C", MonoisotopicMass - mass, NeutralMass)) %>%
  select(-AdductName, -mass, -priority)

# Remove old columns Adduct1Mass, Adduct2Mass, Adduct1Index, Adduct2Index
annotation_df_expanded <- annotation_df_expanded %>%
  select(-Adduct1Mass, -Adduct2Mass, -Adduct1Index, -Adduct2Index)

# Order columns
column_order <- c("PutativeAdduct", "Level", "AdductMass", "AdductIndex", setdiff(names(annotation_df_expanded), c("PutativeAdduct", "Level", "AdductMass", "AdductIndex")))
annotation_df_expanded <- annotation_df_expanded[, column_order]
annotation_df_expanded <- annotation_df_expanded %>%  rename(AdductPair = Adducts)

#Every mass in C level will be assessed considerAdducts#Every mass in C level will be assessed considering all the adducts included in the initial element table
#C level are M+0 excluded in A or B levels

##D level:
#Rules:
#Isotopes with a treshold upper than the selected excluded
#masses in A, B and C levels excluded

#ISOTOPE ANOTATION FOR FILTERING:

isotopes<-rMSI2:::isotopeAnnotation(rMSIprocPeakMatrix, params)

#Isotopes in dataframe format:

result2<-isotopes[names(isotopes) != "isotopicPeaks"] 
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
isotopes_df_format$scoretreshold<- ifelse(isotopes_df_format$mz_MN %in% isotopes_mass_vector, "YES", "NO")


#FILTERING FOR D GROUP:

all_masses <-rMSIprocPeakMatrix$mass

isotopic_peaks_th<- all_masses[isotopes$isotopicPeaks]

monoisotopic_masses <- all_masses[isotopes$monoisotopicPeaks]
masses_without_adducts <- all_masses[!all_masses %in% monoisotopic_masses] #A, B and C groups
masses_without_adducts_isotopes<- masses_without_adducts[!masses_without_adducts %in% isotopic_peaks_th] #Discarded in adduct annotation
masses_without_adducts_isotopes_b<- masses_without_adducts_isotopes[!masses_without_adducts_isotopes %in% all_masses[aductes$B$Adduct1Index]] #removing group B mz
masses_without_adducts_isotopes_b<- masses_without_adducts_isotopes_b[!masses_without_adducts_isotopes_b%in% all_masses[aductes$B$Adduct2Index]]
#as B group are monoisotopic peaks + some peak in peakmatrix, we remove both

#Adding D group:
# mz in pkm that have not been annotated yet (filtering previously: isotopes and monoisotopic peaks already annotated)

all_columns <- names(annotation_df_expanded)
new_annotations <- data.frame(matrix(ncol = length(all_columns), nrow = 0))
colnames(new_annotations) <- all_columns

for (mass in masses_without_adducts_isotopes_b) {
  for (i in 1:nrow(adducts_C_df)) {
    mono_index <- which(rMSIprocPeakMatrix$mass==mass)
    
    new_row <- data.frame(
      Level = "D",
      MonoisotopicMass = mass,
      MonoisotopicIndex = ifelse(length(mono_index) > 0, mono_index, NA),
      PutativeAdduct = paste0("[M", adducts_C_df$name[i], "]"),
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

annotation_df_expanded_D <- rbind(annotation_df_expanded, new_annotations)

#Calculate the neutral mass for group D:
# Extract adduct name
annotation_df_expanded_D <- annotation_df_expanded_D %>%
  mutate(AdductName = gsub("\\[M(\\+.*)\\]", "\\1", PutativeAdduct))
# Join with adduct table
annotation_df_expanded_D <- annotation_df_expanded_D %>%
  left_join(adducts_C_df, by = c("AdductName" = "name"))
# Calculate neutral mass
annotation_df_expanded_D <- annotation_df_expanded_D %>%
  mutate(NeutralMass = ifelse(Level == "D", MonoisotopicMass - mass, NeutralMass))
# Remove columns
annotation_df_expanded_D <- annotation_df_expanded_D %>%
  select(-AdductName, -mass, -priority)

#Remove rows with negative values in NeutralMass
annotation_df_expanded_D <- annotation_df_expanded_D %>%
  filter(NeutralMass >= 0)



#Probabilities as score for filtering adducts: first calculate the probability, then weight the probabilities and then normalize them


#Filtering A, B:
annotation_A_B <- annotation_df_expanded_D %>%
  filter(Level %in% c("A", "B"))

# Density function
density_models <- annotation_A_B %>%
  group_by(PutativeAdduct) %>% summarise(density = list(density(NeutralMass, na.rm = TRUE)), .groups = 'drop')

get_density_value <- function(mass, density_model) {
  if (is.null(density_model) || length(density_model$x) < 2) {
    return(NA)  # NA si no hi ha suficients dades
  }
  approx(density_model$x, density_model$y, xout = mass, rule = 2)$y
}

# Filtering C, D
annotation_C_D <- annotation_df_expanded_D %>%
  filter(Level %in% c("C", "D"))

#Frequency

adduct_freq <- annotation_df_expanded_D %>%
  group_by(PutativeAdduct) %>%
  summarise(freq = n(), .groups = 'drop')

# Probability for C, D:
annotation_C_D <- annotation_C_D %>%
  rowwise() %>%
  mutate(WeightedProbability = {
    adduct_name <- gsub("\\[M\\+", "", PutativeAdduct)
    adduct_name <- gsub("\\]", "", adduct_name)
    
    # Filtering masses:
    adduct_mass <- params$peakAnnotation$adductElementsTable$mass[
      params$peakAnnotation$adductElementsTable$name == paste0("+", adduct_name)
    ]
    
    if (length(adduct_mass) == 0) {
      warning(paste("No se encontró masa para el aducto:", PutativeAdduct))
      return(NA)
    }
    
    adduct_density <- density_models$density[match(PutativeAdduct, density_models$PutativeAdduct)]
    
    if (length(adduct_density) == 0) {
      warning(paste("No se encontró modelo de densidad para el aducto:", PutativeAdduct))
      return(NA)
    }
    
    # probability
    prob <- get_density_value(NeutralMass, adduct_density[[1]])
    
    #weight
    adduct_weight <- adduct_freq$freq[match(PutativeAdduct, adduct_freq$PutativeAdduct)]
    weighted_prob <- prob * adduct_weight
    weighted_prob
    
  }) %>%
  ungroup()

#Normalization

total_weighted_probability <- sum(annotation_C_D$WeightedProbability, na.rm = TRUE)  # All the probabilities weighted

annotation_C_D <- annotation_C_D %>%
  mutate(NormalizedProbability = ifelse(is.na(WeightedProbability), NA, WeightedProbability / total_weighted_probability))   # Normalization of total addition

# Combinar
annotation_df_expanded_vf <- annotation_df_expanded_D %>%
  filter(!Level %in% c("C", "D")) %>%
  bind_rows(annotation_C_D)   #OBJECTE FINAL

annotation_df_expanded_vf <- annotation_df_expanded_vf %>%
  select(1:4, MonoisotopicMass, everything())


#Modify AdductMass <- MonoisotopicMass and AdductIndex <- MonoisotopicIndex
annotation_df_expanded_vf$MonoisotopicMass[annotation_df_expanded_vf$Level %in% c("A", "B")] <- 
annotation_df_expanded_vf$AdductMass[annotation_df_expanded_vf$Level %in% c("A", "B")]
annotation_df_expanded_vf$MonoisotopicIndex[annotation_df_expanded_vf$Level %in% c("A", "B")] <- 
annotation_df_expanded_vf$AdductIndex[annotation_df_expanded_vf$Level %in% c("A", "B")]

#Reorganization
annotation_df_expanded_vf <- annotation_df_expanded_vf[, !(names(annotation_df_expanded_vf) %in% c("AdductMass", "AdductIndex"))]
annotation_df_expanded_vf <- annotation_df_expanded_vf %>% select(1:3, MonoisotopicIndex, everything())


##density plot 

density_model_repr <- density_models %>%
  mutate(
    x = map(density, ~ .x$x), 
    y = map(density, ~ .x$y)   
  ) %>%
  select(-density) %>%         
  unnest(cols = c(x, y))        


ggplot(density_model_repr, aes(x = x, y = y, color = PutativeAdduct)) +
  geom_line(linewidth = 1) +  # Representar líneas de densidad
  labs(
    title = "Density function for each adduct based on A and B levels",
    x = "NeutralMass",
    y = "Density of probabilities",
    color = "PutativeAdduct"
  ) +
  theme_minimal()



##Anotar databases


MS1_2ID <- read.csv("MS1_2ID.csv")

# limits 5 ppm
annotation_df_expanded_vf_a <- annotation_df_expanded_vf %>%
  mutate(LowerMass = NeutralMass * (1 - 5e-6),
         UpperMass = NeutralMass * (1 + 5e-6))

#data.table format
setDT(annotation_df_expanded_vf_a)
setDT(MS1_2ID)

# columns in db
MS1_2ID[, `:=`(Start = MonoisotopicMass, End = MonoisotopicMass)]

# set keys
setkey(annotation_df_expanded_vf_a, LowerMass, UpperMass)
setkey(MS1_2ID, Start, End)

# merge
match_annotation_1 <- foverlaps(MS1_2ID, annotation_df_expanded_vf_a, 
                                by.x = c("Start", "End"),
                                by.y = c("LowerMass", "UpperMass"),
                                type = "within")

match_annotation_1<- match_annotation_1[,-c("Start","End")]
match_annotation_1 <- match_annotation_1[!is.na(match_annotation_1$PutativeAdduct), ]
match_annotation_1<-unique(match_annotation_1)
match_annotation_1 <- match_annotation_1 %>%  select(-AdductPriority, -LowerMass, -UpperMass)




# SPATIAL ANNOTATION

# 1. Filtering A and B levels since C and D have been assigned to all possible adducts included in the adduct table
filtered_data <- match_annotation_1 %>% 
  filter(Level %in% c("A", "B")) %>%  arrange(NeutralMass)


#2. Setting a tolerance to group neutral masses according to a tolerance

assign_mass_groups <- function(masses, tolerance) {
  group <- 1
  groups <- numeric(length(masses))
  groups[1] <- group
  
  for (i in 2:length(masses)) {
    if (abs(masses[i] - masses[i - 1]) <= tolerance) {
      groups[i] <- group
    } else {
      group <- group + 1
      groups[i] <- group
    }
  }
  
  return(groups)
}

filtered_data <- filtered_data %>%
  arrange(NeutralMass) %>%
  mutate(MassGroup = assign_mass_groups(NeutralMass, 0.001))

#3. Filtering neutral masses
neutral_masses_groups <- filtered_data %>%
  group_by(MassGroup) %>%
  nest()

#4. Correlation
# Intensities normalization
normalizedintensities <- rMSIprocPeakMatrix$intensity / rMSIprocPeakMatrix$normalizations$TIC
correlations_list <- list()

# Function applied to each group
for (group in 1:nrow(neutral_masses_groups)) {
  # Index for each group
  group_data <- neutral_masses_groups$data[[group]]
  # Extract neutral mass
  sorted_group_data <- group_data[order(group_data$Level), ]
  neutral_mass <- sorted_group_data$NeutralMass[1]
  # Group adducts (unique adduct)
  adducts_in_group <- unique(group_data$PutativeAdduct)
  intensities_by_adduct <- list()
  for (adduct in adducts_in_group) {
    # Filter for each UNIQUE adduct
    adduct_data <- group_data %>% filter(PutativeAdduct == adduct) %>% slice_sample(n = 1)
    # Extract intensities
    intensity_values <- normalizedintensities[, adduct_data$MonoisotopicIndex]
    # Intensities by adduct
    intensities_by_adduct[[adduct]] <- intensity_values
  }
  # Correlation calculation
  if (length(intensities_by_adduct) > 1) {
    intensity_matrix <- do.call(cbind, intensities_by_adduct)
    cor_matrix <- cor(intensity_matrix, use = "complete.obs")
    # Name of neutral mass to each correlation matrix
    correlations_list[[as.character(neutral_mass)]] <- cor_matrix
  }
}

#Interpretation
for (neutral_mass in names(correlation_matrices)){}
correlation_matrices <- correlations_list
  cor_matrix <- correlation_matrices$`189.04141`
  pheatmap(cor_matrix, 
           main = paste("Correlations for Neutral Mass:", "189.04141" ),
           cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50))
  
  
  
### Tables with adducts and frequencies for each neutral mass
  

  generate_adducts_table <- function(neutral_masses_groups) {
    adducts_tables_list <- list()
    for (i in 1:nrow(neutral_masses_groups)) {
      annotations_df <- neutral_masses_groups$data[[i]]
      if ("PutativeAdduct" %in% colnames(annotations_df)) {
        neutral_mass <- annotations_df$NeutralMass[1]
        adducts <- annotations_df$PutativeAdduct
        adducts_df <- data.frame(Adduct = adducts)
        adducts_freq <- adducts_df %>%
          group_by(Adduct) %>%
          summarise(Frequency = n()) %>%
          arrange(desc(Frequency))
        adducts_tables_list[[as.character(neutral_mass)]] <- adducts_freq
      } else {
        warning(paste("La columna 'PutativeAdduct' no se encuentra en el tibble para el grupo:", neutral_masses_groups$MassGroup[i]))
      }
    }
    return(adducts_tables_list)
  }
  adducts_tables_list <- generate_adducts_table(neutral_masses_groups)
  for (neutral_mass in names(adducts_tables_list)) {
    cat("\nAdduct table for neutral mass:", neutral_mass, "\n")
    print(adducts_tables_list[[neutral_mass]])
  }
  
  