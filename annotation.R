

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
result<-rMSI2:::isotopeAnnotation(rMSIprocPeakMatrix, params)
aductes <- rMSI2::peakAnnotation(rMSIprocPeakMatrix, params=params)


## MSIAnnotator

#settings
params<-rMSI2::ProcessingParameters()
params$peakAnnotation$ppmMassTolerance<-15
NH4<- data.frame(name="+NH4", mass=18.033823, priority=0)
params$peakAnnotation$adductElementsTable<- rbind(params$peakAnnotation$adductElementsTable,NH4)
aductes <- rMSI2::peakAnnotation(rMSIprocPeakMatrix, params=params)

#format

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


adducts_C_df <- data.frame(
  name = c("+K", "+Na", "+H", "+NH4"),
  mass = c(38.96, 22.99, 1.01, 18.03),
  priority = c(0, 0, 0, 0)
)

expand_adducts <- function(df, adducts_C_df) {
  df_expanded <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    if (row$Level %in% c("A", "B")) {
      adduct_pairs <- unlist(strsplit(as.character(row$Adducts), " & "))
      expanded_rows <- lapply(adduct_pairs, function(adduct) {
        new_row <- row
        new_row$PutativeAdduct <- adduct
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

column_order <- c("PutativeAdduct", "Level", setdiff(names(annotation_df_expanded), c("PutativeAdduct", "Level")))
annotation_df_expanded <- annotation_df_expanded[, column_order]

#Fins aquí m'assigna al nivell C tots els aductes de la taula

## Masses faltants de la peakmatrix

new_masses <- rMSIprocPeakMatrix$mass

missing_masses <- new_masses[!new_masses %in% annotation_df_expanded$MonoisotopicMass]

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

# Filtrar niveles A y B
annotation_A_B <- annotation_df_expanded %>%
  filter(Level %in% c("A", "B"))

# Calcular la densidad para cada aducto usando NeutralMass
density_models <- annotation_A_B %>%
  group_by(PutativeAdduct) %>%
  summarise(density = list(density(NeutralMass, na.rm = TRUE)), .groups = 'drop')

# Función para obtener la densidad
get_density_value <- function(mass, density_model) {
  if (is.null(density_model) || length(density_model$x) < 2) {
    return(NA)  # Devuelve NA si no hay datos suficientes para la densidad
  }
  approx(density_model$x, density_model$y, xout = mass, rule = 2)$y
}

# Filtrar niveles C y D
annotation_C_D <- annotation_df_expanded %>%
  filter(Level %in% c("C", "D"))

# Calcular las probabilidades basadas en la densidad para cada aducto en los niveles C y D
annotation_C_D <- annotation_C_D %>%
  rowwise() %>%
  mutate(Probability = {
    adduct_name <- gsub("\\[M\\+", "", PutativeAdduct)  # Extraer el nombre del aducto
    adduct_name <- gsub("\\]", "", adduct_name)  # Eliminar el corchete de cierre
    
    # Encontrar la masa del aducto correspondiente
    adduct_mass <- params$peakAnnotation$adductElementsTable$mass[
      params$peakAnnotation$adductElementsTable$name == paste0("+", adduct_name)
    ]
    
    # Si no se encuentra la masa del aducto, emitir una advertencia y continuar
    if (length(adduct_mass) == 0) {
      warning(paste("No se encontró masa para el aducto:", PutativeAdduct))
      return(NA)
    }
    
    # Calcular la masa combinada
    combined_mass <- MonoisotopicMass + adduct_mass
    
    # Obtener el modelo de densidad para el aducto
    adduct_density <- density_models$density[match(PutativeAdduct, density_models$PutativeAdduct)]
    
    # Si no se encuentra el modelo de densidad, emitir una advertencia y continuar
    if (length(adduct_density) == 0) {
      warning(paste("No se encontró modelo de densidad para el aducto:", PutativeAdduct))
      return(NA)
    }
    
    # Calcular la probabilidad de densidad para la masa combinada
    prob <- get_density_value(combined_mass, adduct_density[[1]])
    
    # Devolver la probabilidad o NA si no se pudo calcular
    prob
  }) %>%
  ungroup()

# Combinar los resultados con el dataframe original sin modificar A y B
annotation_df_expanded_v3 <- annotation_df_expanded %>%
  filter(!Level %in% c("C", "D")) %>%
  bind_rows(annotation_C_D)




###################ADUCTES

MetaboCoreUtils::adducts(polarity = c("positive","negative"))


#####Database

cmpdb <- CompDb(all_minus_CompDb.sqlite)

