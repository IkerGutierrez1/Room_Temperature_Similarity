load_data <- function(path){
  df <- read.csv2(path)
  
  df <- df %>%
    filter(!is.na(timestamp))
  
  df <- interpolate_NA(df)
  
  return (df)
}


interpolate_NA <- function(df){
  df <- df %>%
    mutate(across(-timestamp, ~ as.numeric(as.character(.))))
  
  #Función para hacer una interpolación de los datos que tenga NA, asumen que timestamp
  #esta correctamente preprocesado
  col_NA <- colSums(is.na(df))
  print(col_NA)
  #Hay 13 valores de yTi que son NA y 195 de Te
  for(col in names(col_NA)){
    if(col_NA[col] > 0){
      print(paste("Columna:", col, "- Interpolando", col_NA[col], "valores NA"))
      df <- df %>%
        mutate(!!col := na.approx(.[[col]], na.rm = FALSE))
    }
  }
  
  # Eliminar filas donde todas las columnas (excepto 'timestamp' y 't') son NA, 
  #la interpolacion no funciona si no hay valoresm validos antes de NA, por lo que
  #los NA de la primera hora del dataset se eliminan aqui
  df <- df %>%
    filter(!rowSums(is.na(select(., -c(timestamp)))) == (ncol(df) - 2))
  return (df)
}
