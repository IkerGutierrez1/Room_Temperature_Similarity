library(dplyr)
library(zoo)
library(ggplot2)
library(tidyr)

#Load functions
source("functions.R")

path <- "data/df_RoomT.csv"

df <- load_data(path)

df_no_time <- df[-1,-1] #Eliminar columna de NAs
cols <- names(df_no_time)
cols <- sub("\\.Sensor.*", "", cols)

#Test
X <- as.matrix(df_no_time)

Xt <- t(X)

X_centered <- scale(X, center = TRUE, scale = FALSE)
X_standar <- scale(X)

media <- attr(X_standar, "scaled:center")
desviacion_estandar <- attr(X_standar, "scaled:scale")

result_SVD <- svd(X_centered)


U <- result_SVD$u #vecotores propios a la izquierda (espacio de observaciones 7391x6)
D <- result_SVD$d #Valor absoluto de valores propios (valor propio) #Singular values
V <- result_SVD$v #vectores propios a la derecha (espacio de vatriables 6x6)



k <- nrow(V)
m <- ncol(V)

n <- nrow(U)
eigen_values <- ((D)^2)/(n-1)

W <- numeric(m)
Y <- matrix(0, nrow = k, ncol = m)

#Se podira hacer en un solo ciclo probablemente
for (i in 1:k){
  wkk <- 0
  for (j in 1:m) {
    wkk <- wkk + ((D[j])^1)*(V[i,j])^2
  }
  W[i] <- wkk
}


for (i in 1:k){
  for (j in 1:m) {
    ykm <- (((D[j])^1)*((V[i,j])^2))/W[i]
    Y[i,j] <- ykm
  }
}

#Mehotd 2
for (i in 1:k){
  wkk <- 0
  for (j in 1:m) {
    wkk <- wkk + (((D[j])^2))/(n-1)*(V[i,j])^2
  }
  W[i] <- wkk
}

for (i in 1:k){
  for (j in 1:m) {
    ykm <- (((((D[j])^2))/(n-1))*((V[i,j])^2))/W[i]
    Y[i,j] <- ykm
  }
}
colors <- c("#4C4C8A", "#6A6A9A", "#7F7FB3", "#A3A3C2", "#B2B2D3", "#C2C2E6")
#colors <- c("#C45A0D", "#A02020", "#4A90E2", "#4E9F3D", "#D4AC0E", "#8B4513")
#png("outputs/SVD_Spectrum.png", width = 1080, height = 720)
png("SVD_Spectrum.png", width = 1080, height = 720)

# Crear un barplot
barplot(Y, beside = TRUE, 
        col = colors,  # Cambia el color según tus preferencias
        xlab = "Principal components",
        ylab = "SVD-Spectrum",
        names.arg = 1:6,  # Etiquetas en el eje x
        ylim = c(0, 1),  # Establece el rango del eje y
        axes = FALSE,
        legend.text = TRUE)  # Agrega leyenda si lo deseas

axis(2, at = seq(0, 1, by = 0.1), las = 1)
grid()
# Agregar leyenda
legend("topright",  # Posición de la leyenda (puedes cambiar a "topleft", "bottomright", etc.)
       legend = cols,  # Nombres correspondientes
       fill = colors,
       #cex = 1,
       #text.width = 0.6,
       cex = 1.2,  # Aumenta el tamaño del texto
       text.width = strwidth(max(cols)) * 1.2,  # Aumenta el ancho de la caja
       box.lwd = 0.5,
       bg = "lightgray")  

dev.off()

#Cumulative PDF
pdf("outputs/cpgrams.pdf")

for(col in cols){
  column <- paste0(col,".Sensor__room_temperature")
  cpgram(df_no_time[[column]],ci.col = "white")
}

dev.off()

#Cumulative 6 pngs
for(col in cols){
  png(filename = paste0("outputs/cumulative/", column, "_cpgram.png"), width = 800, height = 600)
  
  column <- paste0(col,".Sensor__room_temperature")
  cpgram(df_no_time[[column]],ci.col = "white")
  
  dev.off()
}





#Otra forma de calucular
dfxdfT <- t(as.matrix(df_no_time)) %*% (as.matrix(df_no_time))
vp <- eigen(dfxdfT) #La raiz de estos correponde con D



# Inicializa un data.frame vacío para almacenar los datos
combined_data <- data.frame()

for (col in cols) {
  column <- paste0(col, ".Sensor__room_temperature")
  
  # Agrega los datos al data.frame combinado
  temp_data <- df_no_time[[column]]
  combined_data <- rbind(combined_data, data.frame(Temperatura = temp_data, Column = column))
}

png("outputs/density.png", width = 1080, height = 720)

# Crea el gráfico con todas las curvas
ggplot(combined_data, aes(x = Temperatura, color = Column)) +
  geom_density() + 
  labs(title = "T Density", x = "Temperature", y = "Density") +
  theme_minimal()
# Cierra el archivo PDF
dev.off()




#Reconstruct original
#Los 6 componentes principales (Los dos metodos tienen que ser lo mismo)
k <- 1
S <- diag(D)
Sk <- S[1:k,1:k]
Uk <- U[,1:k]
principal <- Uk%*%Sk #Si k = 1 no funciona multiplicar un escalar  (componente principal)
principal <- Uk*Sk

T_build <- as.matrix(principal)+sum(media)/length(media) #Forma de conseguir una temperatura

Vk <- t(V)[1:k,]
Xk <- as.matrix(principal)%*%Vk 


Xk <- Xk+media


#Move to dataframe whith timestamps
Xk_df <- data.frame(Xk)
Xk_mean <- data.frame(T_build = rowMeans(Xk_df[,1:6]))
names(Xk_df) <- cols
Xk_df$timestamp <- df$timestamp[-1]
Xk_mean$timestamp <- df$timestamp[-1]
T_build_df <- data.frame(T_build)
T_build_df$timestamp <- df$timestamp[-1]

#Scatter de Xk_df y X para comparar las reconstrucciones
#Save dfs
write.csv2(T_build_df, "k1_T_build.csv", row.names = FALSE)


# Añadir una columna para identificar el origen de los datos
df_graph <- df[-1,]

colnames(df_graph) <- paste("df_",colnames(df_graph), sep = "")
df_combined <- bind_cols(df_graph, Xk_df)

df_filtered <- df_combined %>%
  filter(timestamp >= "2023-08-07" & timestamp <= "2023-08-08")

df_filtered$timestamp <- as.POSIXct(df_filtered$timestamp, format = "%Y-%m-%d %H:%M:%S")






#All plots
base_save_dir = "outputs/plots/SWD_standarized/"
is_standarized = TRUE
is_centered = FALSE

start_date = "2023-08-07"
end_date = "2023-08-08"
end_date = "2023-08-14"

date_breaks = "8 hours"
date_breaks = "4 hours"

plot_type = "daily"
plot_type = "weakly"


for (k in 1:6){
  save_dir = paste0(base_save_dir,"k=",k,"/")
  if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
  }
  S <- diag(D)
  Sk <- S[1:k,1:k]
  Uk <- U[,1:k]
  if (k != 1){
    principal <- Uk%*%Sk} #Si k = 1 no funciona multiplicar un escalar  (componente principal)}
  else{
    principal <- Uk*Sk}
  
  Vk <- t(V)[1:k,]
  Xk <- as.matrix(principal)%*%Vk 
  
  #Linea solo si no se ha estandarizado al principio
  if (is_standarized){
    Xk <- Xk*desviacion_estandar+media}
  
  else if(is_centered){
    Xk <- Xk + media
  }
  
  Xk_df <- data.frame(Xk)
  names(Xk_df) <- cols
  Xk_df$timestamp <- df$timestamp[-1]
  
  #No hace falta que se ejecute en cada ciclo
  df_graph <- df[-1,]
  
  colnames(df_graph) <- paste("df_",colnames(df_graph), sep = "")
  df_combined <- bind_cols(df_graph, Xk_df)
  
  df_filtered <- df_combined %>%
    filter(timestamp >= start_date & timestamp <= end_date)
  
  df_filtered$timestamp <- as.POSIXct(df_filtered$timestamp, format = "%Y-%m-%d %H:%M:%S")
  
  #Scatter
  file_path = paste0(save_dir,"scatter.PNG")
  p <- ggplot(df_combined, aes(x = df_RoomB.Sensor__room_temperature, y = RoomB.Sensor__room_temperature)) +
    geom_point() +  # Esta capa agrega los puntos al gráfico
    labs(
      x = "Original",
      y = "SWD"
    ) +
    theme_bw()
  
  ggsave(file_path, plot = p, width = 8, height = 6, dpi = 300)
  
  for (col in names(Xk_df)){
    if (!grepl("timestamp", col)) {
      file_path = paste0(save_dir,plot_type,"/")
      if (!dir.exists(file_path)){
        dir.create(file_path, recursive = TRUE)
      }
      file_path = paste0(file_path,col,".PNG")
      print(file_path)
      original_col <- paste0("df_",col)
  
      df_long <- df_filtered %>%
        pivot_longer(cols = c(!!col, !!original_col),
                     names_to = "Source",
                     values_to = "Temperature")
      
      colors <- c(col = "red", 
                  original_col = "blue")
      
      # Cambiar los nombres de las columnas (que están como strings) por sus valores dinámicos
      names(colors) <- c(col, original_col)
      
      
      # Crear la gráfica
      g <- ggplot(df_long, aes(x = timestamp, y = Temperature, color = Source)) + 
        geom_point(size = 2, alpha = 0.7) +
        geom_line() +
        
        # Título y etiquetas
        labs(
          title = col,
          x = "Timestamp",
          y = "Temperature (°C)"
        ) +
        
        # Personalización de la apariencia
        theme_bw() +  # Fondo blanco
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
          axis.title = element_text(size = 12, face = "bold"),  
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
        ) +
        
        # Ajustar las fechas del eje X
        scale_x_datetime(date_breaks = date_breaks, date_labels = "%Y-%m-%d %H:%M") +
        
        # Limitar el eje Y entre 18 y 30
        ylim(18, 30) +
        
        scale_color_manual(values = colors, 
                           labels = c("Original", "SWD"))
      
      
      # Guardar la gráfica como archivo PNG
      ggsave(file_path, plot = g, width = 8, height = 6, dpi = 300)
      
  }}}





#Write errors
fileConn <- file("results.txt", open = "w")

writeLines("Ciclo\t\tRoom\tMAE\tRMSE", con = fileConn)

is_standarized = FALSE

for (k in 1:6){
  S <- diag(D)
  Sk <- S[1:k,1:k]
  Uk <- U[,1:k]
  if (k != 1){
    principal <- Uk%*%Sk} #Si k = 1 no funciona multiplicar un escalar  (componente principal)}
  else{
    principal <- Uk*Sk}
    
  Vk <- t(V)[1:k,]
  Xk <- as.matrix(principal)%*%Vk 
  
  #Linea solo si no se ha estandarizado al principio
  if (is_standarized){
  Xk <- Xk*desviacion_estandar+media}

  err_M <- X-Xk
  
  err_Abs <- abs(X-Xk)
  
  
  mae_per_variable <- colMeans(err_Abs)
  
  err_M_squared <- err_M^2
  
  mse <- colMeans(err_M_squared)
  rmse <- sqrt(mse)
  
  for (j in 1:length(mae_per_variable)){
    writeLines(paste(k, names(mae_per_variable)[j],mae_per_variable[j],rmse[j],sep="\t"), con = fileConn)
  }
  writeLines("", con = fileConn)
}
close(fileConn)


colnames <- names(df)
file_name = "outputs/"

for (col in colnames){
  
  x_var <- df[[col]][-1]
  y_var <- Xk_df[[col]]
  
  temp_df <- data.frame(x = x_var, y = y_var)
  
  p <- ggplot(temp_df, aes(x = x, y = y)) +
    geom_point(color = "red", size = 2) +
    labs(title = "Real vsd PCA",
         x = "Real measurment",
         y = "Reconstructed by PCA") +
    theme_bw() +  # Aplicar tema con fondo blanco
    theme(
      plot.background = element_rect(fill = "white"), # Fondo blanco para todo el gráfico
      panel.grid.major = element_line(color = "grey90"), # Líneas de cuadrícula más suaves
      panel.grid.minor = element_blank()               # Sin líneas menores
    )
  
  file_name <- paste0("outputs/scatter/scatter_", col, ".png")
  ggsave(file_name, plot = p, width = 8, height = 8, units = "in", dpi = 100)
}


df_long <- df %>%
  pivot_longer(cols = -timestamp, names_to = "variable", values_to = "value")

df_filtered <- df_long %>%
  filter(timestamp >= as.POSIXct("2023-11-13 00:00:00") & timestamp <= as.POSIXct("2023-11-13 23:59:59"))

# Crear el gráfico
p <- ggplot(df_filtered, aes(x = timestamp, y = value, color = variable)) +
  geom_point() +  # Usar geom_point() si prefieres puntos en lugar de líneas
  labs(title = "Valores de las columnas frente a timestamp",
       x = "Timestamp",
       y = "Valor") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Eliminar el título de la leyenda

# Mostrar el gráfico
print(p)
