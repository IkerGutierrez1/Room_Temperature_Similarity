library(dplyr)
library(zoo)
library(ggplot2)

#Load functions
source("functions.R")

path <- "data/df_RoomT.csv"

df <- load_data(path)

df_no_time <- df[-1,-1]
cols <- names(df_no_time)
cols <- sub("\\.Sensor.*", "", cols)

result_SVD <- svd(df_no_time)


U <- result_SVD$u #vecotores propios a la izquierda (espacio de observaciones 7391x6)
D <- result_SVD$d #Valor absoluto de valores propios (valor propio)
V <- result_SVD$v #vectores propios a la derecha (espacio de vatriables 6x6)

k <- nrow(V)
m <- ncol(V)

W <- numeric(m)
Y <- matrix(0, nrow = k, ncol = m)

for (i in 1:k){
  wkk <- 0
  for (j in 1:m) {
    wkk <- wkk + ((D[j])^2)*(V[i,j])^2
  }
  W[i] <- wkk
}

#Revisar esto
for (i in 1:k){
  for (j in 1:m) {
    ykm <- (((D[j])^2)*((V[i,j])^2))/W[i]
    Y[i,j] <- ykm
  }
}

for (i in 1:k){
  wkk <- 0
  for (j in 1:m) {
    wkk <- wkk + ((D[j])^1)*(V[i,j])^2
  }
  W[i] <- wkk
}

#Revisar esto
for (i in 1:k){
  for (j in 1:m) {
    ykm <- (((D[j])^1)*((V[i,j])^2))/W[i]
    Y[i,j] <- ykm
  }
}

colors <- c("#4C4C8A", "#6A6A9A", "#7F7FB3", "#A3A3C2", "#B2B2D3", "#C2C2E6")
#colors <- c("#C45A0D", "#A02020", "#4A90E2", "#4E9F3D", "#D4AC0E", "#8B4513")
png("outputs/SVD_Spectrum.png", width = 1080, height = 720)

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

#Los pone en graficas separadas estaria bien cambiarlo a todo en 1
pdf("outputs/cpgrams.pdf")

for(col in cols){
  column <- paste0(col,".Sensor__room_temperature")
  cpgram(df_no_time[[column]])
}

dev.off()



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



