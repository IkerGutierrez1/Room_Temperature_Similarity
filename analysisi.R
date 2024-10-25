library(dplyr)
library(zoo)

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
       cex = 0.7,
       text.width = 0.6,
       box.lwd = 0.5,
       bg = "lightgray")  


acf(df_no_time$RoomA.Sensor__room_temperature)
cpgram(df_no_time$RoomA.Sensor__room_temperature)
cpgram(df_no_time$RoomB.Sensor__room_temperature)

resultados <- lapply(df_no_time, function(col) {
  cpgram(col, plot = FALSE)
})





#Otra forma de calucular
dfxdfT <- t(as.matrix(df_no_time)) %*% (as.matrix(df_no_time))
vp <- eigen(dfxdfT) #La raiz de estos correponde con D
