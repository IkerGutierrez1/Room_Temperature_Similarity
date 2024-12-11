#Load functions
source("functions.R")

#Read data
path <- "data/df_RoomT.csv"
df <- load_data(path)

#Delte col with NA
df_no_time <- df[-1,-1] 
cols <- names(df_no_time)

#SVD
X <- as.matrix(df_mean_no_t)
X_centered <- scale(X, center = TRUE, scale = FALSE)

means <- attr(X_centered, "scaled:center")

result_SVD <- svd(X_centered)
U <- result_SVD$u #Left singular vec
D <- result_SVD$d #Singular values
V <- result_SVD$v #Right singular vec

#SVD Spectrum calculation
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

#Reconstruct singals and global signal
k <- 1
S <- diag(D)
Sk <- S[1:k,1:k]
Uk <- U[,1:k]
if (k != 1) {
  principal <- Uk %*% Sk  # Matricial prodcut k no 1
} else {
  principal <- Uk * Sk   # PProduct k 1
}

#Add mean to obtain temperature measurments for globals signal
T_build <- as.matrix(principal)+sum(mean)/length(mean) 

#Reconstruction of all signals
Vk <- t(V)[1:k,]
Xk <- as.matrix(principal)%*%Vk 

#Add mean to recover temperatures
Xk <- Xk+means


#Move to dataframe whith timestamps for plot
Xk_df <- data.frame(Xk)
Xk_mean <- data.frame(T_build = rowMeans(Xk_df[,1:6]))
names(Xk_df) <- cols

Xk_df$timestamp <- df$timestamp[-1]
Xk_mean$timestamp <- df$timestamp[-1]

T_build_df <- data.frame(T_build)
T_build_df$timestamp <- df$timestamp[-1]


#Save dfs
write.csv2(T_build_df, "k1_T_build.csv", row.names = FALSE)

#--------------------------
#Plot spectrum
colors <- c("#4C4C8A", "#6A6A9A", "#7F7FB3", "#A3A3C2", "#B2B2D3", "#C2C2E6")
png("SVD_Spectrum.png", width = 1080, height = 720)
barplot(Y, beside = TRUE, 
        col = colors, 
        xlab = "Principal components",
        ylab = "SVD-Spectrum",
        names.arg = 1:6,  
        ylim = c(0, 1), 
        axes = FALSE,
        legend.text = TRUE)  
axis(2, at = seq(0, 1, by = 0.1), las = 1)
grid()
legend("topright", 
       legend = cols,  
       fill = colors,
       cex = 1.2,  
       text.width = strwidth(max(cols)) * 1.2, 
       box.lwd = 0.5,
       bg = "lightgray")  

dev.off()

#Cumulative PDF plots
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
#-----------------------------------------

#Plots for all the reconstructions + scatter plots 
#---------------------------------------------------------
#All plots of reconstruction for all rooms for all principal components and scatter plost of reconstructed
#real data
base_save_dir = "outputs/plots/SWD_test/"
is_standarized = FALSE
is_centered = FALSE

start_date = "2023-08-07"
start_date = "2023-03-07"
end_date = "2023-08-08"
end_date = "2023-08-14"

date_breaks = "8 hours"
date_breaks = "4 hours"

plot_type = "daily"
plot_type = "weakly"


for (k in 1:1){
  save_dir = paste0(base_save_dir,"k=",k,"/")
  if (!dir.exists(save_dir)){
    dir.create(save_dir, recursive = TRUE)
  }
  S <- diag(D)
  Sk <- S[1:k,1:k]
  Uk <- U[,1:k]
  if (k != 1){
    principal <- Uk%*%Sk} 
  else{
    principal <- Uk*Sk}
  
  Vk <- t(V)[1:k,]
  Xk <- as.matrix(principal)%*%Vk 
  
  if (is_standarized){
    Xk <- Xk*desviacion_estandar+media}
  
  else if(is_centered){
    Xk <- Xk + media
  }
  
  
  Xk_df <- data.frame(Xk)
  Bk_df <- data.frame(principal)
  names(Xk_df) <- cols
  Xk_df$timestamp <- df$timestamp[-1]
  
  df_graph <- df[-1,]
  
  colnames(df_graph) <- paste("df_",colnames(df_graph), sep = "")
  df_combined <- bind_cols(df_graph, Xk_df)
  
  df_filtered <- df_combined %>%
    filter(timestamp >= start_date & timestamp <= end_date)
  
  df_filtered$timestamp <- as.POSIXct(df_filtered$timestamp, format = "%Y-%m-%d %H:%M:%S")
  
  #Scatter
  file_path = paste0(save_dir,"scatter.PNG")
  p <- ggplot(df_combined, aes(x = df_RoomB.Sensor__room_temperature, y = RoomB.Sensor__room_temperature)) +
    geom_point() +  
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
      
      names(colors) <- c(col, original_col)
      
      
      g <- ggplot(df_long, aes(x = timestamp, y = Temperature, color = Source)) + 
        geom_point(size = 2, alpha = 0.7) +
        geom_line() +
        
        labs(
          title = col,
          x = "Timestamp",
          y = "Temperature (°C)"
        ) +
        
        theme_bw() +  
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
          axis.title = element_text(size = 12, face = "bold"),  
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
        ) +
        
        scale_x_datetime(date_breaks = date_breaks, date_labels = "%Y-%m-%d %H:%M") +
        
        ylim(18, 30) +
        
        scale_color_manual(values = colors, 
                           labels = c("Original", "SWD"))
      
      
      ggsave(file_path, plot = g, width = 8, height = 6, dpi = 300)
      
    }}}





#Calculte and save residuals
fileConn <- file("results.txt", open = "w")

writeLines("Ciclo\t\tRoom\tMAE\tRMSE", con = fileConn)

is_standarized = FALSE

for (k in 1:6){
  S <- diag(D)
  Sk <- S[1:k,1:k]
  Uk <- U[,1:k]
  if (k != 1){
    principal <- Uk%*%Sk} 
  else{
    principal <- Uk*Sk}
  
  Vk <- t(V)[1:k,]
  Xk <- as.matrix(principal)%*%Vk 
  
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
    theme_bw() +  
    theme(
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey90"), 
      panel.grid.minor = element_blank()               
    )
  
  file_name <- paste0("outputs/scatter/scatter_", col, ".png")
  ggsave(file_name, plot = p, width = 8, height = 8, units = "in", dpi = 100)
}


df_long <- df %>%
  pivot_longer(cols = -timestamp, names_to = "variable", values_to = "value")

df_filtered <- df_long %>%
  filter(timestamp >= as.POSIXct("2023-11-13 00:00:00") & timestamp <= as.POSIXct("2023-11-13 23:59:59"))


p <- ggplot(df_filtered, aes(x = timestamp, y = value, color = variable)) +
  geom_point() +  
  labs(title = "Valores de las columnas frente a timestamp",
       x = "Timestamp",
       y = "Valor") +
  theme_minimal() +
  theme(legend.title = element_blank()) 

print(p)
#------------------------------------------

#Centralized with previous day mean
#Add colum with previous day mean, SVD should be perform with the fist part of the code introducing df_mean_no_t to SVD
df_mean_prev <- df %>%
  mutate(date = as.Date(timestamp)) %>%
  group_by(date) %>%
  summarise(
    across(RoomA.Sensor__room_temperature:RoomF.Sensor__room_temperature, 
           mean, .names = "mean_{.col}"), 
    .groups = "drop"
  ) %>%
  mutate(
    across(starts_with("mean_"), lag, .names = "previous_day_{.col}")
  ) %>%
  right_join(
    df %>% mutate(date = as.Date(timestamp)), 
    by = "date"
  ) %>%
  select(-starts_with("mean_"))

df_mean <- df_mean_prev %>%
  mutate(across(
    starts_with("Room"), 
    ~ . - get(paste0("previous_day_mean_", cur_column())),
    .names = "{.col}"
  ))

df_mean <- na.omit(df_mean)
df_daily_mean <- df_mean %>%
  select(starts_with("previous_day_mean"))
df_mean <- df_mean %>%
  select(-contains("previous_day_mean"), -date)

df_mean_no_t <- df_mean[,-1]


#Reconstruct data
Xk <- Xk + df_daily_mean 
Xk_df <- as.data.frame(Xk)
Xk_df$timestamp <- df_mean$timestamp


#Pre process
cols <- names(df_no_time)
new_col_names <- paste0(cols, "_reconstructed")
new_col_names <- c(new_col_names,"timestamp")

colnames(Xk_df) <- new_col_names


df_mean_prev <- df_mean_prev %>%
  select(-contains("previous_day_mean"), -date)
df_mean_prev <- df_mean_prev %>%
  mutate(timestamp = if_else(
    grepl("^\\d{4}-\\d{2}-\\d{2}$", timestamp),  
    paste0(timestamp, " 00:00:00"),              
    timestamp                                    
  ))
Xk_df <- Xk_df %>%
  mutate(timestamp = if_else(
    grepl("^\\d{4}-\\d{2}-\\d{2}$", timestamp), 
    paste0(timestamp, " 00:00:00"),              
    timestamp                                    
  ))
df_test <- left_join(df_mean_prev, Xk_df, by = "timestamp") 
df_test <- na.omit(df_test)

start_date <- as.POSIXct("2023-08-07 00:00:00")  
end_date <- as.POSIXct("2023-08-13 23:59:59")   

# Filtrar las filas dentro del rango de fechas
df_test <- df_test %>%
  filter(timestamp >= start_date & timestamp <= end_date)


p <- ggplot(data = df_test, aes(x = timestamp)) +
  # Puntos originales
  geom_point(aes(y = RoomB.Sensor__room_temperature, color = "Original")) +
  # Línea que conecta los puntos originales
  geom_line(aes(y = RoomB.Sensor__room_temperature, color = "Original", group = 1)) +  # group = 1 asegura que todas las observaciones se conecten
  # Puntos reconstruidos
  geom_point(aes(y = RoomB.Sensor__room_temperature_reconstructed, color = "Reconstructed")) +
  # Línea que conecta los puntos reconstruidos
  geom_line(aes(y = RoomB.Sensor__room_temperature_reconstructed, color = "Reconstructed", group = 1)) +  # group = 1 asegura que todas las observaciones se conecten
  # Etiquetas y tema
  labs(
    title = "RoomB original vs reconstructed",
    x = "Timestamp",
    y = "Valores",
    color = "Leyenda"
  ) +
  theme_minimal()

# Mostrar el gráfico
p

Xk_principal <- principal + rowMeans(df_daily_mean) #Single data reconstruction, to obtain temperature data the mean of the previous day mean

df_principal_plot <- data.frame(
  timestamp = df_mean$timestamp,
  reconstructed = Xk_principal
)

start_date <- as.POSIXct("2023-08-07 00:00:00")  
end_date <- as.POSIXct("2023-08-13 23:59:59")   

# Filtrar las filas dentro del rango de fechas
df_principal_plot <- df_principal_plot %>%
  filter(timestamp >= start_date & timestamp <= end_date)

p_princial <- ggplot(data = df_principal_plot, aes(x = timestamp, y = reconstructed)) +
  geom_point(color = "red") +               # Puntos individuales
  geom_line(color = "red", group = 1) +   # Línea que conecta los puntos
  labs(
    title = "Resultados frente a Timestamps",
    x = "Timestamp",
    y = "Resultado"
  ) +
  theme_minimal()

p_princial

