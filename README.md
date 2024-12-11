# Room Temperature Signal Similarity

This repository contains a R project used for the analysis of the similarity between different room tempratures data, differente alternatives to create a global signal for the entire building is analysis. It includes the R scripts and the plots created to visualize the data. 

## Data sources

The data was collected from six office in Aalborg University in Denmark from February 2023 to December 2023. It has a frecuency of five minutes and contains measurments of occupany, ventilation, heating, lighting and general enviromental factors. The dataset contains data from six office rooms. 
The data was preprocessed to change the frecuency to 1 hour.

## Analysis aproach

The analsysis consited in perfomring a Pincipal Component Analysis (PCA) using Singular Value Decomposition (SVD), a spectrum of the decomposition is calculted ploting the energy associated with each component.
A cumulative periodogram is also calulted for each signal and they are compared betwwen them. 

For the reconstruction, two aproaches are followed, since the signals are similar enough the mean of all the singla is calucalted and that is the new global singal, the other aproaches is a data dimensionality reduction 
perfomr with SVD.

## Repository structure

- data contains the .csv file with the input data

- outputs contains plots of the analysis

- analysis.R Is the script to perform the analysis. It first reads data and performs the SVD, calculates the coefficients for the spcetrum W and Y and reconstruct the data in Xk_df and T_build_df. From line 77 the scptrum and cumulative plots are obtained. Line 125 is the code corresponding
to the plot of the comparisson from the reconstruction with differnte k for all the rooms. From line 333 onward a different centralization of the data is tried, using the mean from the previous day.

- functions.R Only contains 2 function used fot loading the data.
