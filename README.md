This folder contains all the files used to run various test regarding the paper.
Per now, the directory is not sorted.

Code is written in R. Runs on R 4.3.3
Necessary packages are (list to be completed):
- R.matlab
- geoR
- MASS
- ggplot2
- fields
- akima

prior trend: AlvheimPrior.csv
data: Alvheim_horizontal_classify_with_griddata.mat

1) Loading of the real data
To load the prior trend and true data: simply run the all_code/Alvheim_data.R script

2) Generating values for AVO intercept and Gradient using the full rock physics forward model:
and generate the AVO intercept and gradient which the MARS is compared to, the following files can be ran:
in the folder forward_model, there is an additional read me with instructions regarding input values  needed to run the files etc.

Additional files are added for various plotting and some other computations done for the paper
