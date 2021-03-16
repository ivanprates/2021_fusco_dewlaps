##########################
### R script by Ivan Prates (ivanprates.org).
### Smithsonian National Museum of Natural History, Washington, DC, USA.
### November 2019.

### The goals of this script are:
### To run genetic clustering analyses using SNMF;
### To make bar plots based on the ancestry coefficients.

## Part 1: Getting ready.

## Installing LEA:
#source("https://bioconductor.org/BiocManager.R")
#biocLite("LEA")
## Obs: try "http://" if "https://" URLs are not supported.
#devtools::install_github("bcm-uga/LEA") ## Alternatively.

## Installing other packages we'll need:
#install.packages("cowplot")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("magrittr")
#install.packages("scico") # For color palettes.
#install.packages("tidyr")
#install.packages("viridis") # For color palettes.

## Loading packages:
library(cowplot)
library(dplyr)
library(ggplot2)
library(LEA)
library(magrittr)
library(tidyr)
library(scico)
library(viridis)

## Setting the number of repetitions:
## For test purposes, used 3. Used 100 in final analyses (reported in the manuscript).
r = 3 
#r = 100

# Setting the number of computer cores we want to use:
c = 4 

## Picking general path, depending on which computer I'm working at:
path = "~/Dropbox/Science/MYPAPERS_ongoing/2020_fuscoauratus/2020_gh_fuscoauratus/"
#path = "~/Dropbox/2018_fuscoauratus/"
#path = "C:/Users/RhinellaX/Dropbox/science/MYPAPERS_ongoing/2018_fuscoauratus/"

## Creating a directory to place LEA results:
dir.create(paste0(path, "LEA"))

## Creating a directory to place resulting SNMF bar plots:
dir.create(paste0(path, "LEA/SNMF_plots"))

## Creating a directory to place resulting qmatrices:
dir.create(paste0(path, "LEA/SNMF_qmatrices"))

## Creating a directory to place resulting workspaces:
dir.create(paste0(path, "LEA/SNMF_workspaces"))

## Creating a directory for the outputs of SNMF analyses:
dir.create(paste0(path, "LEA/SNMF_analyses"))

## Setting this new directory as the working directory:
setwd(paste0(path, "LEA/SNMF_analyses"))

## Setting a single input dataset:
input_file = "fusco_R12_c95_n162_m120"

## Reading genetic data:
## Missing sites should be coded as 9! Don't forget to replace -1 to 9 if using VCF Tools.
gendata = read.table(paste0(path, "VCFtools_SNMF/", input_file, "/", input_file, "_1SNP-locus.012"), sep = "\t", row.names = 1) # read genotype file from vcftools
dim(gendata)

## Writing down data in LFMM format (LEA saves it all outside R, instead of as objects):
write.lfmm(gendata, paste0(input_file, "_genotypes.lfmm"))

## Part 2: Running genetic clustering analyses using sNMF:

## Setting a, sNMF's regularization parameter:
## For test purposes, setting a single value.
a = 400

## For the actual analyses, setting alpha as a loop so we can iterate over different values: 
#for (a in c(1,25,50,100,200,400,800,1600,3200)){

## Running sNMF:
project.snmf = snmf(paste0(input_file, "_genotypes.lfmm"), CPU = c, ploidy = 2,
                    entropy = TRUE, alpha = a, repetitions = r, project = "new",
                    #K = 4:6 ## Values used for testing.
                    K = 1:12 ## Values used in manuscript.
                    )

## Showing project structure:
show(project.snmf)

## Showing summary of project results:
snmf_summary = summary(project.snmf)
snmf_summary

## Selecting criterion to determine the best K:
#crossEntropy = snmf_summary$crossEntropy[1,] # Using mininum cross entropy across runs.
crossEntropy = snmf_summary$crossEntropy[2,] # Using mean cross entropy across runs.

## Selecting best K based on minimum (or mean) cross entropy among runs:
bestK = names(which.min(crossEntropy))
bestK %<>% gsub(pattern = "K = ", replacement = "", .) # removes gaps from sequences. "." denotes object when using "%<>%"
bestK = as.integer(bestK)
#bestK = 6 ## Resulting best K (was K = 6).
bestK ## View best K

## Setting a couple variables to help us name outputs down the road:
mincrossentropy = min(crossEntropy)
minentropy = round(mincrossentropy, digits = 2) # Rounding digits.

## Getting the cross entropy of all runs for K = best K:
ce = cross.entropy(project.snmf, K = bestK)
ce

## Selecting the run with the lowest cross-entropy for K = best K:
bestrun = which.min(ce)
bestrun

## Getting the entropy value that corresponds to the best run for the best K:
e = round(ce[bestrun], digits = 2) # Rounding digits.
e

## Plotting cross-entropy criterion of all runs of the project:
## Note that this will plot minimum, not mean, entropy scores, so bestK in plots may differ from when using mean values.
png(filename = paste0(path, "LEA/SNMF_plots/", input_file, "_a", a, "_e", minentropy, "_K", bestK, ".png"))
plot(project.snmf, cex = 1.2, col = "blue", pch = 19) 
dev.off() # Saving plot as .png image.

## Part 3: Preparing data to plot ancestry coefficients and add sample info:

## Getting the Q matrix for the best run and K:
qmatrix = Q(project.snmf, K = bestK, run = bestrun)
qmatrix = as.data.frame(qmatrix)
colnames(qmatrix) = c(paste0(rep("Cluster_", bestK), 1:bestK)) # Change column names for cluster name. Will repeat "cluster" K times, then paste0 with second element from 1 to K times.

## Adding individual ID names for plotting:
individuals = read.table(file = paste0(path, "VCFtools_SNMF/", input_file, "/", input_file, "_1SNP-locus.012.indv")) # Read list of sample IDs from vcftools.
qmatrix$ID = individuals$V1 # Change row names to sample IDs.

## "Melt" dataframe using the gather function so we can assign samples to clusters based on max qscores:
qmatrix_melt = gather(qmatrix, key = cluster_assign, value = coeff, 1:bestK)

## Assign specimens to cluster based on the highest qscore (coeff) values:
cluster_assign = qmatrix_melt %>% group_by(ID) %>% top_n(n = 1, wt = coeff)

## Ordering rows to allow combination with a data frame with sample information:
cluster_assign = arrange(cluster_assign, ID)

## Do samples in qmatrix and sample_info file match?
table(qmatrix$ID == cluster_assign$ID) # Must all be "TRUE".

## Now add cluster assigment to qmatrix:
qmatrix$cluster_assign = cluster_assign$cluster_assign

## We can now round the qscore values to 2 digits for clarity:
qmatrix %<>% mutate_at(paste0("Cluster_", 1:bestK), round, digits = 2)

## Now, let's add locality information (including lat-longs) to the qmatrix:
sample_info = read.csv(file = paste0(path, "sample_information/2019-05_ddRAD_samples_n164_fuscoauratus.csv"), header = TRUE)

# Removing two samples not used in SNMF due to high levels of missing data.
sample_info = subset(sample_info, ID != "H0805")
sample_info = subset(sample_info, ID != "H4109")

## Do samples in qmatrix and sample_info file match?
table(qmatrix$ID %in% sample_info$ID) # Must all be "TRUE".

## Adding locality info for the qmatrix. This will be useful to plot maps later.
qmatrix$state = sample_info$state
qmatrix$country = sample_info$country
qmatrix$locality = sample_info$locality
qmatrix$latitude = sample_info$latitude
qmatrix$longitude = sample_info$longitude
qmatrix$dewlap_group = sample_info$dewlap_group

## Arranging order of samples in qmatrix to make beautiful plots:
## WAY 1: Using the sample order that resulted from a phylogenetic analysis using RaxML.
sample_order = read.csv(file = paste0(path, "LEA/SNMF_plot_order/SNMF_plot_order_n162_2020-02-26.csv"), header = TRUE)
colnames(sample_order) = c("ID", "plot_order")
sample_order = arrange(sample_order, ID)

## Do samples in qmatrix and sample_order file match?
table(qmatrix$ID == sample_order$ID) # Must all be "TRUE".

## Adding a row with the sample order into the qmatrix based on the order given by the sample_info file:
qmatrix$plot_order = sample_order$plot_order

## Saving this individual-based qmatrix:
write.csv(qmatrix, file = paste0(path, "LEA/SNMF_qmatrices/qmatrix_", input_file, "_a", a, "_K", bestK, ".csv"), quote = FALSE, row.names = FALSE)

## We also want to calculate the average qscores by locality to plot qscore pie charts on a nice map.
## Calculating mean qscores per locality:
mean_qscores_by_locality = qmatrix %>% group_by(country, state, locality, latitude, longitude) %>% summarise_at(c(paste0("Cluster_", 1:bestK)), mean) 

## Also in this case we want to round the qscores to 2 digits:
mean_qscores_by_locality %<>% mutate_at(paste0("Cluster_", 1:bestK), round, digits = 2)

## Writing this locality-level qmatrix:
write.csv(mean_qscores_by_locality, file = paste0(path, "LEA/SNMF_qmatrices/qmatrix_mean_qscores_by_locality_", input_file, "_a", a, "_K", bestK, ".csv"), row.names = FALSE)

## Now let's do the final adjustment to prepare the data for ggplot.

## Replacing unknown dewlap colors by NA:
qmatrix$dewlap_group %<>% gsub(pattern = "U", replacement = NA)

## Ordering rows by pre-defined plotting order (e.g., sample order given by fineRADstructure or RaxML):
qmatrix = arrange(qmatrix, plot_order)

## "Melt" dataframe using the gather function, as required by ggplot:
qmatrix_remelt = gather(qmatrix, key = sNMF_cluster, value = coeff, 1:bestK)

## To ensure the order of samples in plots, transform columns into factors, defining levels:
qmatrix_remelt$ID = factor(qmatrix_remelt$ID, levels = unique(qmatrix$ID))
qmatrix_remelt$sNMF_cluster = factor(qmatrix_remelt$sNMF_cluster, levels = unique(qmatrix$cluster_assign))

## PART 4: Plotting ancestry coefficients as barplots with ggplot:

## Let's first Select a color palette to make bar plots on the individual qscores:
palette_name = "tokyo_mod"
myPalette = c("#3C1835", "#652570", "#495D8F", "#BEDA9F", "#7EAC92", "#FFFFD8") ## Setting manually.
myPalette = rev(myPalette)

## Plotting with ggplot:
SNMF_plot = 
  
  ggplot(data = qmatrix_remelt, aes(x = ID)) +
  
  ## Adding bars that represent ancestry coefficients:
  geom_bar(aes(y = coeff, fill = sNMF_cluster), 
           color = "gray95", ## Color of bar lines.
           size = 0, ## Thickness of bar lines.
           stat = "identity", position = "fill", show.legend = FALSE) +
  
  ## Adding points that correspond to dewlap color below bars:
  geom_point(aes(x = ID, y = -0.02, color = dewlap_group), size = 5) +
  
  ## Adding points that correpond to localities below bars:
  ## Comment out "scale_color_manual" below for this to work.
  #geom_point(aes(x = ID, y = -0.02, color = locality), size = 2.5) +
  #guides(color = FALSE) +
  
  ## Adding dark empty circles that will encircle the points that show dewlap color:
  #geom_point(aes(x = ID, y = -0.02), shape = 21, color = "gray30", size = 2) +
  
  ## Filling bars by cluster assigment:
  scale_fill_manual("sNMF Cluster", values = myPalette[c(1:bestK)]) + 
  
  ## Filling circles by dewlap color:
  scale_color_manual(values = c("#c5beab", "#ec95b4", "#f4ee32"), guide = "none") +
  
  ## Adjusting labels:
  labs(y = "Ancestry coefficient", x = "") +
  
  ## Adjusting limits of the y axis:
  scale_y_continuous(limits = c(-0.03,1), expand = c(0,0)) +
  
  ## Changing overall font size:
  theme_minimal(base_size = 16) +
  
  ## Changing theme parameters:
  theme(
    
    ## Adding individual IDs below bars:  
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7, margin = margin(-1.5, 0, 0, 0)),
    #axis.text.x = element_blank(), ## Removing IDs from plot.
    
    ## Making plot wider:
    aspect.ratio = 1/5,
    
    ## Adjusting other parameters:
    panel.grid.major.x = element_blank(),
    panel.grid = element_blank(), 
    axis.text.y = element_text(angle = 90, margin = margin(0, -2, 0, 0))
  )

## Checking plot:
SNMF_plot

## Saving plot:
save_plot(filename = paste0(path, "LEA/SNMF_plots/", input_file, "_a", a, "_e", minentropy, "_K", 
                            bestK, "_", palette_name, "_dewlap_group.pdf"), 
          plot = SNMF_plot, base_width = 30, base_height = 15)
      
## Closing loops:
#} # First loop (a values)
#} # Second loop (input files)

## PART 6: Wrapping up:

## Saving our progress as an R workspace image:
save.image(file = paste0(path, "LEA/SNMF_workspaces/", input_file, "_a", a, ".RData"))

## Clearing working space
#rm(list = ls())

## Done!
