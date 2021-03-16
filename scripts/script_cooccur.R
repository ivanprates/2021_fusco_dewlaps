##########################
### R script by Ivan Prates (ivanprates.org).
### Smithsonian National Museum of Natural History, Washington, DC, USA.
### November 2019.

### The goals of this R script are:
### 1. To test for reduced co-occurrence between the three major dewlap phenotypes of Anolis fuscoauratus and other sympatric anole species. 
### 2. To plot the results of co-occurrence analyses.

## PART 1: Getting ready.

## Install packages:
#install.packages("cooccur")

## Load packages:
library(cooccur)

## Selecting and setting working directory:
path = "~/Dropbox/Science/MYPAPERS_ongoing/2020_fuscoauratus/2020_gh_fuscoauratus/sympatric_anoles/"
#path = "C:/Users/RhinellaX/Dropbox/Science/MYPAPERS_ongoing/2018_fuscoauratus/sympatric_species"
setwd(path)

## Importing list data:
input_data = read.csv(file = paste0(path, "sympatric_Anolis_matrix_fusco_as_0-1_2020-02-29.csv"), header = TRUE, row.names = "locality")
dewlap_presabs_matrix = t(input_data[, 6:25])
#View(dewlap_presabs_matrix)

## Selecting test to run (each corresponding to a combination of species):
#test_number = "test_1" ## Using all species; this analyses is basically the same below because coccur excludes species with too few occurrences.
test_number = "test_2" ## Using only species with more than five occurrence records.
#test_number = "test_3" ## Separating species into two groups based on relative reflectance: more and less reflective (regardless of clade or microhabitat). 
#test_number = "test_4" ## Separating species into two groups based on clade: Dactyloa and Norops (regardless of dewlap color)
#test_number = "test_5" ## Separating species into three groups based on clade and relative dewlap reflectance: Dactyloa, Draconura less reflective, and Draconura more reflective.

## PART 2: Selecting the comparisons to be made.

## Test 1:
## Dewlaps of A. fuscoauratus versus co-ocurrence with all other anole species:
test_1 = dewlap_presabs_matrix[c("Anolis_fuscoauratus_GRAY", "Anolis_fuscoauratus_PINK", "Anolis_fuscoauratus_YELLOW",
                                          "Anolis_auratus",
                                          "Anolis_chrysolepis",
                                          "Anolis_dissimilis",
                                          "Anolis_nasofrontalis",
                                          "Anolis_ortonii",
                                          "Anolis_planiceps",
                                          #"Anolis_pseudotigrinus",
                                          "Anolis_punctatus",
                                          "Anolis_scypheus",
                                          "Anolis_tandai",
                                          "Anolis_trachyderma",
                                          "Anolis_transversalis"
                                          ), ]

## Test 2:
## Dewlaps of A. fuscoauratus versus co-ocurrence with other anole species that have at least five records:
test_2 = dewlap_presabs_matrix[c("Anolis_fuscoauratus_GRAY", "Anolis_fuscoauratus_PINK", "Anolis_fuscoauratus_YELLOW",
                                          "Anolis_ortonii",
                                          "Anolis_punctatus",
                                          "Anolis_tandai",
                                          "Anolis_trachyderma",
                                          "Anolis_transversalis"
                                          ), ]


## Test 3:
## Dewlaps of A. fuscoauratus versus more (white, yellow) or less (black, blue, red) reflective colors in other anoles, regardless of habitat:
test_3 = dewlap_presabs_matrix[c("Anolis_fuscoauratus_GRAY", "Anolis_fuscoauratus_PINK", "Anolis_fuscoauratus_YELLOW",
                                          "Anolis_more_reflective",
                                          "Anolis_less_reflective"
                                          ), ]

## Test 4:
## Dewlaps of A. fuscoauratus versus more or less reflective color in other anoles, segregated by clade:
test_4 = dewlap_presabs_matrix[c("Anolis_fuscoauratus_GRAY", "Anolis_fuscoauratus_PINK", "Anolis_fuscoauratus_YELLOW",
                                          "Draconura", 
                                          "Dactyloa" 
                                          ), ]

## Test 5:
## Dewlaps of A. fuscoauratus versus more or less reflective color in other anoles, segregated by clade, with Draconura segregated by color reflectiveness:
test_5 = dewlap_presabs_matrix[c("Anolis_fuscoauratus_GRAY", "Anolis_fuscoauratus_PINK", "Anolis_fuscoauratus_YELLOW",
                                          "Draconura_less_reflective", 
                                          "Draconura_more_reflective", 
                                          "Dactyloa" 
                                          ), ]

## PART 3: Testing for reduced co-occurrences.

## Selecting data based on desired test:
if (test_number == "test_1") {
  test_matrix = test_1
} else if (test_number == "test_2") {
  test_matrix = test_2
} else if (test_number == "test_3") {
  test_matrix = test_3
} else if (test_number == "test_4") {
  test_matrix = test_4
} else if (test_number == "test_5") {
  test_matrix = test_5
}

## Running the cooccur algorithm:
cooccur_dewlaps = cooccur(test_matrix, type = "spp_site", spp_names = TRUE,
                          prob = "hyper") ## "hyper" is the approach of Griffith et al. 2016, which uses a more succint and faster combinatorial approach. 

## Creating new pdf file to save pairwise co-occurrence plot:
pdf(file = paste0(path, "/cooccur_", test_number, "_pairwise.pdf"), height = 10, width = 12)

## Visualizing pairwise co-occurrence plot:
plot(cooccur_dewlaps, plotrand = TRUE) # add "plotrand = TRUE" to include randomly associated species species

## Saving pairwise co-occurrence plot:
dev.off()

## Printing the total positive, negative, and random species pairs classified:
## Will also give the number of species pairs that were not classifiable due to low statistical power.
summary(cooccur_dewlaps)
  
## Checking significance levels for negative and positive co-occurrence patterns:
cooccur_output = prob.table(cooccur_dewlaps)
cooccur_output

## Saving these outputs as a text file:
capture.output(cooccur_output, file = paste0(test_number, ".txt"))

## Inspect detailed results for specific dewlap phenotypes:
pair(mod = cooccur_dewlaps, spp = "Anolis_fuscoauratus_GRAY", all = FALSE)
pair(mod = cooccur_dewlaps, spp = "Anolis_fuscoauratus_PINK", all = FALSE)
pair(mod = cooccur_dewlaps, spp = "Anolis_fuscoauratus_YELLOW", all = FALSE)

## Producing a table of the percentage of each species total pairings that were classified as positive, negative, and random:
pair.attributes(cooccur_dewlaps)

## Visualing these percentages in a bar plot:
pair.profile(mod = cooccur_dewlaps)

## Creating new pdf file to save plot of observed versus expected co-occurrences:
pdf(file = paste0(path, "/cooccur_", test_number, "_observed_vs_expected.pdf"), height = 10, width = 12)

## Plotting the observed versus expected co-occurrences:
obs.v.exp(cooccur_dewlaps)

## Saving observed versus expected co-occurrences plot:
dev.off()

effect.sizes(cooccur_dewlaps)

## Done!

