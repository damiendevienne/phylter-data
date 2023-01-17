# Synteny analysis for the publication presenting PhylteR 

## `test_phylter_synteny_5_2.py`
The python script `test_phylter_synteny_5_2.py` takes as input 
- the file named `genomes` that contains, for all the 14 genomes analyzed, the name of the genes, their associated scaffold and their coordinates and strand n these scaffolds. 
- the name of two species (to get the list of 14 species, the script above can be run without any argument)
- a list of outliers such as those given in https://github.com/damiendevienne/phylter-data/tree/main/Carnivora/data/phylter-results (for phylter) or https://github.com/damiendevienne/phylter-data/tree/main/Carnivora/data/treeshrink-results (for treeshrink).

It ouputs the number of genes common to the two species compared, the number of syntenic outliers, the number of outlier sequences (from the list of outliers given) that belong to the species of interest, and the number of syntenic outliers among these outliers. Finally, it also computes the p-value associated the probability of observing this number of syntenic outliers in the list of outliers gibven the total number of syntenic outliers and the total number of genes.

## `code-synteny.R`
The R code `code-synteny.R` allows calling the previous script for all possible pairs of species, formatting the results, saving it as csv files (present in thjis folder also) and draw the figures presented in the manuscript.
