# Reconciliation analysis to compute duplication score associated to each sequence in the Carnivora dataset

A duplication score was computed for each sequence of the Carnivora dataset. These scores are available in the file `Outlier_Duplication.tar.gz`.

This score was obtained by *reconciling* each individual gene tree with the Carnivora species tree (https://github.com/damiendevienne/phylter-data/blob/main/Carnivora/Supermatrix_14463_genes_53_spp_UFBS_TESTNEW_Codon.tre) and computing the mean number of event infered along the path from the root of the species tree to each species for each gene. 

Reconciliations were computed with ALE (Amalgamated likelihood estimation, available at https://github.com/ssolo/ALE) and the duplication score was computed with the homemade R script `Duplication_score.R`. Here are the details: 

- Reconciliation with ALE
Here is the way ALE was ran:

```sh
ALEml_undated Supermatrix_14463_genes_53_spp_UFBS_TESTNEW_Codon.tre [Gene_tree] \
  sample=100 O_R=10000 separators="" tau=0
```
This output files with extensions `.uml_rec` (contains a summary of event by tree branches) and `.uTs` (contains a summary of transfer events).

ALE renames the species tree internal nodes with numerical ID. that are needed to identify duplication event from `.uml_rec` files. They can be found in the species tree in each uml_rec files (they are the same across all reconciliations using the exact same species tree as input).

```sh
grep -oh "(.*;" `ls *rec | head -n 1` | head -n 1  > aletree
```

- Duplication score with `Duplication_score.R`
The R script is to be run in the directory containing the `.uml_rec` and `aletree` files.
It output `Outlier_Duplication.txt` a data-frame with the duplication score by leaf and by gene.
