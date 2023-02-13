# Duplication analysis with Reconciliation for the publication presenting PhylteR


## Reconciliation with ALE

You can find **Amalgamated likelihood estimation (ALE)** [here](https://github.com/ssolo/ALE).

For each gene found in  the carnivora dataset we reconcile the gene tree with the species tree. Here we force the gene origination to the root of the species tree and set the transfer rate to 0

```sh
ALEml_undated Supermatrix_14463_genes_53_spp_UFBS_TESTNEW_Codon.tre [Gene_tree] \
  sample=100 O_R=10000 separators="" tau=0
```
This output files with the extension `.uml_rec` (contains a summary of event by tree branches) and `.uTs` (contains a summary of transfer events).


ALE renames the species tree internal nodes with numerical ID. They are needed to identify duplication event from `.uml_rec` files. They can be found in the species tree in each uml_rec files (they are the same across all reconciliations using the exact same species tree as input).

```sh
grep -oh "(.*;" `ls *rec | head -n 1` | head -n 1  > aletree
```


## `Duplication_score.R`
The R script is to be run in the directory containing the `.uml_rec` and `aletree` files.
It output `Outlier_Duplication.txt` a data-frame with the duplication score by leaf and by gene.
