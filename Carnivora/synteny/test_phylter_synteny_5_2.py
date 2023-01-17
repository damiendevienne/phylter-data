#
#  DONNEES : on regarde une liste de marqueurs orthologues sur deux espèces. Chaque marqueur a donc une coordonnée (chromosome, start, end,strand) dans chaque espèce
#
#  DONNEES : on a une liste de marqueurs "à enlever", résultat de Phylter
#
#  EN ENTREE : les noms de deux espèces
#
#  PRINCIPE : on identifie les gènes qui ont une rupture de synténie des deux côtés (les "syntenic outliers")
#  on compte le nombre de syntenic outliers qui sont identifiés par Phylter
#  on calcule une p-value sous l'hypothèse que les syntenic outliers et phylter sont indépendants
#
#  DEFINITION D'UNE RUPTURE DE SYNTENIE : deux marqueurs consécutifs dans une espèce et pas dans l'autre (la direction compte)
#
#  METHODE : on trie les marqueurs en fonction d'une espèce, et pour chaque couple consécutifs, on regarde s'ils sont consécutifs dans l'autre.
#
#
#  PARAMETRES : deux fichiers pour les génomes des deux espèces (l'orthologie est déduite du nom du gène) [+ le fichier phylter + les noms des especes]
#
#
#
#
#
#
#
#

import sys,string,os,functools,random
from math import comb

parameters = sys.argv[1:]

# deux comportements : sans parametre, affiche la liste des espèces et le nombre de scaffolds
# avec les noms de deux espèces (genre_spe), affiche le nombre de ruptures de synténie et une p-value
if len(parameters) != 1 and len(parameters) != 4 and len(parameters) != 5:
    print("usage: python3 test_phylter_synteny.py genomes_file [result_phylter spe1 spe2 [number of replicates=100]]")
    print("      one parameter for a summary list of available species with chromosome or scaffold numbers")
    print("      three parameters for a comparison of species couple, an additional one for the number of replicates")
    exit()

# retourne -1, 0, 1
def cmp(a, b):
    return (a > b) - (a < b) 

# compare deux coordonnées au format [id,chromosome,start,stop]
def compare(x,y):
    if x[0] != y[0]:
        return cmp(x[0],y[0])
    else:
        return cmp(min(int(x[1]),int(x[2])),min(int(y[1]),int(y[2])))

def compute_syntenic_outliers(markers_keys,species1,species2):
    
    result = []

    markers_keys.sort(key = functools.cmp_to_key(lambda x,y: compare(markers[x][species1],markers[y][species1])))
    
    k = 0
    while k < len(markers_keys):
            markers[markers_keys[k]][species1][4] = k+1
            k = k + 1
            
    markers_keys.sort(key = functools.cmp_to_key(lambda x,y: compare(markers[x][species2],markers[y][species2])))

    previous_break = True
    k = 0
    while k < len(markers_keys):
        if ((k == len(markers_keys) - 1) or
            (markers[markers_keys[k]][species2][0] != markers[markers_keys[k+1]][species2][0]) or
            (markers[markers_keys[k]][species1][0] != markers[markers_keys[k+1]][species1][0]) or
            (markers[markers_keys[k+1]][species1][4] * markers[markers_keys[k+1]][species1][3] * markers[markers_keys[k+1]][species2][3] -
             markers[markers_keys[k]][species1][4] * markers[markers_keys[k]][species1][3] * markers[markers_keys[k]][species2][3] != 1)):
            if previous_break:
                result.append(markers_keys[k])
            previous_break = True
        else:
            previous_break = False
        k = k + 1
                    
    return result


print("reading markers")
# format markers : [id] (pour la famille) [species] (pour l'espèce), puis [chromosome,start,stop,strand]
markers = {}
list_species = []

file_list = open(parameters[0],"r").readlines()[1:]

for line in file_list:
    
    words = line.strip().split()
    family_name = words[0].split("_")[2]
    species = "_".join(words[0].split("_")[:2])
    scaffold = words[2].split(":")[0]
    start = int(words[2].split(":")[1][1:].split("-")[0])
    stop = int(words[2].split(":")[1][1:].split("-")[1])
    if words[2].split(":")[1][0] == "+":
        strand = 1
    else:
        strand = -1

    if not family_name in markers:
        markers[family_name] = {}
    markers[family_name][species] = [scaffold,start,stop,strand,0]
        
    if not species in list_species:
        list_species.append(species)


if len(parameters) == 1:
    print("list of available species with chromosome or scaffold numbers")
    for species in list_species:
        list_chromosomes = []
        for m in markers:
            if species in markers[m]:
                chromosome = markers[m][species][0]
                if not chromosome in list_chromosomes+["UNKNOWN","null"]:
                    list_chromosomes.append(chromosome)
        print("   ",species,len(list_chromosomes))
else:
    species1 = parameters[2]
    species2 = parameters[3]
    file_ouput = open(species1+"-"+species2+".markers","w")
    
    if not species1 in list_species or not species2 in list_species:
        print("unknown species: check the available species names by running the program with no parameter")
        print("usage: python3 test_phylter_synteny.py [result_phylter spe1 spe2]")
        exit()
    print("counting the number of synteny breaks between genomes of",species1,"and",species2)
    
    
    print("reading phylter output")
    file_phylter = open(parameters[1],"r").readlines()
    list_of_phylter_markers = {}
    line = 0
    while line < len(file_phylter):
        if file_phylter[line][0] != "#":
            marker_name = file_phylter[line].strip().split()[0].split("_")[0]
            species = file_phylter[line].strip().split()[1]
            if marker_name in list_of_phylter_markers:
                list_of_phylter_markers[marker_name].append(species)
            else:
                list_of_phylter_markers[marker_name] = [species]
        line = line + 1
            
    markers_keys = list(markers.keys())
    k = 0
    while k < len(markers_keys):
        if not (species1 in markers[markers_keys[k]] and species2 in markers[markers_keys[k]]):
            del markers_keys[k]
        else:
            k = k + 1
            
    markers_keys_filtered = []
    k = 0
    while k < len(markers_keys):
        #print(k,markers_keys[k],markers_keys[k] in list_of_phylter_markers)
        if (markers_keys[k] in list_of_phylter_markers and
            (species1 in list_of_phylter_markers[markers_keys[k]] or 
             species2 in list_of_phylter_markers[markers_keys[k]])):
            markers_keys_filtered.append(markers_keys[k])
        k = k + 1
    
    print("   there are ",len(markers_keys),"markers common to the two species")
    print("   among them,",len(markers_keys_filtered),"markers are removed by phylter")
    
    
    syntenic_outliers = compute_syntenic_outliers(markers_keys,species1,species2)
    print("   there are ",len(syntenic_outliers),"syntenic outliers")
    
    intersection = 0
    for k in syntenic_outliers:
        if k in markers_keys_filtered:
            intersection = intersection + 1
    
    n = len(markers_keys)
    k = len(syntenic_outliers)
    l = len(markers_keys_filtered)
    o = intersection
    numerateur=0
    for i in range(o, min(l,k)+1): 
        #print(i)
        numerateur=numerateur+comb(k,i)*comb(n-k, l-i)
    pval = numerateur/comb(n,l)
    print("p-value under independence hypothesis",pval)

    print(n,k,l,o,pval)