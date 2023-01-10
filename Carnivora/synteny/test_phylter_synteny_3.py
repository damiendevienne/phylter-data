#
#  DONNES : on regarde une liste de marqueurs orthologues sur deux espèces. Chaque marqueur a donc une coordonnée (chromosome, start, end,strand) dans chaque espèce
#
#  PRINCIPE : on compte le nombre de ruptures de synténie entre les deux espèces dans trois configurations
#  - tous les marqueurs
#  - en enlevant les marqueurs filtrés par phylter
#  - en enlevant le même nombre de marqueurs pris au hasard dans la liste (100 réplicats)
#
#  DEFINITION D'UNE RUPTURE DE SYNTENIE : deux marqueurs consécutifs dans une espèce et pas dans l'autre (la direction compte)
#
#  METHODE : on trie les marqueurs en fonction d'une espèce, et pour chaque couple consécutifs, on regarde s'ils sont consécutifs dans l'autre.
#
#
#  PARAMETRES : deux fichiers pour les génomes des deux espèces (l'orthologie est déduite du nom du gène) + le fichier phylter
#
#
#
#
#
#
#
#

import sys,string,os,functools,random

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

def compute_breaks(markers_keys,species1,species2,output=False):
#    print("list of",len(markers_keys),"markers")
    distrib_longueurs = []
    markers_keys.sort(key = functools.cmp_to_key(lambda x,y: compare(markers[x][species1],markers[y][species1])))
    
    k = 0
    while k < len(markers_keys):
            markers[markers_keys[k]][species1][4] = k+1
            #print(k,markers[markers_keys[k]][species1])
            k = k + 1
            
    markers_keys.sort(key = functools.cmp_to_key(lambda x,y: compare(markers[x][species2],markers[y][species2])))
    
    longueur = 1
    breaks = 0
    double_breaks = 0
    previous_break = True
    k = 0
    while k < len(markers_keys) - 1:
        if (markers[markers_keys[k]][species2][0] != markers[markers_keys[k+1]][species2][0] or
            #markers[markers_keys[k]][species1][0] != markers[markers_keys[k+1]][species1][0] or
            markers[markers_keys[k+1]][species1][4] * markers[markers_keys[k+1]][species1][3] * markers[markers_keys[k+1]][species2][3] -
            markers[markers_keys[k]][species1][4] * markers[markers_keys[k]][species1][3] * markers[markers_keys[k]][species2][3] != 1):
        #if (markers[markers_keys[k]][species2][1] != markers[markers_keys[k+1]][species2][1] or
        #    markers[markers_keys[k]][species1][1] != markers[markers_keys[k+1]][species1][1] or
        #    abs(markers[markers_keys[k+1]][species1][5] - markers[markers_keys[k]][species1][5]) != 1):
            breaks = breaks + 1
            if previous_break:
                double_breaks = double_breaks + 1
            previous_break = True
            distrib_longueurs.append(longueur)
            longueur = 1
            #print("break",markers[markers_keys[k]][species1],markers[markers_keys[k+1]][species1],markers[markers_keys[k]][species2],markers[markers_keys[k+1]][species2])
        else:
            previous_break = False
            longueur = longueur + 1
            #print("non break")
        k = k + 1
                    
    return breaks,double_breaks,distrib_longueurs


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
    if len(parameters) == 5:
        REPLICATES = int(parameters[4])
    else:
        REPLICATES = 100
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
            
    markers_keys.sort(key = functools.cmp_to_key(lambda x,y: compare(markers[x][species1],markers[y][species1])))
    
    k = 0
    while k < len(markers_keys):
            markers[markers_keys[k]][species1][4] = [k+1,(markers_keys[k] in list_of_phylter_markers and
            (species1 in list_of_phylter_markers[markers_keys[k]] or 
             species2 in list_of_phylter_markers[markers_keys[k]]))]
            #print(k,markers[markers_keys[k]][species1])
            k = k + 1
            
    markers_keys.sort(key = functools.cmp_to_key(lambda x,y: compare(markers[x][species2],markers[y][species2])))
    for m in markers_keys:
        file_ouput.write(species1+" "+
                         markers[m][species1][0]+" "+
                         str(markers[m][species1][1])+" "+
                         str(markers[m][species1][2])+" "+
                         str(markers[m][species1][3])+" "+
                         species2+" "+
                         markers[m][species2][0]+" "+
                         str(markers[m][species2][1])+" "+
                         str(markers[m][species2][2])+" "+
                         str(markers[m][species2][3])+" "+
                         str(markers[m][species1][4][1])+"\n")
            
    markers_keys_filtered = markers_keys.copy()
    k = 0
    while k < len(markers_keys_filtered):
        if (markers_keys_filtered[k] in list_of_phylter_markers and
            (species1 in list_of_phylter_markers[markers_keys_filtered[k]] or 
             species2 in list_of_phylter_markers[markers_keys_filtered[k]])):
            del markers_keys_filtered[k]
        else:
            k = k + 1
    
    number_removed = len(markers_keys) - len(markers_keys_filtered)
    
    print("   there are ",len(markers_keys),"markers common to the two species")
    print("   among them,",number_removed,"markers are removed by phylter")
    
    
    observed_breaks,observed_double_breaks,distrib = compute_breaks(markers_keys,species1,species2)
    print("   there are ",observed_breaks,"synteny breaks and",observed_double_breaks,"double breaks")
    file_output = open("without_filter","w")
    for d in distrib:
        file_output.write(str(d)+"\n")

    observed_breaks2,observed_double_breaks2,distrib2 = compute_breaks(markers_keys_filtered,species1,species2)
    print("   there remains", (observed_breaks2,observed_double_breaks2),"synteny breaks and double breaks")
    file_output = open("with_filter","w")
    for d in distrib2:
        file_output.write(str(d)+"\n")
        
    print("computing a p-value with",REPLICATES,"replicates")
    
    number_breaks_smaller_than_observed = 0
    number_boudle_breaks__smaller_than_observed = 0
    for nb_experiences in range(REPLICATES):
        markers_keys_filtered = markers_keys.copy()
        positions_to_filter = []
        while len(positions_to_filter) < number_removed:
            proposition = int(random.random()*len(markers_keys_filtered))
            if not proposition in positions_to_filter:
                positions_to_filter.append(proposition)
        k = 0
        pos = 0
        while k < len(markers_keys_filtered):
            if pos in positions_to_filter:
                del markers_keys_filtered[k]
            else:
                k = k + 1
            pos = pos + 1
    
        breaks,double_breaks,distrib2=compute_breaks(markers_keys_filtered,species1,species2)
        #print(breaks,double_breaks)
        if breaks <= observed_breaks2:
            number_breaks_smaller_than_observed = number_breaks_smaller_than_observed + 1
        if double_breaks <= observed_double_breaks2:
            number_boudle_breaks__smaller_than_observed = number_boudle_breaks__smaller_than_observed + 1
            
    print("   p-values for breaks " +str(number_breaks_smaller_than_observed)+"/"+str(REPLICATES), "p-value for double-breaks "+str(number_boudle_breaks__smaller_than_observed)+"/"+str(REPLICATES))
    
    print(len(markers_keys),"\t",number_removed,"\t",observed_breaks,"\t",observed_double_breaks, "\t",observed_breaks2,"\t",observed_double_breaks2,"\t", number_breaks_smaller_than_observed/REPLICATES, "\t",number_boudle_breaks__smaller_than_observed/REPLICATES)
#for k in markers_keys:
        #print(markers[k][species2],markers[k][species1])
