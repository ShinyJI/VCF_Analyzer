import sys, re
import matplotlib.pyplot as plt

plt.rcParams['toolbar'] = 'None'

filePath = r"C:\Users\Jack\Desktop\Cours\HMIN113M_Systeme\Projet_variants\human_CEU.vcf"

DATE = ""

# Les cas possibles pour chaque colonne obligatoire : dictionnaire de liste de dictionnaires
# TODO rajouter les cas de base ? (Champs de base pour ALT et les autres -> A,T,C,G ?)
COLONNES = {'CHROM' : [], 'POS' : [], 'ID' : [], 'REF' : [], 'ALT' : [],
            'QUAL' : [], 'FILTER' : [], 'INFO' : [], 'FORMAT' : []}

# Les nom des colonnes facultatives stockés dans un ensemble (pas de doublons)
ADDITIONAL_COLUMNS = set()

# Les informations sur les variants stockées dans une liste de dictionnaires
# De la forme : VARIANTS = [{CHROM : [v1.1, v1.2], POS : [v2.1], ...}, {CHROM: [v2.1]}, ... ]
VARIANTS = list()

# Dictionnaire synchronisant les colonnes avec leur ordre de lecture (Equivalent d'Hashmap)
SYNCHRO_COLUMNS = {0 : 'CHROM', 1 : 'POS', 2 : 'ID', 3 : 'REF', 4 : 'ALT',
                   5 : 'QUAL', 6 : 'FILTER', 7 : 'INFO'}

# TODO
# Transforme la date en une date lisible
def setDate(chaine):
    global DATE
    
    date = chaine


# Check si le fichier est bien un .vcf
def checkName(filePath):
    return filePath.endswith(".vcf")


# Ouvre et renvoie le fichier 
def getFile(filePath):
    try :
        f = open(filePath)
        return f
    except :
        print("Erreur lors de l'ouverture du fichier")
        sys.exit(1)

        
# Ajoute au dictionnaire des cas un cas donné
# Un cas est initialement de la forme : <ID=X,Number=Y,Type=Z,Description=W, ... >
# Il est transformé en la forme : {ID : X, Number : Y, Type : Z, Description : W}
# Champs par défaut si manquants Type -> FLAG, Number -> 0, Description -> ''
def addCase(caseName, case):
    global COLONNES
    
    ID = re.search('ID=[^,>]*', case).group(0)[3:]
    
    try:
        Description = re.search('Description="[^,"]*', case).group(0)[13:]
    except:
        Description = ''

    try :
        Number = re.search('Number=[^,>]*', case).group(0)[7:]
    except :
        Number = 0
    
    try:
        Type = re.search('Type=[^,>]*', case).group(0)[5:]
    except :
        Type = 'FLAG'

    dic = {'ID' : ID, 'Number' : Number, 'Type' : Type, 'Description' : Description}
    
    if dic not in COLONNES[caseName]:   # Vérification doublons :
        COLONNES[caseName].append(dic)  # Dans certains fichiers, les en-têtes sont répétées


# Rempli les listes avec les données d'une ligne
def fillInfos(line):
    global COLONNES
    
    caseName = re.search('^[^=]*', line).group(0)

    if caseName in COLONNES.keys(): # Si c'est une information sur les champs des colonnes
        addCase(caseName, line)
        
    elif caseName == 'fileDate':    # Si c'est la date
        setDate(line[9:])

    #TODO : Rajouter des cas : version, reference, source, contig ?
        
    else :  # Champs non utilisés
        pass


# Rempli l'ensemble des colonnes supplémentaires et synchronise leur ordre
def fillAdditionnalColumns(line):
    global ADDITIONAL_COLUMNS, SYNCHRO_COLUMNS
    
    try :
        column = re.search('FORMAT.*', line).group(0)   # Nom des colonnes apres INFO
        column_list = column.split('\t')    # Transformation en liste
        for column_name in column_list:
            ADDITIONAL_COLUMNS.add(column_name) # Ajout de la colonne à l'ensemble des colonnes
            if column_name not in SYNCHRO_COLUMNS.values(): # Synchronisation
                SYNCHRO_COLUMNS[len(SYNCHRO_COLUMNS)] = column_name 
            
    except :    # Si il n'y a pas de colonnes supplémentaires
        pass


# Rempli la liste de variants avec les informations d'une ligne (1 variant)
def fillVariants(line):
    global SYNCHRO_COLUMNS, VARIANTS
    
    dic = {}
    
    field_list = line.split('\t')
    
    column_number = 0
    for column in field_list :
        elements_list = column.split(';')
        dic[SYNCHRO_COLUMNS[column_number]] = elements_list
        column_number += 1

    VARIANTS.append(dic)


# Rempli les dictionnaires et listes avec les données en entrée
def fillAllDatas(f):
    
    for line in f : # Pour chaque ligne de f
        if line.startswith('##'):   # Lignes d'informations
            fillInfos(line[2:])
        
        elif line.startswith('#'): # Ligne représentant les colonnes
            fillAdditionnalColumns(line)
        
        else:   # Les lignes représentant un variant
            fillVariants(line)
            pass

#########################################################################################################################
#########################################################################################################################
                                                # Fonctions d'analyse

# Recupération des qualités des variants
def QualityByChrom():
    global VARIANTS

    dic = {}    # Dictionnaire {CHROM1 : [QUAL1, QUAL2, ... ], CHROM2 : [QUAL3, ... ], ... }
    maxQual = 0
    
    for variant in VARIANTS:
        CHROM = variant['CHROM'][0] # Chromosome du variant
        
        try:
            QUAL = int(variant['QUAL'][0])   # Qualité du variant
        except:
            QUAL = 0
        
        if QUAL > maxQual : # Mise à jour de la qualité maximale
            maxQual = QUAL
            
        try :
            dic[CHROM].append(QUAL)
        except :
            dic[CHROM] = [QUAL]

    if maxQual > 0 :    # Si la qualité est renseignée
        
        analyseQuality(maxQual, dic)  # Analyse

        x = dic.keys()
        height = [sum(qual)/len(qual) for qual in dic.values()]
        ylabel = 'Qualité'
        xlabel = 'Chromosome'
        title = 'Moyenne de la qualité des variants pour chaque chromosome'
        color = (0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0)
        edgecolor = 'blue'
        
        plotBar(x, height, xlabel, ylabel, title, color, edgecolor)



# Analyse la qualité des variants
def analyseQuality(maxQual, qualityDic):
    
    dic = {}
    for CHROM, QUAL in qualityDic.items():
        nbrOccurence = QUAL.count(maxQual)
        try :
            dic[CHROM] += nbrOccurence/len(QUAL)*100
        except :
            dic[CHROM] = nbrOccurence/len(QUAL)*100

    
    x = dic.keys()
    height = dic.values()
    ylabel = 'Nombre de variant (en %)'
    xlabel = 'Chromosome'
    title = 'Nombre de variants ayant la meilleure qualité (' + str(maxQual)+  ') pour chaque chromosome, en %'
    color = 'orange'
    edgecolor = 'red'
    
    plotBar(x, height, xlabel, ylabel, title, color, edgecolor)


# Occurence des base nucléiques pour tous les variants
def nbrREFTotal():
    global VARIANTS
    
    dic = {} # Dictionnaire des occurences des gènes
    
    for variant in VARIANTS :
        for base in variant['REF']:
            try:
                dic[base] += 1
            except :
                dic[base] = 1

    if sum(dic.values()) > 0 :  # Si il y a des valeurs à afficher
            
        labels = dic.keys()
        sizes = [c for c in dic.values()]
        title = "Nombre d'occurences de chaques base nucléique pour tous les variants"
        
        plotPieChart(sizes, labels, title)


# Occurence des bases nucléiques des variants pour un chromosome donné
def nbrREFByChrom(chrom):
    global VARIANTS
    
    dic = {} # Dictionnaire des occurences des gènes
    
    for variant in VARIANTS :
        if variant['CHROM'][0] == chrom:
            for base in variant['REF']:
                try:
                    dic[base] += 1
                except :
                    dic[base] = 1

    if sum(dic.values()) > 0 :  # Si il y a des valeurs à afficher
            
        labels = dic.keys()
        sizes = [c for c in dic.values()]
        title = "Nombre d'occurences de chaques base nucléique \npour les variants du chromosome " + chrom
        
        plotPieChart(sizes, labels, title)

        
# Fonction de plot de pie chart
def plotPieChart(sizes, labels, title):

    plt.figure()
    plt.pie(sizes, labels = labels, autopct='%1.1f%%', startangle=140)
    plt.title(title, bbox={'facecolor':'0.8', 'pad':5})


# Fonction de plot de diagramme
def plotBar(x, height, xlabel, ylabel, title, color, edgecolor):

    width = 1.0
    plt.figure()
    plt.bar(x, height, width, color=color, edgecolor = edgecolor)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)


def main():

    # ETAPE 1 : Vérification et ouverture du fichier
    if checkName(filePath):
        f = getFile(filePath)

    # ETAPE 2 : Stockage des informations
    fillAllDatas(f)

    #ETAPE 3 : Analyse des variants
    
    # Idées : 
    #         - Combien de délétions/insertions ?
    #         - Taille
    #         - Analyse sur les tags existants
    
    #TODO : Analyser les chromosomes en fonction des infos
    
    QualityByChrom() # Analyse en fonction de la qualité
    nbrREFTotal()    # Analyse des bases azotées les plus rencontrées
    nbrREFByChrom('1')  # Idem mais pour un chromosome en particulier
    

def printInfo():
    global COLONNES
    
    for (nom, liste) in COLONNES.items():
        print(nom , ":")
        for dicoInfo in liste:
            for (infoName, infoValue) in dicoInfo.items():
                print(infoName, ":", infoValue, end=' ')
            print('')
            

def printChromName():
    global VARIANTS

    d = set()

    for variant in VARIANTS:
        d.add(variant['CHROM'][0])

    print(d)
    return d

main()
#print(SYNCHRO_COLUMNS)
#printInfo()
#print(VARIANTS[0])
plt.show(block = False)
printChromName()
