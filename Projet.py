import sys, re

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
        column = re.search('FORMAT.*', line).group(0)
        column_list = column.split('\t')
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

    
def main():

    # ETAPE 1 : Vérification et ouverture du fichier
    if(checkName(filePath)):
        f = getFile(filePath)

    # ETAPE 2 : Stockage des informations
    fillAllDatas(f)

    #ETAPE 3 : Analyse des variants
    
    # Idées : - Moyenne des qualités
    #         - Variant meilleure qualité
    #         - Combien de délétions/insertions ?
    #         - Chromosome le plus rencontré (A, T, C, G) ?
    #         - Taille
    #         - Analyse sur les tags existants
    
    #TODO : Analyser les chromosomes en fonction des infos
    
def printInfo():
    global COLONNES
    
    for (nom, liste) in COLONNES.items():
        print(nom , ":")
        for dicoInfo in liste:
            for (infoName, infoValue) in dicoInfo.items():
                print(infoName, ":", infoValue, end=' ')
            print('')

main()
#print(SYNCHRO_COLUMNS)
#printInfo()
