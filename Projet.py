import sys, re

filePath = r"C:\Users\Jack\Desktop\Cours\HMIN113M_Systeme\Projet_variants\human_CEU_MEI.vcf"

date = ""
# Les cas possibles pour chaque colonne obligatoire : dictionnaire de liste de dictionnaires
# TODO rajouter les cas de base ?
COLONNES = {'CHROM' : [], 'POS' : [], 'ID' : [], 'REF' : [], 'ALT' : [],
            'QUAL' : [], 'FILTER' : [], 'INFO' : [], 'FORMAT' : []}


# TODO
# Transforme la date en une date lisible
def setDate(chaine):
    global date
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
# Un cas est initialement de la forme : <ID=X,Number=Y,Type=Z,Description=W>
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
    COLONNES[caseName].append(dic)


# Rempli les listes avec les données d'une ligne
def fillInfos(line):
    global COLONNES
    
    caseName = re.search('^[^=]*', line).group(0)

    if caseName in COLONNES.keys(): # Si c'est une information sur les champs des colonnes
        addCase(caseName, line)
        
    elif caseName == 'fileDate':    # Si c'est la date
        setDate(line[9:])

    else :  # Champs inutile ou à rajouter
        pass

        
def printInfo():
    global COLONNES
    for (nom, liste) in COLONNES.items():
        print(nom , ":")
        for dicoInfo in liste:
            for (infoName, infoValue) in dicoInfo.items():
                print(infoName, ":", infoValue, end=', ')
            print('')


def main():
    
    if(checkName(filePath)):
        f = getFile(filePath)

    for line in f : # Pour chaque ligne de f
        if line.startswith('##'):   # Lignes d'informations
            fillInfos(line[2:])
        
        elif line.startswith('#'): # Ligne représentant les colonnes
            #TODO : Rajouter les colonnes secondaires (pas toujours là)
            pass
        
        else:   # Les lignes représentant un chromosome
            #TODO : Stocker dans un dictionnaire les informations en fonction des infos (Number)
            pass
    printInfo()


    #TODO : Analyser les chromosomes en fonction des infos
        
main()
print(date)

