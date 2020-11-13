import sys, re
import matplotlib.pyplot as plt

plt.rcParams['toolbar'] = 'None'

filePath = r"C:\Users\Jack\Desktop\Cours\HMIN113M_Systeme\Projet_variants\HUMAN_CEU_MEI.vcf"

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


def init_colonnes():
    '''
        Initialisation des colonnes avec les balises réservées
    '''
    
    global COLONNES
    
    FILTER_INIT = [{'ID' : 'PASS', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'Le variant passe le filtre'},
                   {'ID' : '.', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'Non renseigné'}]

    INFO_INIT = [{'ID' : 'AA', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'ancestral allele'},
                 {'ID' : 'AC', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'allele count in genotypes, for each ALT allele, in the same order as listed'},
                 {'ID' : 'AF', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes'},
                 {'ID' : 'AN', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'total number of alleles in called genotypes'},
                 {'ID' : 'BQ', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'RMS base quality at this position'},
                 {'ID' : 'CIGAR', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'cigar string describing how to align an alternate allele to the reference allele'},
                 {'ID' : 'DB', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'dbSNP membership'},
                 {'ID' : 'DP', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'combined depth across samples, e.g. DP=154'},
                 {'ID' : 'END', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'end position of the variant described in this record (esp. for CNVs)'},
                 {'ID' : 'H2', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'membership in hapmap2'},
                 {'ID' : 'MQ', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'RMS mapping quality, e.g. MQ=52'},
                 {'ID' : 'MQ0', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'Number of MAPQ == 0 reads covering this record'},
                 {'ID' : 'NS', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'Number of samples with data'},
                 {'ID' : 'SB', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'strand bias at this position'},
                 {'ID' : 'SOMATIC', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'indicates that the record is a somatic mutation, for cancer genomics'},
                 {'ID' : 'VALIDATED', 'Number' : 0, 'Type' : 'FLAG', 'Description' : 'validated by follow-up experiment'}]


    # Ajout des valeurs si elles n'y sont pas déjà 

    for new_case in FILTER_INIT :
        isIn = False
        ID = new_case['ID']
        for existing_case in COLONNES['FILTER']:
            if ID == existing_case['ID']:
                isIn = True
                break
        if not isIn :
            COLONNES['FILTER'].append(new_case)
        
    for new_case in INFO_INIT :
        isIn = False
        ID = new_case['ID']
        for existing_case in COLONNES['INFO']:
            if ID == existing_case['ID']:
                isIn = True
                break
        if not isIn :
            COLONNES['INFO'].append(new_case)


# TODO
def setDate(chaine):
    '''
        Transforme la date en une date lisible
    '''
    
    global DATE
    
    date = chaine


def checkName(filePath):
    '''
        Vérifie si le fichier est bien un .vcf
    '''
    
    return filePath.endswith(".vcf")


def getFile(filePath):
    '''
        Ouvre et renvoie le fichier 
    '''
    
    try :
        f = open(filePath)
        return f
    except :
        print("Erreur lors de l'ouverture du fichier")
        sys.exit(1)

        
def addCase(caseName, case):
    '''
         Ajoute au dictionnaire des cas un cas donné,  un cas est initialement de la forme : <ID=X,Number=Y,Type=Z,Description=W, ... >
         Il est transformé en la forme : {ID : X, Number : Y, Type : Z, Description : W}
         Champs par défaut si manquants Type -> FLAG, Number -> 0, Description -> ''
    '''
    
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


def fillInfos(line):
    '''
        Rempli les listes avec les données d'une ligne
    '''
    
    global COLONNES
    
    caseName = re.search('^[^=]*', line).group(0)

    if caseName in COLONNES.keys(): # Si c'est une information sur les champs des colonnes
        addCase(caseName, line)
        
    elif caseName == 'fileDate':    # Si c'est la date
        setDate(line[9:])

    #TODO : Rajouter des cas : version, reference, source, contig ?
        
    else :  # Champs non utilisés
        pass


def fillAdditionnalColumns(line):
    '''
        Rempli l'ensemble des colonnes supplémentaires et synchronise leur ordre
    '''
    
    global ADDITIONAL_COLUMNS, SYNCHRO_COLUMNS
    
    try :
        column = re.search('FORMAT.*', line).group(0)   # Nom des colonnes apres INFO
        column_list = column.split('\t')    # Transformation en liste
        
        if ' ' in column_list[0]:   # Si les colonnes sont séparées par des espaces et non des tabulations
            column_list = column.split(' ')
            
        for column_name in column_list:
            if column_name != '' :
                column_name.replace('\n', '') # Supression des sauts de ligne
                ADDITIONAL_COLUMNS.add(column_name) # Ajout de la colonne à l'ensemble des colonnes
                if column_name not in SYNCHRO_COLUMNS.values(): # Synchronisation
                    SYNCHRO_COLUMNS[len(SYNCHRO_COLUMNS)] = column_name 
            
    except :    # Si il n'y a pas de colonnes supplémentaires
        pass


def fillVariants(line):
    '''
        Rempli la liste de variants avec les informations d'une ligne (1 variant)
        Synchronise le dictionnaire des colonnes
    '''
    
    global SYNCHRO_COLUMNS, VARIANTS
    
    dic = {}
    field_list = line.split('\t')

    if len(field_list) < 8 :    # Si les colonnes sont séparées par des espaces et non des tabulations
        field_list = line.split(' ')
    
    column_number = 0
    for column in field_list :  # Pour chaque champ
        if column != '' :
            column = column.replace('\n', '') # Supression des sauts de ligne
            elements_list = column.split(';')
            dic[SYNCHRO_COLUMNS[column_number]] = elements_list
            column_number += 1

    VARIANTS.append(dic)


def fillAllDatas(f):
    '''
        Rempli les dictionnaires et listes avec les données en entrée
    '''
    
    for line in f : # Pour chaque ligne de f
        if line.startswith('##'):   # Lignes d'informations
            fillInfos(line[2:])
        
        elif line.startswith('#'): # Ligne représentant les colonnes
            fillAdditionnalColumns(line)
        
        else:   # Les lignes représentant un variant
            fillVariants(line)
            pass


######################################################################################################################################################################
######################################################################################################################################################################
                                                # Fonctions d'analyse

def getInsertionsAndDeletions(chaine1, chaine2):
    '''
        Retourne les insertions et délétions lors du passage de chaine1 à chaine2
        Utilisé dans l'analyse des champs REF et ALT
    '''

    dic = {'insertions' : [], 'deletions' : []}
    l1 = [base for base in chaine1]
    l2 = [base for base in chaine2]

    while l1 :  # Tant que l1 n'est pas vide
        base = l1[0]
        
        if base in l2 :
            l2.remove(base)
        else :
            dic['deletions'].append(base)
            
        l1 = l1[1:]

    dic['insertions'].extend(l2)    # Ajout des bases non trouvées

    return dic


def isMadeOfNuc(chaine):
    '''
        Renvoie True si une chaine de caractères est composée
        de nucléotides seulement (ATCGN)
        Utilisé dans l'analyse des champs REF et ALT
    '''

    return bool(re.match('^[ATCGN]*$', chaine))


def insertAndDel2Dic(dicID):
    '''
        Transforme un dictionnaire {'insertions' : [], 'deletions' : []},
        en 2 dictionnaires avec les occurences pour chaque nucléotide
        Utilisé dans l'analyse des champs REF et ALT
    '''

    dic_ins = {}
    dic_del = {}

    for ins in list(dicID.values())[0]:
        try :
            dic_ins[ins] += 1
        except :
            dic_ins[ins] = 1

    for dele in list(dicID.values())[1]:
        try :
            dic_del[dele] += 1
        except :
            dic_del[dele] = 1

    return dic_ins, dic_del


def mergeDic(dic1, dic2):
    '''
        Ajoute le contenu de dic2 à dic1, crée une nouvelle clé
        si elle n'exsitait pas dans dic1
        Utilisé dans l'analyse des champs REF et ALT
    '''

    for key, value in dic2.items():
        try :
            dic1[key] = dic1[key] + value
        except :
            dic1[key] = value

    return dic1



def descriptionBarString(dic):
    '''
        Transforme un dictionnaire {ID : desc, ID : desc, ... }
        en un String -> 'ID : desc \nID : desc ... }
        Utilisé dans l'analyse des champs REF et ALT
    '''
    
    res = ''
    for ID, desc in dic.items():
        if desc :
            res += ID + ' : ' + desc + '\n'
    res = res[:-1]
    return res


def analyzeInsertAndDelTotal():
    '''
        Analyze les insertions et délétions pour tous les variants de
        chaque chromosome, en se basant sur les colonnes REF et ALT
        Affichage de 3 plot
    '''

    global VARIANTS, COLONNES

    dic_ins = {}    # Forme {A : 1, C : 14, N : 7, ... }
    dic_del = {}    # Idem
    dic_balises = {}    # Forme {ID_Balise1 : 7, ID_Balise2 : 4, ... }
    dic_descriptions = {}
    
    for variant in VARIANTS:
        ref = variant['REF'][0] # Base de référence
        ALT = variant['ALT']    # Liste des bases alternatives

        for alt in ALT :
            alt_split = alt.split(',')  # Pour les cas : 'T,C' -> ['T', 'C']
            for new_alt in alt_split:
                if isMadeOfNuc(new_alt) and isMadeOfNuc(ref):  # Si ce sont des bases nucléiques
                    dic_temp = getInsertionsAndDeletions(ref, new_alt)
                    dic_ins_temp, dic_del_temp = insertAndDel2Dic(dic_temp)
                    dic_ins = mergeDic(dic_ins, dic_ins_temp)   # Mise à jour des dictionnaires
                    dic_del = mergeDic(dic_del, dic_del_temp)

                else :  # Si c'est une balise
                    for case in COLONNES['ALT']:
                        new_alt = new_alt.replace('>', '')
                        new_alt = new_alt.replace('<', '')
                        if case['ID'] == new_alt :  # Si elle est répertoriée
                            try :
                                dic_balises[new_alt] += 1   # Ajout au dictionnaire des balises
                            except :
                                dic_balises[new_alt] = 1    

                            description = case['Description']
                            if (new_alt, description) not in dic_descriptions.items():  # Ajout aux descriptions pour légender
                                dic_descriptions[new_alt] = description
    
    # Plot des 3 dictionnaires si jugé nécessaire

    if sum(dic_ins.values()) > 0 :
            
        x_ins = dic_ins.keys()
        height_ins = dic_ins.values()
        ylabel_ins = "Nombre d'insertions"
        xlabel_ins = 'Nucléotide'
        title_ins = 'Insertions des nucléotides pour tous les variants'
        
        color_ins = 'green'
        edgecolor_ins = 'blue'
        
        plotBar(x_ins, height_ins, xlabel_ins, ylabel_ins, title_ins, color_ins, edgecolor_ins)

    if sum(dic_del.values()) > 0 :
        x_del = dic_del.keys()
        height_del = dic_del.values()
        ylabel_del = 'Nombre de délétions'
        xlabel_del = 'Nucléotides'
        title_del = 'Délétions des nucléotides pour tous les variants'
        color_del = 'red'
        edgecolor_del = 'yellow'
        
        plotBar(x_del, height_del, xlabel_del, ylabel_del, title_del, color_del, edgecolor_del)
                
    if sum(dic_balises.values()) > 0 :
        x_balises = dic_balises.keys()
        height_balises = dic_balises.values()
        ylabel_balises = 'Nombre de délétions et insertions'
        xlabel_balises = 'Nucléotides'
        title_balises = 'Délétions et insertions des nucléotides pour tous les variants'
        color_balises = 'white'
        edgecolor_balises = 'orange'
        text_balises = descriptionBarString(dic_descriptions)
        
        plotBar(x_balises, height_balises, xlabel_balises, ylabel_balises, title_balises, color_balises, edgecolor_balises, text_balises)


def analyzeInsertAndDelChrom(chrom):
    '''
        Analyze les insertions et délétions pour tous les variants
        d'un chromosome donné, en se basant sur les colonnes REF et ALT
        Affichage de 3 plot
    '''
    
    global VARIANTS, COLONNES

    dic_ins = {}    # Forme {A : 1, C : 14, N : 7, ... }
    dic_del = {}    # Idem
    dic_balises = {}    # Forme {ID_Balise1 : 7, ID_Balise2 : 4, ... }
    dic_descriptions = {}

    for variant in VARIANTS:
        if variant['CHROM'][0] == chrom :  # Selection pour le bon chromosome
            ref = variant['REF'][0] # Base de référence
            ALT = variant['ALT']    # Liste des bases alternatives

            for alt in ALT :    # Pour chaque élément de la colonne ALT
                alt_split = alt.split(',')  # Pour les cas : 'T,C' -> ['T', 'C']
                for new_alt in alt_split:
                    if isMadeOfNuc(new_alt) and isMadeOfNuc(ref):  # Si ce sont des bases nucléiques
                        dic_temp = getInsertionsAndDeletions(ref, new_alt)
                        dic_ins_temp, dic_del_temp = insertAndDel2Dic(dic_temp)
                        dic_ins = mergeDic(dic_ins, dic_ins_temp)   # Mise à jour des dictionnaires
                        dic_del = mergeDic(dic_del, dic_del_temp)

                    else :  # Si c'est une balise
                        for case in COLONNES['ALT']:
                            new_alt = new_alt.replace('>', '')
                            new_alt = new_alt.replace('<', '')
                            if case['ID'] == new_alt :  # Si elle est répertoriée
                                try :
                                    dic_balises[new_alt] += 1   # Ajout au dictionnaire des balises
                                except :
                                    dic_balises[new_alt] = 1    

                                description = case['Description']
                                if (new_alt, description) not in dic_descriptions.items():  # Ajout aux descriptions pour légender
                                    dic_descriptions[new_alt] = description
        
    # Plot des 3 dictionnaires si jugé nécessaire

    if sum(dic_ins.values()) > 0 :
            
        x_ins = dic_ins.keys()
        height_ins = dic_ins.values()
        ylabel_ins = "Nombre d'insertions"
        xlabel_ins = 'Nucléotide'
        title_ins = 'Insertions des nucléotides pour tous les variants du chromosome ' + str(chrom)
        
        color_ins = 'green'
        edgecolor_ins = 'blue'
        
        plotBar(x_ins, height_ins, xlabel_ins, ylabel_ins, title_ins, color_ins, edgecolor_ins)

    if sum(dic_del.values()) > 0 :
        x_del = dic_del.keys()
        height_del = dic_del.values()
        ylabel_del = 'Nombre de délétions'
        xlabel_del = 'Nucléotides'
        title_del = 'Délétions des nucléotides pour tous les variants du chromosome ' + str(chrom)
        color_del = 'red'
        edgecolor_del = 'yellow'
        
        plotBar(x_del, height_del, xlabel_del, ylabel_del, title_del, color_del, edgecolor_del)
                
    if sum(dic_balises.values()) > 0 :
        x_balises = dic_balises.keys()
        height_balises = dic_balises.values()
        ylabel_balises = 'Nombre de délétions et insertions'
        xlabel_balises = 'Nucléotides'
        title_balises = 'Délétions et insertions des nucléotides pour les variants du chromosome ' + str(chrom)
        color_balises = 'white'
        edgecolor_balises = 'orange'
        text_balises = descriptionBarString(dic_descriptions)
        
        plotBar(x_balises, height_balises, xlabel_balises, ylabel_balises, title_balises, color_balises, edgecolor_balises, text_balises)
    

def qualityTotal():
    '''
        Analyse de la qualité des variants pour tous les chromosomes, triés par chromosome
        Affichage de barplot
    '''
    
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

    if maxQual > 0 :    # Si tous les variants n'ont pas une qualité de 0

        x = dic.keys()
        height = [sum(qual)/len(qual) for qual in dic.values()]
        ylabel = 'Qualité'
        xlabel = 'Chromosome'
        title = 'Moyenne de la qualité des variants pour chaque chromosome'
        color = (0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0)
        edgecolor = 'blue'
        
        plotBar(x, height, xlabel, ylabel, title, color, edgecolor)


def qualityChrom(chrom):
    '''
        Analyse de la qualité des variants pour un chromosome donné
        Affichage de barplot
    '''
    
    global VARIANTS

    dic = {}    # Forme {1 : qualité, 2 : qualité, ... }
    i = 1
    for variant in VARIANTS :
        if variant['CHROM'][0] == chrom :
            try :
                dic[i] = int(variant['QUAL'][0])
            except :
                dic[i] = 0
            i += 1

    if(sum(dic.values()) > 0) : # Si tous les variants n'ont pas une qualité de 0
        
        x = dic.keys()
        height = dic.values()
        ylabel = 'Qualité'
        xlabel = 'Variant'
        title = 'Qualité des variants du chromosome ' + str(chrom)
        color = (0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0)
        edgecolor = 'black'
        
        plotBar(x, height, xlabel, ylabel, title, color, edgecolor)
    

def REFbaseTypeTotal():
    '''
        Analyse des ccurences des base nucléiques de l'allèle de référence
        Affichage de pieChart
    '''
    
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
        sizes = dic.values()
        title = "Nombre d'occurences de chaques base nucléique du génome de référence"
        
        plotPieChart(sizes, labels, title)


def REFbaseTypeChrom(chrom):
    '''
        Occurence des bases nucléiques de l'allèle de référence pour un chromosome donné
        Affichage de pieChart
    '''
    
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
        title = "Nombre d'occurences de chaques base nucléique \ndu génome de référence pour le chromosome " + chrom
        
        plotPieChart(sizes, labels, title)


def analyzeFilterTotal():
    '''
        Analyse du pourcentage de variants qui ont passé le filtre, trié par chromosome
        Affichage de pieChart
    '''
    
    global VARIANTS, COLONNES

    dic  = {}
    dic_legends = {}    # {Nom du (des) filtre utilisé : description}
    
    for variant in VARIANTS :

        FILTER = variant['FILTER']  # Qualité du variant    
        name = ''
        
        for filter_name in FILTER : # Première boucle pour créer le nom
            name += filter_name + ' et '
        name = name[0:-4]   # Enlève le 'et' à la fin

        if name not in dic_legends.keys():
            
            for filter_name in FILTER : # Seconde boucle pour les descriptions
                
                try:
                    dic_legends[name] += getFilterDescription(filter_name) + '\n'
                except :
                    dic_legends[name] = getFilterDescription(filter_name) + '\n'

            dic_legends[name] = dic_legends[name][:-1]  # Enlève le '\n' à la fin
        
        try :
            dic[name] += 1
        except :
            dic[name] = 1

    labels = dic.keys()
    sizes = dic.values()
    title = 'Filtre des variants'
    legends = dic_legends.values()
    
    plotPieChart(sizes, labels, title, legends)


def getFilterDescription(ID):
    '''
        Renvoie la description d'un FILTER donné
        Utilisé dans l'analyse des FILTER
    '''
    
    global COLONNES
    
    for FILTER in COLONNES['FILTER']:
        if FILTER['ID'] == ID :
            return FILTER['Description']


def analyzeFilterChrom(chrom):
    '''
        Analyse du filtrage des variants pour un chromosome donné
    '''
    
    global VARIANTS, COLONNES

    dic  = {}
    dic_legends = {}    # {Nom du (des) filtre utilisé : description}
    
    for variant in VARIANTS :
        if variant['CHROM'][0] == chrom :

            FILTER = variant['FILTER']  # Qualité du variant    
            name = ''
            
            for filter_name in FILTER : # Première boucle pour créer le nom
                name += filter_name + ' et '
            name = name[0:-4]   # Enlève le 'et' à la fin

            if name not in dic_legends.keys():
                
                for filter_name in FILTER : # Seconde boucle pour les descriptions
                    
                    try:
                        dic_legends[name] += getFilterDescription(filter_name) + '\n'
                    except :
                        dic_legends[name] = getFilterDescription(filter_name) + '\n'

                dic_legends[name] = dic_legends[name][:-1]  # Enlève le '\n' à la fin
            
            try :
                dic[name] += 1
            except :
                dic[name] = 1

    labels = dic.keys()
    sizes = dic.values()
    title = 'Filtre pour les variants du chromosome ' + str(chrom)
    legends = dic_legends.values()
    
    plotPieChart(sizes, labels, title, legends)


def plotPieChart(sizes, labels, title, legends = []):
    '''
        Fonction de plot de pie chart
    '''
    
    plt.figure(figsize=(8,5))
    patches, texts, junk = plt.pie(sizes, labels = labels, autopct='%1.1f%%', startangle=180)
    plt.title(title, bbox={'facecolor':'0.8', 'pad':5})

    if legends :
        plt.legend(patches, legends)


def plotBar(x, height, xlabel, ylabel, title, color, edgecolor, text = ''):
    '''
        Fonction de plot de diagramme
    '''

    width = 1.0
    fig, ax = plt.subplots(figsize=(7,5))
    plt.bar(x, height, width, color=color, edgecolor = edgecolor)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)

    if text :
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=props)


def main():
    global VARIANT, SYNCHRO_COLUMNS, COLONNES
    
    # ETAPE 1 : Vérification et ouverture du fichier
    if checkName(filePath):
        f = getFile(filePath)

    # ETAPE 2 : Stockage des informations
    fillAllDatas(f)
    init_colonnes()

    
    #printInfo()
    #print(COLONNES)
    #printVariantFILTER()
    #print(SYNCHRO_COLUMNS)
    #printVariants()
    

    #ETAPE 3 : Analyse des variants
    
    # Idées : 
    #         - Taille
    #         - Analyse sur les tags existants
    
    #TODO : Analyser les chromosomes en fonction des infos

    qualityTotal()  # Analyse de la qualité pour tous les variants de tous les chromosomes
    #qualityChrom('1') # Analyse de la qualité des variants d'un chromosome donné
    
    REFbaseTypeTotal()    # Analyse des bases azotées les plus rencontrées pour l'allèle de référence
    #REFbaseTypeChrom('1')  # Idem mais pour un chromosome en particulier pour l'allèle de référence

    analyzeFilterTotal()   # Analyse du filtre pour tous les variants
    #analyzeFilterChrom('1')    # Analyse du filtre pour un chromosome donné
    
    analyzeInsertAndDelTotal()  # Analyse des insertions et délétions
    #analyzeInsertAndDelChrom('1')  # Analyse des insertions et délétions pour un chromosome donné

##    for chrom in printChromName():
##        qualityChrom(chrom)
##        analyzeFilterChrom(chrom)
##        REFbaseTypeChrom(chrom)
##        analyzeInsertAndDelChrom(chrom)
    

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

def printVariantFILTER():
    global VARIANTS
    
    for variant in VARIANTS:
        print(variant['FILTER'])

def printVariants():
    for variant in VARIANTS :
        print(variant)

main()

plt.show(block = False)
