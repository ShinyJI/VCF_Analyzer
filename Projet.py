import os
import re
import sys

import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import *
from tkinter import ttk, filedialog, messagebox
from PIL import Image, ImageTk
from random import choice


plt.rcParams['toolbar'] = 'None'

# Les cas possibles pour chaque colonne obligatoire : dictionnaire de liste de dictionnaires
COLONNES = {
    'CHROM': [],
    'POS': [],
    'ID': [],
    'REF': [],
    'ALT': [],
    'QUAL': [],
    'FILTER': [],
    'INFO': [],
    'FORMAT': [],
}

# Les nom des colonnes facultatives stockés dans un ensemble (pas de doublons)
ADDITIONAL_COLUMNS = set()

# Les informations sur les variants stockées dans une liste de dictionnaires
# De la forme : VARIANTS = [{CHROM : [v1.1, v1.2], POS : [v2.1], ...}, {CHROM: [v2.1]}, ... ]
VARIANTS = list()

# Dictionnaire synchronisant les colonnes avec leur ordre de lecture (Equivalent d'Hashmap)
SYNCHRO_COLUMNS = {
    0: 'CHROM',
    1: 'POS',
    2: 'ID',
    3: 'REF',
    4: 'ALT',
    5: 'QUAL',
    6: 'FILTER',
    7: 'INFO',
}


############## INTERFACE GRAPHIQUE ##############


selected_analyse = ''
selected_chrom = ''
selected_quality = 0
filename = ''


def getChromList():
    """
    Retourne la liste des chromosomes
    """

    global VARIANTS

    L = []

    for variant in VARIANTS:
        chrom = variant['CHROM'][0]
        if chrom not in L:
            L.append(chrom)
    return L


# Création de la fenêtre
root = Tk()
root.title('VFC Analizer')
root.geometry('800x500')
root.resizable(False, False)
root.iconbitmap(os.getcwd() + r'\icone_vcf.ico')


# Création des frames
frame_1 = Frame(root, width=800, height=200)  # Image
frame_1.grid(row=0, column=0, pady=10)
image_path = os.getcwd() + r'\image_adn.png'


frame_2 = Frame(root, width=800, height=100)  # Browse Button + current_file_name
frame_2.grid(row=1, column=0, pady=15)
shown_file_name = Label(frame_2, text='')

frame_3 = Frame(
    root, width=800, height=100
)  # 'Analyse' + 'chrom' + 2 dropdowns + 'Qualité minimale...' + dropdown
frame_3.grid(row=2, column=0)

frame_4 = Frame(root, width=800, height=100)  # Analyse Button + AnalyseRandom Button
frame_4.grid(row=3, column=0)

color = 'pink'
root.configure(bg=color)
frame_2.configure(bg=color)
frame_3.configure(bg=color)
frame_4.configure(bg=color)


# Fermeture du programme si on ferme la fenêtre
def on_closing():
    root.destroy()
    
root.protocol("WM_DELETE_WINDOW", on_closing)

def makeWindow():
    """
    Instancie les éléments de la fenêtre graphique
    """

    global root

    makeImage()

    makeBrowseButton()
    UpdateShownFilename()

    makeAnalyseText()
    makeChromText()

    makeAnalyseButton()
    makeAnalyseRandomButton()


def makeImage():

    global frame_1, image_path

    image = Image.open(image_path)
    image = image.resize((796, 200), Image.ANTIALIAS)
    photo = ImageTk.PhotoImage(image)

    label = Label(frame_1, image=photo)
    label.image = photo
    label.pack()


def makeAnalyseText():

    global frame_3

    text = 'Analyse : '

    analyseText = Label(frame_3, text=text)
    analyseText.grid(row=0, column=0, pady=20)


def makeChromText():

    global frame_3

    text = 'Chromosome : '

    analyseText = Label(frame_3, text=text)
    analyseText.grid(row=0, column=2, pady=20)


# Event lors de la sélection d'un chromosome
def updateSelectedChrom(event):
    global selected_chrom, dropdown_chrom

    selected_chrom = dropdown_chrom.get()


CHROM_LIST = ['']
dropdown_chrom = ttk.Combobox(frame_3, value=CHROM_LIST, state="readonly")
dropdown_chrom.current(0)
dropdown_chrom.bind('<<ComboboxSelected>>', updateSelectedChrom)
dropdown_chrom.grid(column=3, row=0, pady=20)


# Event lors de la sélection d'une analyse
def updateSelectedAnalyse(event):
    global selected_analyse, dropdown_analyse

    selected_analyse = dropdown_analyse.get()


ANALYSES_LIST = ['', 'Qualité', 'Génome de référence', 'Filtre', 'Insertions/Délétions']
dropdown_analyse = ttk.Combobox(frame_3, value=ANALYSES_LIST, state="readonly")
dropdown_analyse.current(0)
dropdown_analyse.bind('<<ComboboxSelected>>', updateSelectedAnalyse)
dropdown_analyse.grid(column=1, row=0, pady=20)


def makeQualityText():
    global frame_3

    text = 'Analyse sur les variants de qualité supérieure ou égale à  : '

    analyseText = Label(frame_3, text=text)
    analyseText.grid(row=1, column=1, padx=10)


# Event lors de la sélection d'une qualité minimale
def updateSelectedQuality(event):
    global selected_quality

    selected_quality = dropdown_quality.get()


QUALITY_LIST = [i for i in range(0, 101)]
dropdown_quality = ttk.Combobox(frame_3, value=QUALITY_LIST, state="readonly")
dropdown_quality.current(0)
dropdown_quality.grid(column=2, row=1)
dropdown_quality.bind('<<ComboboxSelected>>', updateSelectedQuality)
makeQualityText()


# Event du bouton d'analyse random
def clickAnalyseRandom():
    """
    Lance une analyse au hasard avec un chromosome au hasard
    """
    global filename, analyseToFunction

    if not filename:
        openNoFileSelectedWindow()
        return

    random_chrom = choice(CHROM_LIST)
    random_analyse = choice(ANALYSES_LIST[1:])

    if random_chrom == 'GLOBAL':
        analyseToFunction[random_analyse][0]()

    else:
        analyseToFunction[random_analyse][1](random_chrom)


def makeAnalyseRandomButton():
    """
    Crée le bouton d'une analyse random
    """

    global frame_4

    buttonAnalyseRandom = Button(
        frame_4, text='Analyse Random', command=clickAnalyseRandom
    )
    buttonAnalyseRandom.grid(column=0, row=0, padx=5, pady=20)


def clickAnalyse():
    """
    Evènement lors du click sur Analyse
    Lance l'analyse désirée si elle est possible
    """

    global selected_analyse, selected_chrom, analyseToFunction, filename

    if not filename:
        openNoFileSelectedWindow()

    elif selected_analyse:  # Si une analyse est sélectionnée

        if (
            selected_chrom in getChromList()
        ):  # Si un chromosome est selectionné, analyse par chromosome
            analyseToFunction[selected_analyse][1](selected_chrom)

        else:  # Analyse globale
            analyseToFunction[selected_analyse][0]()

    else:
        openNoAnalyseSelectedWindow()


def makeAnalyseButton():
    """
    Crée le bouton d'une analyse random
    """

    global frame_4

    buttonAnalyse = Button(frame_4, text='Analyse', command=clickAnalyse)
    buttonAnalyse.grid(column=1, row=0, padx=5, pady=20)


def clickBrowse():
    """
    Cherche un fichier, le vérifie et l'ouvre
    """

    global filename

    choice = filedialog.askopenfilename(
        initialdir="/",
        title="Select A File",
        filetype=(("vcf files", "*.vcf"), ("all files", "*.*")),
    )

    if not choice:  # Pas de fichié selectionné
        pass

    elif checkName(choice):  # Si le nom du fichier est correct
        file = getFile(choice)
        if file:  # Si le fichier s'ouvre correctement
            setGlobalEmpty()  # Réinitialise
            init_all(file)  # Rempli
            updateDropDown()  # Met à jour
            filename = choice
            UpdateShownFilename()

    else:
        openFichierIncorrectWindow()


def getNameAndExt(filename):
    """
    Renvoie le nom et l'extension d'une fichier à partir de son chemin
    """

    return os.path.basename(filename)


def UpdateShownFilename():
    """
    Remplace et affiche le nom du fichier choisi sur la fenêtre
    """

    global frame_2, shown_file_name, filename, color

    text = 'Fichier : ' + getNameAndExt(filename)

    shown_file_name.destroy()  # Enlève l'ancien affichage

    shown_file_name = Label(
        frame_2, text=text, font="Time 12 bold", bg=color
    )  # Crée et affiche un nouveau
    shown_file_name.pack(padx=5, pady=10)


def makeBrowseButton():
    """
    Crée le bouton de recherche de fichier
    """

    global frame_2

    myButtonBrowse = Button(
        frame_2, text='Sélectionner un fichier', command=clickBrowse
    )
    myButtonBrowse.pack(padx=5, pady=10)


def updateDropDown():
    """
    Mise à jour du dropdown.
    Crée un nouveau bouton au même endroit au lieu de remplacer les valeurs,
    pour des raisons esthétiques
    """

    global dropdown_chrom, CHROM_LIST, frame_3

    CHROM_LIST = getChromList()
    CHROM_LIST.insert(0, 'GLOBAL')
    dropdown_chrom = ttk.Combobox(frame_3, value=CHROM_LIST, state="readonly")
    dropdown_chrom.current(0)
    dropdown_chrom.bind('<<ComboboxSelected>>', updateSelectedChrom)
    dropdown_chrom.grid(column=3, row=0)


############## INITIALISATION ET STOCKAGE ##############


def init_colonnes():
    """
    Initialisation des colonnes avec les balises réservées
    """

    global COLONNES

    FILTER_INIT = [
        {
            'ID': 'PASS',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'Le variant passe le filtre',
        },
        {'ID': '.', 'Number': 0, 'Type': 'FLAG', 'Description': 'Non renseigné'},
    ]

    INFO_INIT = [
        {'ID': 'AA', 'Number': 0, 'Type': 'FLAG', 'Description': 'ancestral allele'},
        {
            'ID': 'AC',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'allele count in genotypes, for each ALT allele, in the same order as listed',
        },
        {
            'ID': 'AF',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes',
        },
        {
            'ID': 'AN',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'total number of alleles in called genotypes',
        },
        {
            'ID': 'BQ',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'RMS base quality at this position',
        },
        {
            'ID': 'CIGAR',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'cigar string describing how to align an alternate allele to the reference allele',
        },
        {'ID': 'DB', 'Number': 0, 'Type': 'FLAG', 'Description': 'dbSNP membership'},
        {
            'ID': 'DP',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'combined depth across samples, e.g. DP=154',
        },
        {
            'ID': 'END',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'end position of the variant described in this record (esp. for CNVs)',
        },
        {
            'ID': 'H2',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'membership in hapmap2',
        },
        {
            'ID': 'MQ',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'RMS mapping quality, e.g. MQ=52',
        },
        {
            'ID': 'MQ0',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'Number of MAPQ == 0 reads covering this record',
        },
        {
            'ID': 'NS',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'Number of samples with data',
        },
        {
            'ID': 'SB',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'strand bias at this position',
        },
        {
            'ID': 'SOMATIC',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'indicates that the record is a somatic mutation, for cancer genomics',
        },
        {
            'ID': 'VALIDATED',
            'Number': 0,
            'Type': 'FLAG',
            'Description': 'validated by follow-up experiment',
        },
    ]

    # Ajout des valeurs si elles n'y sont pas déjà

    for new_case in FILTER_INIT:
        isIn = False
        ID = new_case['ID']
        for existing_case in COLONNES['FILTER']:
            if ID == existing_case['ID']:
                isIn = True
                break
        if not isIn:
            COLONNES['FILTER'].append(new_case)

    for new_case in INFO_INIT:
        isIn = False
        ID = new_case['ID']
        for existing_case in COLONNES['INFO']:
            if ID == existing_case['ID']:
                isIn = True
                break
        if not isIn:
            COLONNES['INFO'].append(new_case)


def setGlobalEmpty():
    """
    Vide les dictionnaires et listes globaux
    Fonction utilisée lors de la sélection d'un nouveau fichier VCF
    """

    global COLONNES, ADDITIONAL_COLUMNS, VARIANTS, SYNCHRO_COLUMNS, selected_chrom

    COLONNES = {
        'CHROM': [],
        'POS': [],
        'ID': [],
        'REF': [],
        'ALT': [],
        'QUAL': [],
        'FILTER': [],
        'INFO': [],
        'FORMAT': [],
    }

    ADDITIONAL_COLUMNS = set()

    VARIANTS = list()

    SYNCHRO_COLUMNS = {
        0: 'CHROM',
        1: 'POS',
        2: 'ID',
        3: 'REF',
        4: 'ALT',
        5: 'QUAL',
        6: 'FILTER',
        7: 'INFO',
    }

    selected_chrom = ''


def init_all(file):
    """
    Initialise les listes et dictionnaires avec les informations du fichier
    """

    fillAllDatas(file)
    init_colonnes()


def checkName(filePath):
    """
    Vérifie si le fichier est bien un .vcf
    """

    return filePath.endswith(".vcf")


def getFile(filePath):
    """
    Ouvre et renvoie le fichier
    """

    try:
        f = open(filePath)
        return f
    except:
        print("Erreur lors de l'ouverture du fichier")
        return null


def addCase(caseName, case):
    """
    Ajoute au dictionnaire des cas un cas donné,  un cas est initialement de la forme : <ID=X,Number=Y,Type=Z,Description=W, ... >
    Il est transformé en la forme : {ID : X, Number : Y, Type : Z, Description : W}
    Champs par défaut si manquants Type -> FLAG, Number -> 0, Description -> ''
    """

    global COLONNES

    ID = re.search('ID=[^,>]*', case).group(0)[3:]

    try:
        Description = re.search('Description="[^,"]*', case).group(0)[13:]
    except:
        Description = ''

    try:
        Number = re.search('Number=[^,>]*', case).group(0)[7:]
    except:
        Number = 0

    try:
        Type = re.search('Type=[^,>]*', case).group(0)[5:]
    except:
        Type = 'FLAG'

    dic = {'ID': ID, 'Number': Number, 'Type': Type, 'Description': Description}

    if dic not in COLONNES[caseName]:  # Vérification doublons :
        COLONNES[caseName].append(
            dic
        )  # Dans certains fichiers, les en-têtes sont répétées


def fillInfos(line):
    """
    Rempli les listes avec les données d'une ligne
    """

    global COLONNES

    caseName = re.search('^[^=]*', line).group(0)

    if (
        caseName in COLONNES.keys()
    ):  # Si c'est une information sur les champs des colonnes
        addCase(caseName, line)

    else:  # Champs non utilisés
        pass


def fillAdditionnalColumns(line):
    """
    Rempli l'ensemble des colonnes supplémentaires et synchronise leur ordre
    """

    global ADDITIONAL_COLUMNS, SYNCHRO_COLUMNS

    try:
        column = re.search('FORMAT.*', line).group(0)  # Nom des colonnes apres INFO
        column_list = column.split('\t')  # Transformation en liste

        if (
            ' ' in column_list[0]
        ):  # Si les colonnes sont séparées par des espaces et non des tabulations
            column_list = column.split(' ')

        for column_name in column_list:
            if column_name != '':
                column_name.replace('\n', '')  # Supression des sauts de ligne
                ADDITIONAL_COLUMNS.add(
                    column_name
                )  # Ajout de la colonne à l'ensemble des colonnes
                if column_name not in SYNCHRO_COLUMNS.values():  # Synchronisation
                    SYNCHRO_COLUMNS[len(SYNCHRO_COLUMNS)] = column_name

    except:  # Si il n'y a pas de colonnes supplémentaires
        pass


def fillVariants(line):
    """
    Rempli la liste de variants avec les informations d'une ligne (1 variant)
    Synchronise le dictionnaire des colonnes
    """

    global SYNCHRO_COLUMNS, VARIANTS

    dic = {}
    field_list = line.split('\t')

    if (
        len(field_list) < 8
    ):  # Si les colonnes sont séparées par des espaces et non des tabulations
        field_list = line.split(' ')

    column_number = 0
    for column in field_list:  # Pour chaque champ
        if column != '':
            column = column.replace('\n', '')  # Supression des sauts de ligne
            elements_list = column.split(';')
            dic[SYNCHRO_COLUMNS[column_number]] = elements_list
            column_number += 1

    VARIANTS.append(dic)


def fillAllDatas(f):
    """
    Rempli les dictionnaires et listes avec les données en entrée
    """

    for line in f:  # Pour chaque ligne de f
        if line.startswith('##'):  # Lignes d'informations
            fillInfos(line[2:])

        elif line.startswith('#'):  # Ligne représentant les colonnes
            fillAdditionnalColumns(line)

        else:  # Les lignes représentant un variant
            fillVariants(line)
            pass


############## FONCTIONS D'ANALYSE ##############


def getQuality(variant):
    """
    Renvoie la qualité d'un variant,
    Renvoie 0 si ce n'est pas un int
    """

    try:
        qual = int(variant['QUAL'][0])
    except:
        qual = 0

    return qual


def getInsertionsAndDeletions(chaine1, chaine2):
    """
    Retourne les insertions et délétions lors du passage de chaine1 à chaine2
    Utilisé dans l'analyse des champs REF et ALT
    """

    dic = {'insertions': [], 'deletions': []}
    l1 = [base for base in chaine1]
    l2 = [base for base in chaine2]

    while l1:  # Tant que l1 n'est pas vide
        base = l1[0]

        if base in l2:
            l2.remove(base)
        else:
            dic['deletions'].append(base)

        l1 = l1[1:]

    dic['insertions'].extend(l2)  # Ajout des bases non trouvées

    return dic


def isMadeOfNuc(chaine):
    """
    Renvoie True si une chaine de caractères est composée
    de nucléotides seulement (ATCGN)
    Utilisé dans l'analyse des champs REF et ALT
    """

    return bool(re.match('^[ATCGN]*$', chaine))


def insertAndDel2Dic(dicID):
    """
    Transforme un dictionnaire {'insertions' : [], 'deletions' : []},
    en 2 dictionnaires avec les occurences pour chaque nucléotide
    Utilisé dans l'analyse des champs REF et ALT
    """

    dic_ins = {}
    dic_del = {}

    for ins in list(dicID.values())[0]:
        try:
            dic_ins[ins] += 1
        except:
            dic_ins[ins] = 1

    for dele in list(dicID.values())[1]:
        try:
            dic_del[dele] += 1
        except:
            dic_del[dele] = 1

    return dic_ins, dic_del


def mergeDic(dic1, dic2):
    """
    Ajoute le contenu de dic2 à dic1, crée une nouvelle clé
    si elle n'exsitait pas dans dic1
    Utilisé dans l'analyse des champs REF et ALT
    """

    for key, value in dic2.items():
        try:
            dic1[key] = dic1[key] + value
        except:
            dic1[key] = value

    return dic1


def descriptionBarString(dic):
    """
    Transforme un dictionnaire {ID : desc, ID : desc, ... }
    en un String -> 'ID : desc \nID : desc ... }
    Utilisé dans l'analyse des champs REF et ALT
    """

    res = ''
    for ID, desc in dic.items():
        if desc:
            res += ID + ' : ' + desc + '\n'
    res = res[:-1]
    return res


def analyzeInsertAndDelTotal():
    """
    Analyze les insertions et délétions pour tous les variants de
    chaque chromosome, en se basant sur les colonnes REF et ALT
    Affichage de 3 plot
    """

    global VARIANTS, COLONNES, selected_quality

    dic_ins = {}  # Forme {A : 1, C : 14, N : 7, ... }
    dic_del = {}  # Idem
    dic_balises = {}  # Forme {ID_Balise1 : 7, ID_Balise2 : 4, ... }
    dic_descriptions = {}

    for variant in VARIANTS:
        ref = variant['REF'][0]  # Base de référence
        ALT = variant['ALT']  # Liste des bases alternatives

        qual = getQuality(variant)

        if qual >= int(selected_quality):

            for alt in ALT:
                alt_split = alt.split(',')  # Pour les cas : 'T,C' -> ['T', 'C']
                for new_alt in alt_split:
                    if isMadeOfNuc(new_alt) and isMadeOfNuc(
                        ref
                    ):  # Si ce sont des bases nucléiques
                        dic_temp = getInsertionsAndDeletions(ref, new_alt)
                        dic_ins_temp, dic_del_temp = insertAndDel2Dic(dic_temp)
                        dic_ins = mergeDic(
                            dic_ins, dic_ins_temp
                        )  # Mise à jour des dictionnaires
                        dic_del = mergeDic(dic_del, dic_del_temp)

                    else:  # Si c'est une balise
                        for case in COLONNES['ALT']:
                            new_alt = new_alt.replace('>', '')
                            new_alt = new_alt.replace('<', '')
                            if case['ID'] == new_alt:  # Si elle est répertoriée
                                try:
                                    dic_balises[
                                        new_alt
                                    ] += 1  # Ajout au dictionnaire des balises
                                except:
                                    dic_balises[new_alt] = 1

                                description = case['Description']
                                if (
                                    new_alt,
                                    description,
                                ) not in dic_descriptions.items():  # Ajout aux descriptions pour légender
                                    dic_descriptions[new_alt] = description

    # Plot des 3 dictionnaires si jugé nécessaire

    plot1, plot2, plot3 = True, True, True

    if sum(dic_ins.values()) > 0:

        x_ins = dic_ins.keys()
        height_ins = dic_ins.values()
        ylabel_ins = "Nombre d'insertions"
        xlabel_ins = 'Nucléotide'
        title_ins = 'Insertions des nucléotides pour tous les variants'

        color_ins = 'green'
        edgecolor_ins = 'blue'

        plotBar(
            x_ins,
            height_ins,
            xlabel_ins,
            ylabel_ins,
            title_ins,
            color_ins,
            edgecolor_ins,
            False,
        )

    else:
        plot1 = False

    if sum(dic_del.values()) > 0:
        x_del = dic_del.keys()
        height_del = dic_del.values()
        ylabel_del = 'Nombre de délétions'
        xlabel_del = 'Nucléotides'
        title_del = 'Délétions des nucléotides pour tous les variants'
        color_del = 'red'
        edgecolor_del = 'yellow'

        plotBar(
            x_del,
            height_del,
            xlabel_del,
            ylabel_del,
            title_del,
            color_del,
            edgecolor_del,
            False,
        )

    else:
        plot2 = False

    if sum(dic_balises.values()) > 0:
        x_balises = dic_balises.keys()
        height_balises = dic_balises.values()
        ylabel_balises = 'Nombre de délétions et insertions'
        xlabel_balises = 'Nucléotides'
        title_balises = 'Délétions et insertions des nucléotides pour tous les variants'
        color_balises = 'white'
        edgecolor_balises = 'orange'
        text_balises = descriptionBarString(dic_descriptions)

        plotBar(
            x_balises,
            height_balises,
            xlabel_balises,
            ylabel_balises,
            title_balises,
            color_balises,
            edgecolor_balises,
            False,
            text_balises,
        )

    else:
        plot3 = False

    if not (plot1 or plot2 or plot3):  # Si rien n'a été affiché
        openNoInfoWindow()


def analyzeInsertAndDelChrom(chrom):
    """
    Analyze les insertions et délétions pour tous les variants
    d'un chromosome donné, en se basant sur les colonnes REF et ALT
    Affichage de 3 plot
    """

    global VARIANTS, COLONNES, selected_quality

    dic_ins = {}  # Forme {A : 1, C : 14, N : 7, ... }
    dic_del = {}  # Idem
    dic_balises = {}  # Forme {ID_Balise1 : 7, ID_Balise2 : 4, ... }
    dic_descriptions = {}

    for variant in VARIANTS:

        qual = getQuality(variant)

        if qual >= int(selected_quality):

            if variant['CHROM'][0] == chrom:  # Selection pour le bon chromosome
                ref = variant['REF'][0]  # Base de référence
                ALT = variant['ALT']  # Liste des bases alternatives

                for alt in ALT:  # Pour chaque élément de la colonne ALT
                    alt_split = alt.split(',')  # Pour les cas : 'T,C' -> ['T', 'C']
                    for new_alt in alt_split:
                        if isMadeOfNuc(new_alt) and isMadeOfNuc(
                            ref
                        ):  # Si ce sont des bases nucléiques
                            dic_temp = getInsertionsAndDeletions(ref, new_alt)
                            dic_ins_temp, dic_del_temp = insertAndDel2Dic(dic_temp)
                            dic_ins = mergeDic(
                                dic_ins, dic_ins_temp
                            )  # Mise à jour des dictionnaires
                            dic_del = mergeDic(dic_del, dic_del_temp)

                        else:  # Si c'est une balise
                            for case in COLONNES['ALT']:
                                new_alt = new_alt.replace('>', '')
                                new_alt = new_alt.replace('<', '')
                                if case['ID'] == new_alt:  # Si elle est répertoriée
                                    try:
                                        dic_balises[
                                            new_alt
                                        ] += 1  # Ajout au dictionnaire des balises
                                    except:
                                        dic_balises[new_alt] = 1

                                    description = case['Description']
                                    if (
                                        new_alt,
                                        description,
                                    ) not in dic_descriptions.items():  # Ajout aux descriptions pour légender
                                        dic_descriptions[new_alt] = description

    plot1 = plot2 = plot3 = True, True, True

    # Plot des 3 dictionnaires si jugé nécessaire

    if sum(dic_ins.values()) > 0:

        x_ins = dic_ins.keys()
        height_ins = dic_ins.values()
        ylabel_ins = "Nombre d'insertions"
        xlabel_ins = 'Nucléotide'
        title_ins = (
            'Insertions des nucléotides pour tous les variants du chromosome '
            + str(chrom)
        )

        color_ins = 'green'
        edgecolor_ins = 'blue'

        plotBar(
            x_ins,
            height_ins,
            xlabel_ins,
            ylabel_ins,
            title_ins,
            color_ins,
            edgecolor_ins,
            False,
        )

    else:
        plot1 = False

    if sum(dic_del.values()) > 0:
        x_del = dic_del.keys()
        height_del = dic_del.values()
        ylabel_del = 'Nombre de délétions'
        xlabel_del = 'Nucléotides'
        title_del = (
            'Délétions des nucléotides pour tous les variants du chromosome '
            + str(chrom)
        )
        color_del = 'red'
        edgecolor_del = 'yellow'

        plotBar(
            x_del,
            height_del,
            xlabel_del,
            ylabel_del,
            title_del,
            color_del,
            edgecolor_del,
            False,
        )

    else:
        plot2 = False

    if sum(dic_balises.values()) > 0:
        x_balises = dic_balises.keys()
        height_balises = dic_balises.values()
        ylabel_balises = 'Nombre de délétions et insertions'
        xlabel_balises = 'Nucléotides'
        title_balises = (
            'Délétions et insertions des nucléotides pour les variants du chromosome '
            + str(chrom)
        )
        color_balises = 'white'
        edgecolor_balises = 'orange'
        text_balises = descriptionBarString(dic_descriptions)

        plotBar(
            x_balises,
            height_balises,
            xlabel_balises,
            ylabel_balises,
            title_balises,
            color_balises,
            edgecolor_balises,
            False,
            text_balises,
        )

    else:
        plot3 = False

    if not (plot1 or plot2 or plot3):  # Si rien n'a été affiché
        openNoInfoWindow()


def qualityTotal():
    """
    Analyse de la qualité des variants pour tous les chromosomes, triés par chromosome
    Affichage de barplot
    """

    global VARIANTS, selected_quality

    dic = (
        {}
    )  # Dictionnaire {CHROM1 : [QUAL1, QUAL2, ... ], CHROM2 : [QUAL3, ... ], ... }
    maxQual = 0

    for variant in VARIANTS:
        CHROM = variant['CHROM'][0]  # Chromosome du variant
        QUAL = getQuality(variant)

        if QUAL >= int(selected_quality):

            if QUAL > maxQual:  # Mise à jour de la qualité maximale
                maxQual = QUAL

            try:
                dic[CHROM].append(QUAL)
            except:
                dic[CHROM] = [QUAL]

    if maxQual > 0:  # Si tous les variants n'ont pas une qualité de 0

        x = dic.keys()
        height = [sum(qual) / len(qual) for qual in dic.values()]
        ylabel = 'Qualité'
        xlabel = 'Chromosome'
        title = 'Moyenne de la qualité des variants pour chaque chromosome'
        color = (0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0)
        edgecolor = 'blue'

        plotBar(x, height, xlabel, ylabel, title, color, edgecolor, True)

    else:
        openNoInfoWindow()


def qualityChrom(chrom):
    """
    Analyse de la qualité des variants pour un chromosome donné
    Affichage de barplot
    """

    global VARIANTS, selected_quality

    dic = {}  # Forme {1 : qualité, 2 : qualité, ... }
    i = 1
    for variant in VARIANTS:
        qual = getQuality(variant)

        if qual >= int(selected_quality):

            if variant['CHROM'][0] == chrom:
                dic[i] = qual
                i += 1

    if sum(dic.values()) > 0:  # Si tous les variants n'ont pas une qualité de 0

        x = dic.keys()
        height = dic.values()
        ylabel = 'Qualité'
        xlabel = 'Variant'
        title = 'Qualité des variants du chromosome ' + str(chrom)
        color = (0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0)
        edgecolor = 'black'

        plotBar(x, height, xlabel, ylabel, title, color, edgecolor, True)

    else:
        openNoInfoWindow()


def REFbaseTypeTotal():
    """
    Analyse des ccurences des base nucléiques de l'allèle de référence
    Affichage de pieChart
    """

    global VARIANTS, selected_quality

    dic = {}  # Dictionnaire des occurences des gènes

    for variant in VARIANTS:
        qual = getQuality(variant)
        if qual >= int(selected_quality):
            for base_long in variant['REF']:
                for base in base_long:
                    if isMadeOfNuc(base):
                        try:
                            dic[base] += 1
                        except:
                            dic[base] = 1

    if sum(dic.values()) > 0:  # Si il y a des valeurs à afficher

        labels = dic.keys()
        sizes = dic.values()
        title = "Nombre d'occurences de chaques base nucléique\ndu génome de référence ayant subi une variation"

        plotPieChart(sizes, labels, title)

    else:
        openNoInfoWindow()


def REFbaseTypeChrom(chrom):
    """
    Occurence des bases nucléiques de l'allèle de référence pour un chromosome donné
    Affichage de pieChart
    """

    global VARIANTS, selected_quality

    dic = {}  # Dictionnaire des occurences des gènes

    for variant in VARIANTS:
        qual = getQuality(variant)
        if qual >= int(selected_quality):
            if variant['CHROM'][0] == chrom:
                for base_long in variant['REF']:
                    for base in base_long:
                        if isMadeOfNuc(base):
                            try:
                                dic[base] += 1
                            except:
                                dic[base] = 1

    if sum(dic.values()) > 0:  # Si il y a des valeurs à afficher

        labels = dic.keys()
        sizes = [c for c in dic.values()]
        title = (
            "Nombre d'occurences de chaques base nucléique \ndu génome de référence pour le chromosome "
            + chrom + "\nayant subi une variation"
        )
        plotPieChart(sizes, labels, title)

    else:
        openNoInfoWindow()


def analyzeFilterTotal():
    """
    Analyse du pourcentage de variants qui ont passé le filtre, trié par chromosome
    Affichage de pieChart
    """

    global VARIANTS, COLONNES, selected_quality

    dic = {}
    dic_legends = {}  # {Nom du (des) filtre utilisé : description}

    for variant in VARIANTS:

        qual = getQuality(variant)

        if qual >= int(selected_quality):

            FILTER = variant['FILTER']  # Etat du filtre du variant
            name = ''

            for filter_name in FILTER:  # Première boucle pour créer le nom
                name += filter_name + ' et '
            name = name[0:-4]  # Enlève le 'et' à la fin

            if name not in dic_legends.keys():

                for filter_name in FILTER:  # Seconde boucle pour les descriptions

                    try:
                        dic_legends[name] += getFilterDescription(filter_name) + '\n'
                    except:
                        dic_legends[name] = getFilterDescription(filter_name) + '\n'

                dic_legends[name] = dic_legends[name][:-1]  # Enlève le '\n' à la fin

            try:
                dic[name] += 1
            except:
                dic[name] = 1

    if dic:
        labels = dic.keys()
        sizes = dic.values()
        title = 'Filtre des variants'
        legends = dic_legends.values()

        plotPieChart(sizes, labels, title, legends)

    else:
        openNoInfoWindow()


def getFilterDescription(ID):
    """
    Renvoie la description d'un FILTER donné
    Utilisé dans l'analyse des FILTER
    """

    global COLONNES

    for FILTER in COLONNES['FILTER']:
        if FILTER['ID'] == ID:
            return FILTER['Description']


def analyzeFilterChrom(chrom):
    """
    Analyse du filtrage des variants pour un chromosome donné
    """

    global VARIANTS, COLONNES, selected_quality

    dic = {}
    dic_legends = {}  # {Nom du (des) filtre utilisé : description}

    for variant in VARIANTS:
        qual = getQuality(variant)

        if qual >= int(selected_quality):

            if variant['CHROM'][0] == chrom:

                FILTER = variant['FILTER']  # Etat du filtre du variant
                name = ''

                for filter_name in FILTER:  # Première boucle pour créer le nom
                    name += filter_name + ' et '
                name = name[0:-4]  # Enlève le 'et' à la fin

                if name not in dic_legends.keys():

                    for filter_name in FILTER:  # Seconde boucle pour les descriptions

                        try:
                            dic_legends[name] += (
                                getFilterDescription(filter_name) + '\n'
                            )
                        except:
                            dic_legends[name] = getFilterDescription(filter_name) + '\n'

                    dic_legends[name] = dic_legends[name][
                        :-1
                    ]  # Enlève le '\n' à la fin

                try:
                    dic[name] += 1
                except:
                    dic[name] = 1

    if dic:
        labels = dic.keys()
        sizes = dic.values()
        title = 'Filtre pour les variants du chromosome ' + str(chrom)
        legends = dic_legends.values()

        plotPieChart(sizes, labels, title, legends)

    else:
        openNoInfoWindow()


############## AFFICHAGES ET FENÊTRES ##############


def plotPieChart(sizes, labels, title, legends=[]):
    """
    Fonction de plot de pie chart
    """

    global root

    root2 = Toplevel(root)

    fig = plt.figure(figsize=(8, 5))

    canvas2 = FigureCanvasTkAgg(fig, master=root2)
    canvas2.draw()

    canvas2.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1.0)

    patches, texts, junk = plt.pie(
        sizes, labels=labels, autopct='%1.1f%%', startangle=180
    )
    plt.title(title, bbox={'facecolor': '0.8', 'pad': 5})

    canvas2.draw()

    if legends:
        plt.legend(patches, legends)

    canvas2.draw()


def plotBar(x, height, xlabel, ylabel, title, color, edgecolor, quality, text=''):
    """
    Fonction de plot de diagramme
    """

    global root

    root2 = Toplevel(root)
    width = 1.0

    fig, ax = plt.subplots(figsize=(7, 5))

    canvas2 = FigureCanvasTkAgg(fig, master=root2)
    canvas2.draw()

    canvas2.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1.0)
    plt.bar(x, height, width, color=color, edgecolor=edgecolor)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)

    if quality:
        plt.axhline(y=20, color='r', linestyle='--', label='salut')
        # TODO
        ax.legend(["Qualité minimale pour un pourcentage d'erreur <= 1 %"])

    if text:
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(
            0.05,
            0.95,
            text,
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment='top',
            bbox=props,
        )

    canvas2.draw()


def openNoInfoWindow():
    """
    Fonction de plot d'une fenêtre indiquant qu'il n'y a pas
    d'informations à afficher
    """

    title = 'No info'
    text = "Pas d'informations à afficher pour cette analyse"
    messagebox.showinfo(title, text)


def openNoAnalyseSelectedWindow():
    """
    Fonction de plot d'une fenêtre indiquant qu'il n'y a pas
    d'analyse sélectionnée
    """

    title = 'No analyse'
    text = "Pas d'analyse sélectionnée"
    messagebox.showinfo(title, text)


def openNoFileSelectedWindow():
    """
    Fonction de plot d'une fenêtre indiquant qu'il n'y a pas
    de fichier sélectionné
    """

    title = 'No file'
    text = "Pas de fichier sélectionné"
    messagebox.showinfo(title, text)


def openFichierIncorrectWindow():
    """
    Fonction de plot d'une fenêtre indiquant que le fichier sélectionné
    est incorrect
    """

    title = 'Incorrect file'
    text = "Le fichier sélectionné est incorrect, veuillez réessayer"
    messagebox.showinfo(title, text)


# Dictionnaire des mots-clés associés aux fonctions
analyseToFunction = {
    'Qualité': [qualityTotal, qualityChrom],
    'Génome de référence': [REFbaseTypeTotal, REFbaseTypeChrom],
    'Filtre': [analyzeFilterTotal, analyzeFilterChrom],
    'Insertions/Délétions': [analyzeInsertAndDelTotal, analyzeInsertAndDelChrom],
}


def main():

    makeWindow()


main()
mainloop()
plt.show(block=False)