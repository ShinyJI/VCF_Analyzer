-------------------
VCF Analyzer README
       2020
-------------------



--------
Sommaire
--------

A. Licence
B. Pourquoi un analyseur de fichiers VCF
C. Prérequis et installation
D. Utilisation



----------
A. Licence
----------

Ce logiciel est la création de Jacques Imbert et Cyril Drago. Il a été créé dans le cadre d'un projet scolaire.

Le VCF Analyzer est open-source et peut être copié et/ou modifié, ainsi qu'être intégré à d'autres applications dans son intégralité.
La distribution de ce produit est autorisée.
Ce produit est interdit à la vente.

La licence ci-dessus doit être incluse dans toute reproduction complète ou partielle de ce logiciel.
De plus, il est aussi demandé de créditer les deux auteurs.


----------------------------------------
B. Pourquoi un analyseur de fichiers VCF
----------------------------------------

Un analyseur de fichier est un outil utile qui permet d'extraire de manière fiable et rapide des informations.
Dans le cadre d'un projet scolaire, nous avons porté notre attention sur les fichiers .VCF, et c'est ainsi qu'est né le VCF Analyzer.
Permettant d'afficher des informations sur un fichier en quelques clics, ceci est avant tout un projet expérimental, mais fonctionnel.


----------------------------
C. Prérequis et installation
----------------------------

Ce logiciel est distribué sous la forme d'un code Python. Il est donc nécessaire d'avoir Python 3.8 ou une version supérieure installée.
De plus, les librairies suivantes sont demandées afin d'exécuter le script : sys, re, os, matplotlib, tkinter, random, PIL.
Aucune installation n'est requise, il suffit de télécharger le script et de l'exécuter.
Il faut par ailleurs avoir, dans le même répertoire que le script, les fichiers image_adn.png et icone_vcf.ico.


--------------
D. Utilisation
--------------

Ce logiciel étant sous la forme d'un script Python, pour l'exécuter, il faut :
	- ouvrir le script dans un IDLE Python tel que PyCharm ou Python IDLE
	- compiler et exécuter le script

Afin d'utiliser le logiciel, il faut d'abord choisir un fichier .VCF à analyser sur votre ordinateur grâce au bouton "Sélectionner un fichier".
Il est ensuite possible de pratiquer plusieurs analyses sur ce fichier, sur l'ensemble des chromosomes ou bien sur un chromosome en particulier présent dans le fichier.
Une "Analyse Random" sert à pratiquer une analyse au hasard sur un chromosome au hasard, ou bien tous les chromosomes.
Il est possible de sélectionner la qualité minimale des variants sur lesquels effectuer les analyses. Si 0 est choisi, tous les variants seront pris en compte (même ceux dont la qualité n'est pas renseignée).
