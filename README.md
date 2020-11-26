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
Dans le cadre d'un projet scolaire, nous avons porté notre attention sur les fichier .VCF, et c'est ainsi qu'est né le VCF Analyzer.
Permettant d'afficher des informations sur un fichier en quelques cliques, ceci est avant-tout un projet expérimental, mais fonctionnel.


----------------------------
C. Prérequis et installation
----------------------------

Ce logiciel est distribué sous la forme d'un code Python. Il est donc nécessaire d'avoir Python 3.8 ou une version supérieure installée.
De plus, les librairies suivantes sont demandées afin d'éxécuter le script : sys, re, os, matplotlib, tkinter, random, PIL.
Aucune installation n'est requise, il suffit de télécharger le script.


--------------
D. Utilisation
--------------

Ce logiciel étant sous la forme d'un script Python, pour l'exécuter, il faut :
	- ouvir le script dans un IDLE Python tel que PyCharm ou Python IDLE
	- compiler et exécuter le scipt

Afin d'utiliser le logiciel, il faut d'abord choisir un fichier .VCF à analyser sur votre ordinateur grâce au bouton "Sélectionner un fichier".
Il est ensuite possible de pratiquer plusieurs analyse sur ce fichier, sur l'ensemble des chromosomes ou bien sur un chromose en particulier présent dans le fichier.
Une "Analyse Random" sert à pratiquer une analyse au hasard sur un chromosome au hasard, ou bien tous les chromosomes. 
