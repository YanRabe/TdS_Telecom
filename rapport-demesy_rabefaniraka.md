---
title: "Rapport - Codage Canal"
author:
- Amaury Demesy
- Yan Rabefaniraka
date: "$1^{er}$ Juin 2025"
output: pdf_document
geometry: margin=2cm
---

**Remarque :** Fichier engendré avec:

~~~
pandoc -o rapport-demesy_rabefaniraka.pdf rapport-demesy_rabefaniraka.md
~~~

# Introduction

Le codage canal est l'ajout selon une loi donnée de redondance d'information au sein de notre transmission. Cela sert à minimser les erreurs de transmission en vérifiant à la réception que la loi de codage utilisée par l'émetteur a bien été respectée.

L'objectif de ce projet est d'étudier l'impact de l'ajout du codage canal sur l'efficacité spectrale et l'efficacité en puissance d'une chaîne de transmission.

Pour ce faire, seront implémentées le **codage en bloc linéaires** avec le *code de Hamming*, puis le **codage convolutif** avec *l'algorithme de Viterbi*.

# Pré-requis (Transmission passe-bas équivalente BPSK)

Les implémentations seront limitées à une modulation binaire (éléments se trouvant dans ${0,1}$). Nous utiliserons plus spécifiquement un mapping ***BPSK*** (ici c'est équivalent à un mapping 2-ASK avec deux valeurs opposées sur l'axe des réels), pour un débit binaire de $3000bits/sec$, avec une fréquence d'échantillonnage de 12kHz. Nous choisirons un filtre de mise en forme et de réception rectangulaire de durée Ts.

## Validation de la chaîne de transmission

Avec $n_0 = 4$, nous remarquons que notre chaîne de transmission sans bruit a bien un TEB nul.

Là où, avec bruit, il se rapproche bien du TEB théorique pour un mapping binaire BPSK.

# Implémentation en codes à bloc linéaires

blabla Hamming

## Cas du décodage dur et son impact

## Cas du décodage souple et son impact

