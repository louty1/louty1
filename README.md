
# Project Title
Projet Needleman-Wunsch - L3 BISM 2023-2024


## Authors:  Louty Sylle KA 

## Mon code

import argparse
#Partie 1: Importation des modules et analyse des arguments en ligne de commande
parser = argparse.ArgumentParser()
parser.add_argument("fichier1", help="Chemin vers le premier fichier")
parser.add_argument("fichier2", help="Chemin vers le deuxième fichier")
parser.add_argument("--gap", type=int, default=-10, help="Valeur du gap pour l'alignement")

argument = parser.parse_args()

#Partie 2: Lecture des séquences à partir des fichiers
with open(argument.fichier1) as f1:
    seq1 = f1.read()

with open(argument.fichier2) as f2:
    seq2 = f2.read()

#Partie 3: Affichage des séquences
print(seq1)
print(seq2)

#Partie 4: Définition du dictionnaire de scores de substitution et de la pénalité de gap
alignment_seq = {
    ('A', 'A'): 2, ('A', 'a'): 2, ('a', 'A'): 2, ('a', 'a'): 2,
    ('C', 'C'): 2, ('C', 'c'): 2, ('c', 'C'): 2, ('c', 'c'): 2,
    ('G', 'G'): 2, ('G', 'g'): 2, ('g', 'G'): 2, ('g', 'g'): 2,
    ('T', 'T'): 2, ('T', 't'): 2, ('t', 'T'): 2, ('t', 't'): 2,
    ('U', 'U'): 2, ('U', 'u'): 2, ('u', 'U'): 2, ('u', 'u'): 2,
    ('A', 'G'): -1, ('G', 'A'): -1, ('C', 'T'): -1, ('T', 'C'): -1,
    ('a', 'g'): -1, ('g', 'a'): -1, ('c', 't'): -1, ('t', 'c'): -1,
    ('A', 'C'): -1, ('C', 'A'): -1, ('A', 'c'): -1, ('a', 'C'): -1, ('a', 'c'): -1,
    ('A', 'T'): -1, ('T', 'A'): -1, ('a', 'T'): -1, ('a', 't'): -1, ('T', 'a'): -1, ('t', 'A'): -1,
    ('C', 'G'): -1, ('G', 'C'): -1, ('C', 'g'): -1, ('c', 'G'): -1, ('c', 'g'): -1,
    ('C', 'U'): -1, ('C', 'u'): -1, ('c', 'U'): -1, ('c', 'u'): -1,
    ('U', 'C'): -1, ('U', 'c'): -1, ('u', 'C'): -1, ('u', 'c'): -1,
    ('U', 'T'): -1, ('U', 't'): -1, ('u', 'T'): -1, ('u', 't'): -1,
}

gap_penalty = argument.gap

#Partie 5: Définition de la fonction Needleman-Wunsch
def needleman_wunsch(seq1, seq2, alignment_seq, gap_penalty):
    m = len(seq1)
    n = len(seq2)

    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        score_matrix[i][0] = gap_penalty * i

    for j in range(1, n + 1):
        score_matrix[0][j] = gap_penalty * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = alignment_seq.get((seq1[i - 1], seq2[j - 1]), 0)

            diag_score = score_matrix[i - 1][j - 1] + match_score
            up_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty

            score_matrix[i][j] = max(diag_score, up_score, left_score)

    score = score_matrix[m][n]

    i, j = m, n
    align1, align2 = '', ''
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + alignment_seq.get((seq1[i - 1], seq2[j - 1]), 0):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = '-' + align2
            j -= 1

    return score, align1, align2

#Partie 6: Appel de la fonction Needleman-Wunsch et affichage des résultats
score, align1, align2 = needleman_wunsch(seq1, seq2, alignment_seq, gap_penalty)

print("Score:", score)
print("Alignement 1:", align1)
print("Alignement 2:", align2)

## Résultat

ACATAAAAATACAATGGAGCCATTACAAAAGAATAAATAAGGTCTTCCTG

TAAAAATACAGCAAAAAAATTAAAATATAAGGGATAACTTTCAAAT


Score: -4
Alignement 1:ACATAAAAATACAATGGAGCCATTACAAAAGAATAAATAAGGTCTTCCTG
Alignement 2:---TAAAAATACAGCAAAAAAATTAAAATATAAGGGATAACTT-TCAAAT


##  Explication

Le code présente un programme Python qui utilise l'algorithme de Needleman-Wunsch pour aligner deux séquences de nucléotides ou d'acides aminés. Voici un aperçu des principales étapes :

1. Importation et analyse des arguments en ligne de commande à l'aide du module `argparse`, permettant à l'utilisateur de spécifier les fichiers contenant les séquences et la valeur du gap pour l'alignement.

2. Lecture des séquences à partir des fichiers spécifiés dans les arguments.

3. Affichage des séquences lues à partir des fichiers.

4. Définition d'un dictionnaire de scores de substitution pour chaque paire de nucléotides ou d'acides aminés possibles, ainsi que la spécification de la pénalité de gap.

5. Implémentation de la fonction `needleman_wunsch` qui effectue l'alignement des séquences en utilisant l'algorithme de Needleman-Wunsch, avec les séquences, le dictionnaire de scores de substitution et la pénalité de gap en tant que paramètres.

6. Appel de la fonction `needleman_wunsch` avec les séquences, le dictionnaire de scores de substitution et la pénalité de gap, suivi de l'affichage du score d'alignement ainsi que des deux séquences alignées.

Ce programme permet donc d'aligner efficacement deux séquences biologiques et d'afficher les résultats d'alignemente.


## Exemple

# Fonction de Needleman-Wunsch

def needleman_wunsch(seq1, seq2, alignment_seq, gap_penalty):
    m = len(seq1)
    n = len(seq2)

    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        score_matrix[i][0] = gap_penalty * i

    for j in range(1, n + 1):
        score_matrix[0][j] = gap_penalty * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = alignment_seq.get((seq1[i - 1], seq2[j - 1]), 0)

            diag_score = score_matrix[i - 1][j - 1] + match_score
            up_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty

            score_matrix[i][j] = max(diag_score, up_score, left_score)

    score = score_matrix[m][n]

    i, j = m, n
    align1, align2 = '', ''
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + alignment_seq.get((seq1[i - 1], seq2[j - 1]), 0):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return score, align1, align2

# Exemple d'utilisation
seq1 = "ACGT"
seq2 = "AGT"
alignment_seq = {'A': {'A': 2, 'G': -1, 'C': -1, 'T': -1},
                 'G': {'A': -1, 'G': 2, 'C': -1, 'T': -1},
                 'C': {'A': -1, 'G': -1, 'C': 2, 'T': -1},
                 'T': {'A': -1, 'G': -1, 'C': -1, 'T': 2}}
gap_penalty = -1

# Appel de la fonction needleman_wunsch avec les paramètres corrects
score, align1, align2 = needleman_wunsch(seq1, seq2, alignment_seq, gap_penalty)

# Affichage des résultats
print("Score:", score)
print("Alignement 1:", align1)
print("Alignement 2:", align2)


Score: -1
Alignement 1: ACGT
Alignement 2: -AGT


## Explication de cet exemple

Dans cet exemple, deux séquences "seq1" et "seq2" sont alignées en utilisant une matrice de substitution prédéfinie et une pénalité de gap de -1. Le score d'alignement ainsi que les séquences alignées sont affichés à la fin.



## Remarque

J'ai eu quelques soucis, notamment pour ouvrir les fichiers contenant les séquences et pour créer le dictionnaire de scores de substitution.
