
# Project Title
Projet Needleman-Wunsch - L3 BISM 2023-2024


## Authors:  Louty Sylle KA 

## Mon code

# Partie 1: Importation des modules et analyse des arguments en ligne de commande

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("fichier1", help="Chemin vers le premier fichier")

parser.add_argument("fichier2", help="Chemin vers le deuxième fichier")

parser.add_argument("--gap", type=int, default=-10, help="Valeur du gap pour l'alignement")

argument = parser.parse_args()

# Partie 2: Lecture des séquences à partir des fichiers

with open(argument.fichier1) as f1:
    seq1 = f1.read()

with open(argument.fichier2) as f2:
    seq2 = f2.read()

# Partie 3: Affichage des séquences

print(seq1)

print(seq2)

# Partie 4: Définition du dictionnaire de scores de substitution et de la pénalité de gap

alignment_seq =
{
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

# Partie 5: Définition de la fonction Needleman-Wunsch


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

Ce programme permet donc d'aligner efficacement deux séquences biologiques et d'afficher les résultats d'alignement.


## Exemple


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
    
seq1 = "ACGT"

seq2 = "AGT"

alignment_seq = {'A': {'A': 2, 'G': -1, 'C': -1, 'T': -1},
                 'G': {'A': -1, 'G': 2, 'C': -1, 'T': -1},
                 'C': {'A': -1, 'G': -1, 'C': 2, 'T': -1},
                 'T': {'A': -1, 'G': -1, 'C': -1, 'T': 2}}

gap_penalty = -1

score, align1, align2 = needleman_wunsch(seq1, seq2, alignment_seq, gap_penalty)

print("Score:", score)

print("Alignement 1:", align1)

print("Alignement 2:", align2)


## Résultat

Score: -1

Alignement 1: ACGT

Alignement 2: -AGT


## Explication de cet exemple

Dans cet exemple, deux séquences "seq1" et "seq2" sont alignées en utilisant une matrice de substitution prédéfinie et une pénalité de gap de -1. Le score d'alignement ainsi que les séquences alignées sont affichés à la fin.


## Matrice de Blossum62

blosum62 = {
 ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
 ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
 ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
 ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
 ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
 ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
 ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
 ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
 ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
 ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
 ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
 ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
 ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
 ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
 ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
 ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
 ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
 ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
 ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
 ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
 ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
 ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
 ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
 ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
 ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
 ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
 ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
 ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
 ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
 ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
 ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
 ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
 ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
 ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
 ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
 ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
 ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
 ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
 ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
 ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
 ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
 ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
 ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
 ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
 ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
 ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
 ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
 ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
 ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
 ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
 ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
 ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
 ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
 ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
 ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
 ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
 ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
 ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
 ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
 ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
 ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
 ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
 ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
 ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
 ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
 ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
 ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
 ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
 ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4


## Remarque

J'ai eu quelques soucis, notamment pour ouvrir les fichiers contenant les séquences et pour créer le dictionnaire de scores de substitution.
