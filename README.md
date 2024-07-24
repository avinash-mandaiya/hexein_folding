# hexein_folding

A toy model for proteins on a 2D hexagonal lattice, where "proteins" consist of only two types of amino acids -- hydrophobic and polar.
These proteins, called hexeins, are folded in an aqueous environment using a constraint satisfaction formulation, which is then solved 
using a _divide and concur_ framework and an iterative fixed-point search method known as RRR. The _divide_ part of the framework comprises
three constraint sets -- motif, perfect matching, and sequentiality.

[//]: # (<img src="https://github.com/user-attachments/assets/385c309e-0a8b-4155-acbf-a2dc32a121b1" width="600" />)

RRR run that finds the fold for the hexein B29 consists of 29 residues in 690 iterations. We use blue, red,
and cyan colors to denote polar, hydrophobic, and solvent residues respectively. All the lattice
sites on this grid are occupied by the 2 types of amino acids or solvents. The plot on the right-hand side represents
the gap between the two constraint sets (divide and concur). The solution is found when the gap/error reaches zero. 

<img src="https://github.com/user-attachments/assets/f0ab7f66-9c71-428c-ae64-6ace7605eaf2" width="50%" height="50%"/>

