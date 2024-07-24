# Hexein Folding

A toy model for proteins on a 2D hexagonal lattice, where "proteins" consist of only two types of amino acids -- hydrophobic and polar. These proteins, called hexeins, are folded in an aqueous environment using a constraint satisfaction formulation. This is then solved using a _divide and concur_ framework and an iterative fixed-point search method known as RRR (Relax-Reflect-Reflect). The _divide_ part of the framework comprises three constraint sets -- motif, perfect matching, and sequentiality.

[//]: # (<img src="https://github.com/user-attachments/assets/385c309e-0a8b-4155-acbf-a2dc32a121b1" width="600" />)

RRR run that finds the fold for the hexein B29 consists of 29 residues in 690 iterations. We use blue, red, and cyan colors to denote polar, hydrophobic, and solvent residues. The color contrast shows confidence in the prediction at a particular iteration.

All the lattice sites on this grid are occupied by the 2 types of amino acids or solvents. The plot on the right-hand side represents the gap between the two constraint sets (divide and concur). The solution is found when the gap/error reaches zero. 

<img src="https://github.com/user-attachments/assets/f0ab7f66-9c71-428c-ae64-6ace7605eaf2" width="50%" height="50%"/>


## Running the program

To run the program, run the following command after updating the input variables in the makefile (TH.mk).

```
make -f TH.mk run
```

## Output Files

Each run would produce three files <ID>.err, <ID>.sol, <ID>.stats

The error file *<ID>.err* contains the error/gap between the divide constraints (constraint set A) and the concur constraint (constraint set B) at different instants during the run.
The sol file *<ID>.sol* contains the values of the sequence and color variables.
The sol file *<ID>.stats* contains the solution and efficiency statistics.
