# ORCA

Download here:
https://www.faccts.de/customer



## Example commands

orca orca_example.inp > orca_example.out


## Visualize the output file with Avogadro

avogadro orca_example.out



## H2O input file example:

```
! RKS B97-3c opt freq
# B97-3c: DFT functional
# opt: geometry optimization
# freq: Frequency analysis

* xyz 0 1
O          0.000000      0.000000      0.000000
H          0.000000      0.759337      0.596044
H          0.000000     -0.759337      0.596044
*
```


## Glycine input file example:

```
! RKS B97-3c opt freq
# B97-3c: DFT functional
# opt: geometry optimization
# freq: Frequency analysis

* xyz 0 1
C
C   1.408137   0.000000   0.000000
H   1.086655   1.042889  -0.000002
H   1.086655  -1.042889   0.000000
H   2.518537   0.000000   0.000001
N
O   0.000000   0.000000  -0.000001
H  -0.850878  -0.000002  -0.000003
O   3.297982   0.549679   0.000003
H   3.391903   1.540471  -0.000003
O   2.386395  -1.162101   0.000002
H   2.101263  -2.115033   0.000004
*
```


## Adenine input file example:

```
! RKS B97-3c opt
# B97-3c: DFT functional
# opt: geometry optimization

* xyz 0 1
C         -0.10173      1.22272      0.00000
C          1.28291      1.14441      0.00000
N          1.89679     -0.08053      0.00000
C          1.03666     -1.11189      0.00000
C         -0.34293     -1.07720      0.00000
N         -0.96316      0.13459      0.00000
N         -0.81747      2.42777      0.00000
N          1.44229     -2.39247      0.00000
C          0.33446     -3.08447      0.00000
N         -0.78310     -2.33647      0.00000
H          2.41724     -2.73351      0.00000
H          0.33742     -4.16747      0.00000
H         -1.80247      2.39377      0.00000
H         -0.31747      3.29777      0.00000
H          1.91291      2.02441      0.00000
*
```