# NWChem

## Installation 

`$ sudo apt install nwchem`

## Optimization Example
```
start water_opt

echo

geometry units angstroms
  O  0.000000  0.000000  0.000000
  H  0.757000  0.586000  0.000000
  H -0.757000  0.586000  0.000000
end

basis
  * library 6-31G**
end

task scf optimize
```

```
$ nwchem water.nw
```

or write output to a file

```
$ nwchem water.nw > water.out
```
