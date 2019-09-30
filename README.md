# anglefinder

A command-line script searching for particular *d*-ligand-*d* (e.g., Cu-O-Cu) angles in a given set of crystal structures (cif files).

### Synopsis:
```
usage: anglefinder.py [-h] -d D -l L --ddmax DDMAX --dlmax DLMAX -a ANGLES
                      ANGLES [--dvalence DVALENCE] [--nligands NLIGANDS]
                      [--nangles NANGLES] [--prefix PREFIX]
                      cifdir

Searches for specific d-ligand-d angles in cif ciles.

positional arguments:
  cifdir                the directory with the cif files.

optional arguments:
  -h, --help            show this help message and exit
  -d D                  d atom type.
  -l L                  ligand atom type.
  --ddmax DDMAX         maximal d-d distance.
  --dlmax DLMAX         maximal d-ligand distance.
  -a ANGLES ANGLES, --angles ANGLES ANGLES
                        minimal and maximal d-ligand-d angles.
  --dvalence DVALENCE   valence of the d atom. Default: None.
  --nligands NLIGANDS   maximal number of ligands. Default: 2.
  --nangles NANGLES     number of angles per line (output). Default: 6.
  --prefix PREFIX       file names for MD and PDF output files.
  ```

### Strategy:
For each cif file, we do the following:

1. Find all sites fully occupied with a given d-element.
2. For every d-site, find its ligands and neighboring d-atoms.
3. For every neighbor, we check if any of the ligands surrounding the central atom is also a ligand of this neighboring atom. If it is the case, we calculate the central-ligand-neighbor angle.
4. If the number of angles is equal to `nligands` or more, we add the respective angles to the preliminary list.
5. If any of the angles in the preliminary list falls into the given range, we put the respective material into the final output list.

The output is given as a markdown and a pdf file (the latter is generated using pandoc with the [eisvogel](https://github.com/Wandmalfarbe/pandoc-latex-template) template. Structures are imported via pymatgen.

### Example:
In the directory `cifs`, find all structures with Cu--O--Cu angles falling in the 92.5° to 94.5° range, with Cu--O distances of 2.1 Å or smaller and with Cu atoms whose separation does not exceed 3.1 Å:

        ./anglefinder.py -d Cu -l O --ddmax 3.1 --dlmax 2.1 -a 92.5 94.5 cifs/

### Requirements
* python 3 with [pymatgen](http://pymatgen.org/) and [pypandoc](https://pypi.org/project/pypandoc/)
* pandoc
