# iso

Python script to calculate theoretical isotopic envelope given intrinsic rates and protection factors.

Parameters to be parsed:
* `--ass` is the file containing assignments (format example in file `moprp.ass`)
* `--pep` is the number of peptide to be considered (integer number)
* `--kint` is the name of the file containing the intrinsic rates of the protein (format example in `all_corr.kint`)
* `--pfact` is the name of the file containing the (ln of) protection factors of the protein (format example `best.pfact`)
* `--times` is the name of the file containing the times at which isotopic envelope has to be evaluated (format example in `moprp.times`)

To run the script with example data given in the repository, digit: 

`python iso.py --ass moprp --pep 5 --kint all_corr --pfact best --times moprp`

