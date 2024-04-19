# About
some quant chem calc test using PSI4.<br>
"proof of concept". Among that, to use a smiles to XYZ format instead of starting from XYZ as is the input format for PSI4.<br>
## heat of formation
Initial idea to calculat energy of a compound to then calculcate _an approximation_ of reaction energy (delta) as a means to check if a reaction is endo- or exothermic.<br>
<br>
Accuracy on this level unfortunately not sufficient, maybe for simple mols like water, where I am close to the reported heat of formation. On the other hand already with e.g.<br>
__nitrobenzene + H2 -> aniline + water__<br>
log units off  .... so one can say that the package works, but some of the approaches (e.g. choosen parameters) are not always correct.
<br>

## torsional energy profile
Another example (in the jupyter notebooks) are simplistic torsional energy calc/visualization, purpose to check the profile of a unhindered compound versus a atropisomeric compound.<br>
The absolute values aren't correct? Also a q of choosen parameters....
<br>

# installation
had some issues with installation, e.g.:<br>
tested, but didn't work:<br>
```bash
conda create -n psi4_env python=3.7
conda activate psi4_env
conda install psi4 psi4
conda install -c conda-forge rdkit
```

a second attempt:<br>
```bash
conda create --name psi4_env python=3.11
pip install psi4-step
pip install rdkit
```

did not work since not supported in 3.11 (for now)<br>
<br>
finally, what worked:<br>
```bash
conda create --name psi4_env python=3.10
conda install psi4 psi4
pip install rdkit
```

(jupyter also required for the notebooks)

### temp files
regularly check you tmp files folder (or add an autodelete function), some of the files can get quite big (giga bytes of data).<br>
in linux it's `/tmp/psi*`

### references
e.g., but not exclusively: https://pubs.acs.org/doi/10.1021/jo1012898<br>
code wise, using a module by Steven Kearnes, license 3-clause BSD, check the conformers.py file for details.

### license
MIT - bascially, if anyone finds this interesting they can do what they want with it (see license.md for details), after all, most of this is based on other open-source stuff.