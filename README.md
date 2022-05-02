# oxidation_lite
Python scripts used to generate figures in
```
Deng, J., Du, Z., Ghosh, D, Karki, B. & Lee, K.K.M., 2020. A magma ocean origin of the divergent redox evolutions of rocky planetary bodies and early atmospheres. Nature Communications, 11(1), 2007.
```
#### Dependency
```
pandas
lmfit
numpy
scipy : interp1d, opt
burnman - v0.9.0
uncertainties
matplotlib
````
#### Note
Tested on python 3 

under `burnman/eos/` ， replace the following files with the ones in `burnman_patches/eos`

```bash
burnman/eos/helper.py

burnman/eos/modified_vinet.py

burnman/eos/__init__.py

burnman/eos/mie_grueneisen_debye.py
```



under `burnman/minerals`, replace (or add if not existed already) these files with the ones in `burnman_patches/minerals`

```bash
burnman/minerals/other.py

burnman/minerals/LF_2012.py
```



Note I made modifications on burnman 0.9.0. There may be some official updates on the scripts mentioned above. In order to play safe, you may want to manually add my added lines in the these scripts. You can use `diff` command to find the difference between my modified scripts and official scripts.

