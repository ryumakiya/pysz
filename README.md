# pysz
This code calculates the thermal SZ power spectrum.
The tSZ Cl calculation is based on the fortran code devoloped by Eiichiro Komatsu:
https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumdks/

Please cite [Komatsu and Seljak 2002:https://arxiv.org/abs/astro-ph/0205468] if you use the code in your paper.

# INSTALL
1. Move to source directory

```
cd pysz/pysz/source
```

2. Compile fortran codes

```
make clean && make
```

if you can not make it, modify Makefile in this directory


3. then go back to top directory of pysz and install python modules

```
cd PYSZDIRECTORY
python setup.py install
```

4. see pysz_usage.ipynb for the usage of the code
