# pysz

This is the python code for the calculation of the thermal SZ anguler power spectrum.
The tSZ power spectrum calculation is based on the fortran code [fastsz](https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumdks/) devoloped by Eiichiro Komatsu,
but there are several important updates e.g., treatment of neutrino mass in the calculation of mass function, mass-concentration relation, etc.
Please see [Makiya et al. 2018](https://arxiv.org/abs/1804.05008) and [Makiya et al. 2020](for the details of the updates.)

Another tSZ power spectrum code, [class_sz](https://github.com/borisbolliet/class_sz_public), is also based on the fastsz.
I have confirmed that the outputs of pysz and class_sz within the numerical uncertainties.

**Please cite [Komatsu and Seljak 2002](https://arxiv.org/abs/astro-ph/0205468) and [Makiya et al. 2018](https://arxiv.org/abs/1804.05008) if you use the code in your paper.**

Any questions and requests are welcome.

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
