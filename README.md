# pysz

This is the python code for the calculation of:
- the thermal SZ anguler power spectrum (used in [Makiya et al. 2018](https://arxiv.org/abs/1804.05008) and [Makiya et al. 2020](https://arxiv.org/abs/1907.07870))
- Compton-Y weighted halo bias b_y (used in Chiang et al. 2020)
- redshift dervative of the Compton-Y, dy/dz (used in Chiang et al. 2020)

The tSZ power spectrum calculation is based on the fortran code [szfast](https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumdks/) devoloped by Eiichiro Komatsu,
but there are several important updates e.g., treatment of neutrino mass in the calculation of mass function, mass-concentration relation, etc.
Please see [Makiya et al. 2018](https://arxiv.org/abs/1804.05008) and [Makiya et al. 2020](https://arxiv.org/abs/1907.07870) for the details of the updates.
**Please cite [Komatsu and Seljak 2002](https://arxiv.org/abs/astro-ph/0205468) and [Makiya et al. 2018](https://arxiv.org/abs/1804.05008) if you use the code in your paper.**



Another tSZ power spectrum code, [class_sz](https://github.com/borisbolliet/class_sz), is also based on the szfast.
I have confirmed that the outputs of pysz and class_sz within the numerical uncertainties.

Any questions and requests are welcome.

# Notes
- You can use [this script](https://github.com/ryumakiya/pysz/blob/master/plot_by.ipynb) to reproduce the Figure 10 of Chiang et al. (2020)

# INSTALL
0. The code requires the python modules of [numpy](https://numpy.org/) and [classy](https://lesgourg.github.io/class_public/class.html).
You should manually install classy.

1. Move to source directory

```
cd pysz/pysz/source
```

2. Compile fortran codes
Please specify the fortran compilar you use in the Makefile.
Then,
```
make clean && make
```
All required fortran modules are at the same directory.


3. then go back to top directory of pysz and install python modules

```
cd PYSZDIRECTORY
python setup.py install
```

4. Check the instlattion
Run the python interpreter and type
```
from pysz import pysz
```
if it doesn't show any errors, installtion is succesfully completed.

# Usage
See pysz_usage.ipynb for the usage of the code.
