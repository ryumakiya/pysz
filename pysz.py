import numpy as np
from ctypes import * 
import sys, os
import platform
from scipy import interpolate
from classy import Class
import time

# sould match with that defined in FORTRAN codes?
MAXSIZE = 4096

# checking platform
LIBDIR = os.path.dirname(__file__)
# LIBDIR = (os.path.expanduser('~'))+'/pymy/pysz'
sys.path.append(LIBDIR)


class tsz_cl:
    def __init__(self):
        # print 'Class for tSZ Cl'
        # self.ptilde = np.loadtxt(LIBDIR+'/aux_files/ptilde.txt')
        self.fort_lib_cl = cdll.LoadLibrary(LIBDIR+"/source/calc_cl")

        self.fort_lib_cl.calc_cl_.argtypes = [
                                    POINTER(c_double), #h0
                                    POINTER(c_double), #obh2
                                    POINTER(c_double), #och2
                                    POINTER(c_double), #mnu
                                    POINTER(c_double), #bias
                                    POINTER(c_int64), #pk_nk
                                    POINTER(c_int64), #pk_nz
                                    np.ctypeslib.ndpointer(dtype=np.double), #karr
                                    np.ctypeslib.ndpointer(dtype=np.double), #pkarr
                                    POINTER(c_int64), #nl
                                    np.ctypeslib.ndpointer(dtype=np.double), #ell
                                    np.ctypeslib.ndpointer(dtype=np.double), #yy
                                    np.ctypeslib.ndpointer(dtype=np.double), #tll
                                    POINTER(c_int64) #flag_nu
                                    ]
        self.fort_lib_cl.calc_cl_.restype = c_void_p

        # Calcualtion setup
        self.kmin = 1e-4
        self.kmax = 10.
        self.zmax = 4.
        self.nk_pk = 200
        self.nz_pk = 101

    def get_tsz_cl(self,ell_arr,params):
        obh2 = params['obh2']
        och2 = params['och2']
        theta = params['theta']
        As = params['As']
        ns = params['ns']
        mnu = params['mnu']
        mass_bias = params['mass_bias']
        flag_nu = params['flag_nu']

        pars = {'output':'mPk','100*theta_s':theta,
                'omega_b':obh2,'omega_cdm':och2,
                'A_s':As,'n_s':ns,\
                'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                'T_ncdm':0.71611,\
                'P_k_max_h/Mpc': self.kmax,'z_max_pk':self.zmax,\
                'deg_ncdm':3.}
        cosmo = Class()
        cosmo.set(pars)
        cosmo.compute()
        h0 = cosmo.h()

        kh_arr = np.logspace(np.log10(self.kmin),np.log10(self.kmax),self.nk_pk)
        kh = np.zeros((self.nz_pk,self.nk_pk))
        pk = np.zeros((self.nz_pk,self.nk_pk))
        pk_zarr = np.linspace(0.,self.zmax,self.nz_pk)
        for i in range(self.nz_pk):
            kh[i,:] = kh_arr
            if flag_nu == 0:
                pk[i,:] = np.array([cosmo.pk(k*h0,pk_zarr[i])*h0**3 for k in kh_arr])
            elif flag_nu == 1:
                pk[i,:] = np.array([cosmo.pk_cb(k*h0,pk_zarr[i])*h0**3 for k in kh_arr])

        ####  set variables for fortran codes ###
        nk = byref(c_int64(self.nk_pk)) # indx_z, indx_k
        nz = byref(c_int64(self.nz_pk))
        
        # params
        h0_in = byref(c_double(h0))
        obh2_in = byref(c_double(obh2))
        och2_in = byref(c_double(och2))
        mnu_in = byref(c_double(mnu))
        mass_bias_in = byref(c_double(mass_bias))
        flag_nu_in = byref(c_int64(flag_nu))
        
        # outputs
        nl = len(ell_arr)
        cl_yy = np.zeros(nl)
        tll = np.zeros((nl,nl))
        nl = c_int64(nl)
    
        t1 = time.time()
        self.fort_lib_cl.calc_cl_(
                h0_in, obh2_in, och2_in, mnu_in,\
                mass_bias_in, nk, nz,\
                np.array(kh),np.array(pk),\
                nl,np.array(ell_arr),\
                cl_yy,tll,\
                flag_nu_in
                )
        print time.time()-t1
        return cl_yy
