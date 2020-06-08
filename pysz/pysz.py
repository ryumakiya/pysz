import numpy as np
from ctypes import * 
import sys, os
from classy import Class

# sould match with that defined in FORTRAN codes?
MAXSIZE = 4096

# checking platform
LIBDIR = os.path.dirname(__file__)
# LIBDIR = (os.path.expanduser('~'))+'/Dropbox/py_mymodules/pysz'
sys.path.append(LIBDIR)


class tsz_cl:
    def __init__(self):
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
                                    np.ctypeslib.ndpointer(dtype=np.double), #pk_z_arr
                                    POINTER(c_int64), #n_ell_ptilde
                                    np.ctypeslib.ndpointer(dtype=np.double), #ell_ptilde
                                    np.ctypeslib.ndpointer(dtype=np.double), #ptilde_arr
                                    POINTER(c_int64), #nl
                                    np.ctypeslib.ndpointer(dtype=np.double), #ell
                                    np.ctypeslib.ndpointer(dtype=np.double), #yy
                                    np.ctypeslib.ndpointer(dtype=np.double), #tll
                                    POINTER(c_int64), #flag_nu
                                    POINTER(c_int64), #flag_tll
                                    POINTER(c_double), #zmin
                                    POINTER(c_double), #zmax
                                    POINTER(c_double), #Mmin 
                                    POINTER(c_double) #Mmax
                                    ]
        self.fort_lib_cl.calc_cl_.restype = c_void_p

        self.fort_lib_by = cdll.LoadLibrary(LIBDIR+"/source/calc_by")
        self.fort_lib_by.calc_by_.argtypes = [
                                    POINTER(c_double), #h0
                                    POINTER(c_double), #obh2
                                    POINTER(c_double), #och2
                                    POINTER(c_double), #mnu
                                    POINTER(c_double), #bias
                                    POINTER(c_int64), #pk_nk
                                    np.ctypeslib.ndpointer(dtype=np.double), #karr
                                    np.ctypeslib.ndpointer(dtype=np.double), #pkarr
                                    POINTER(c_int64), #n_ell_ptilde
                                    np.ctypeslib.ndpointer(dtype=np.double), #ell_ptilde
                                    np.ctypeslib.ndpointer(dtype=np.double), #ptilde_arr
                                    POINTER(c_int64), #nz
                                    np.ctypeslib.ndpointer(dtype=np.double), #z_arr
                                    np.ctypeslib.ndpointer(dtype=np.double), #by_arr
                                    np.ctypeslib.ndpointer(dtype=np.double), #dydz_arr
                                    POINTER(c_int64), #flag_nu
                                    POINTER(c_double), #Mmin 
                                    POINTER(c_double) #Mmax
                                    ]
        self.fort_lib_by.calc_by_.restype = c_void_p

        self.fort_lib_dydzdMh = cdll.LoadLibrary(LIBDIR+"/source/calc_dydzdMh")
        self.fort_lib_dydzdMh.calc_dydzdmh_.argtypes = [
                                    POINTER(c_double), #h0
                                    POINTER(c_double), #obh2
                                    POINTER(c_double), #och2
                                    POINTER(c_double), #mnu
                                    POINTER(c_double), #bias
                                    POINTER(c_int64), #pk_nk
                                    POINTER(c_int64), #pk_nz
                                    np.ctypeslib.ndpointer(dtype=np.double), #karr
                                    np.ctypeslib.ndpointer(dtype=np.double), #pkarr
                                    POINTER(c_int64), #n_ell_ptilde
                                    np.ctypeslib.ndpointer(dtype=np.double), #ell_ptilde
                                    np.ctypeslib.ndpointer(dtype=np.double), #ptilde_arr
                                    POINTER(c_double), # z
                                    POINTER(c_double), # Mh
                                    POINTER(c_double), # dydzdMh
                                    POINTER(c_int64), #flag_nu
                                    ]
        self.fort_lib_dydzdMh.calc_dydzdmh_.restype = c_void_p


        # Calcualtion setup
        self.kmin = 1e-3
        self.kmax = 5.
        self.zmin = 1e-5
        self.zmax = 4. 
        self.dz = 0.04
        self.nk_pk = 200
        self.nz_pk = int((self.zmax-self.zmin)/self.dz)+1
        self.Mmin = 1e11
        self.Mmax = 5e15
        self.ptilde_in = np.loadtxt(LIBDIR+'/data/ptilde.txt')
        self.ell_ptilde = np.array(self.ptilde_in[:,0])
        self.ptilde_arr = np.array(self.ptilde_in[:,1])
        self.n_ptilde = len(self.ell_ptilde)

        # Class
        self.cosmo = Class()

    def get_tsz_cl(self,ell_arr,params,zmin=None,zmax=None,Mmin=None,Mmax=None):
        # integration range, optional
        if zmin == None: zmin = self.zmin
        if zmax == None: zmax = self.zmax
        nz_pk = int((zmax-zmin)/self.dz)+1
        if Mmin == None: Mmin = self.Mmin
        if Mmax == None: Mmax = self.Mmax

        obh2 = params['obh2']
        och2 = params['och2']
        As = params['As']
        ns = params['ns']
        mnu = params['mnu']
        mass_bias = params['mass_bias']
        flag_nu_logic = params['flag_nu']
        flag_tll_logic = params['flag_tll']
        if type(flag_nu_logic) != bool:
            print('flag_nu must be boolean.')
            sys.exit()
        if type(flag_tll_logic) != bool:
            print('flag_tll must be boolean.')
            sys.exit()
        
        if flag_nu_logic:
            flag_nu = 1
        else:
            flag_nu = 0
        if flag_tll_logic:
            flag_tll = 1
        else:
            flag_tll = 0
        
        if 'theta' in params.keys():
            theta = params['theta']
            pars = {'output':'mPk','100*theta_s':theta,
                    'omega_b':obh2,'omega_cdm':och2,
                    'A_s':As,'n_s':ns,\
                    'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                    'T_ncdm':0.71611,\
                    'P_k_max_h/Mpc': self.kmax,'z_max_pk':self.zmax,\
                    'deg_ncdm':3.}
            self.cosmo.set(pars)
            self.cosmo.compute()
            h0 = self.cosmo.h()
        elif 'h0' in params.keys():
            h0 = params['h0']
            pars = {'output':'mPk','h':h0,
                    'omega_b':obh2,'omega_cdm':och2,
                    'A_s':As,'n_s':ns,\
                    'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                    'T_ncdm':0.71611,\
                    'P_k_max_h/Mpc': self.kmax,'z_max_pk':self.zmax,\
                    'deg_ncdm':3.}
            self.cosmo.set(pars)
            self.cosmo.compute()

        kh_arr = np.logspace(np.log10(self.kmin),np.log10(self.kmax),self.nk_pk)
        kh = np.zeros((nz_pk,self.nk_pk))
        pk = np.zeros((nz_pk,self.nk_pk))
        pk_zarr = np.linspace(zmin,zmax,nz_pk+1)
        for i in range(self.nz_pk):
            kh[i,:] = kh_arr
            if flag_nu == 0:
                pk[i,:] = np.array([self.cosmo.pk(k*h0,pk_zarr[i])*h0**3 for k in kh_arr])
            elif flag_nu == 1:
                pk[i,:] = np.array([self.cosmo.pk_cb(k*h0,pk_zarr[i])*h0**3 for k in kh_arr])

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
        flag_tll_in = byref(c_int64(flag_tll))
        zmin_in = byref(c_double(zmin))
        zmax_in = byref(c_double(zmax))
        Mmin_in = byref(c_double(Mmin))
        Mmax_in = byref(c_double(Mmax))
        
        # outputs
        nl = len(ell_arr)
        cl_yy = np.zeros((2,nl))
        tll = np.zeros((nl,nl))
        nl = c_int64(nl)
   
        self.fort_lib_cl.calc_cl_(
                h0_in, obh2_in, och2_in, mnu_in,\
                mass_bias_in, nk, nz,\
                np.array(kh),np.array(pk),np.array(pk_zarr),\
                byref(c_int64(self.n_ptilde)), self.ell_ptilde, self.ptilde_arr,\
                nl,np.array(ell_arr),\
                cl_yy,tll,\
                flag_nu_in,flag_tll_in,\
                zmin_in,zmax_in,\
                Mmin_in,Mmax_in
                )

        self.cosmo.struct_cleanup()
        return cl_yy, tll

    def get_by_dydz(self,z_arr,params,Mmin=None,Mmax=None):
        if Mmin == None: Mmin = self.Mmin
        if Mmax == None: Mmax = self.Mmax
        nz = len(z_arr)
        zmax = np.max(z_arr)
    
        obh2 = params['obh2']
        och2 = params['och2']
        As = params['As']
        ns = params['ns']
        mnu = params['mnu']
        mass_bias = params['mass_bias']
        flag_nu_logic = params['flag_nu']
        if type(flag_nu_logic) != bool:
            print('flag_nu must be boolean.')
            sys.exit()
        
        if flag_nu_logic:
            flag_nu = 1
        else:
            flag_nu = 0
        
        if 'theta' in params.keys():
            theta = params['theta']
            pars = {'output':'mPk','100*theta_s':theta,
                    'omega_b':obh2,'omega_cdm':och2,
                    'A_s':As,'n_s':ns,\
                    'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                    'T_ncdm':0.71611,\
                    'P_k_max_h/Mpc': self.kmax,'z_max_pk':zmax,\
                    'deg_ncdm':3.}
            self.cosmo.set(pars)
            self.cosmo.compute()
            h0 = self.cosmo.h()
        elif 'h0' in params.keys():
            h0 = params['h0']
            pars = {'output':'mPk','h':h0,
                    'omega_b':obh2,'omega_cdm':och2,
                    'A_s':As,'n_s':ns,\
                    'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                    'T_ncdm':0.71611,\
                    'P_k_max_h/Mpc': self.kmax,'z_max_pk':zmax,\
                    'deg_ncdm':3.}
            self.cosmo.set(pars)
            self.cosmo.compute()

        kh_arr = np.logspace(np.log10(self.kmin),np.log10(self.kmax),self.nk_pk)
        kh = np.zeros((nz,self.nk_pk))
        pk = np.zeros((nz,self.nk_pk))
        for i in range(nz):
            kh[i,:] = kh_arr
            if flag_nu == 0:
                pk[i,:] = np.array([self.cosmo.pk(k*h0,z_arr[i])*h0**3 for k in kh_arr])
            elif flag_nu == 1:
                pk[i,:] = np.array([self.cosmo.pk_cb(k*h0,z_arr[i])*h0**3 for k in kh_arr])

        ####  set variables for fortran codes ###
        nk = byref(c_int64(self.nk_pk)) # indx_z, indx_k
        
        # params
        h0_in = byref(c_double(h0))
        obh2_in = byref(c_double(obh2))
        och2_in = byref(c_double(och2))
        mnu_in = byref(c_double(mnu))
        mass_bias_in = byref(c_double(mass_bias))
        flag_nu_in = byref(c_int64(flag_nu))
        
        # outputs
        by_arr = np.zeros(nz)
        dydz_arr = np.zeros(nz)
        nz = c_int64(nz)
        Mmin_in = byref(c_double(Mmin))
        Mmax_in = byref(c_double(Mmax))
 
        self.fort_lib_by.calc_by_(
                h0_in, obh2_in, och2_in, mnu_in, \
                mass_bias_in, nk, \
                np.array(kh),np.array(pk), \
                byref(c_int64(self.n_ptilde)), self.ell_ptilde, self.ptilde_arr,\
                nz,np.array(z_arr), \
                by_arr,dydz_arr, \
                flag_nu_in, \
                Mmin_in, Mmax_in
                )

        self.cosmo.struct_cleanup()
        return by_arr, dydz_arr

    def get_dydzdlogMh(self,z,Mh,params):
        obh2 = params['obh2']
        och2 = params['och2']
        As = params['As']
        ns = params['ns']
        mnu = params['mnu']
        mass_bias = params['mass_bias']
        flag_nu_logic = params['flag_nu']
        if type(flag_nu_logic) != bool:
            print('flag_nu must be boolean.')
            sys.exit()
        
        if flag_nu_logic:
            flag_nu = 1
        else:
            flag_nu = 0
        
        if 'theta' in params.keys():
            theta = params['theta']
            pars = {'output':'mPk','100*theta_s':theta,
                    'omega_b':obh2,'omega_cdm':och2,
                    'A_s':As,'n_s':ns,\
                    'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                    'T_ncdm':0.71611,\
                    'P_k_max_h/Mpc': self.kmax,'z_max_pk':z,\
                    'deg_ncdm':3.}
            self.cosmo.set(pars)
            self.cosmo.compute()
            h0 = self.cosmo.h()
        elif 'h0' in params.keys():
            h0 = params['h0']
            pars = {'output':'mPk','h':h0,
                    'omega_b':obh2,'omega_cdm':och2,
                    'A_s':As,'n_s':ns,\
                    'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                    'T_ncdm':0.71611,\
                    'P_k_max_h/Mpc': self.kmax,'z_max_pk':z,\
                    'deg_ncdm':3.}
            self.cosmo.set(pars)
            self.cosmo.compute()

        kh_arr = np.logspace(np.log10(self.kmin),np.log10(self.kmax),self.nk_pk)
        kh = np.zeros((1,self.nk_pk))
        pk = np.zeros((1,self.nk_pk))
        kh[0,:] = kh_arr
        if flag_nu == 0:
            pk[0,:] = np.array([self.cosmo.pk(k*h0,z)*h0**3 for k in kh_arr])
        elif flag_nu == 1:
            pk[0,:] = np.array([self.cosmo.pk_cb(k*h0,z)*h0**3 for k in kh_arr])

        ####  set variables for fortran codes ###
        nk = byref(c_int64(self.nk_pk)) # indx_z, indx_k
        nz = byref(c_int64(1))
         
        # params
        h0_in = byref(c_double(h0))
        obh2_in = byref(c_double(obh2))
        och2_in = byref(c_double(och2))
        mnu_in = byref(c_double(mnu))
        mass_bias_in = byref(c_double(mass_bias))
        flag_nu_in = byref(c_int64(flag_nu))
        
        # outputs
        z_in = byref(c_double(z))
        Mh_in = byref(c_double(Mh))
        dydzdMh = c_double(0.0)
   
        self.fort_lib_dydzdMh.calc_dydzdmh_(
                h0_in, obh2_in, och2_in, mnu_in, \
                mass_bias_in, nk, nz, \
                np.array(kh),np.array(pk), \
                byref(c_int64(self.n_ptilde)), self.ell_ptilde, self.ptilde_arr,\
                z_in,Mh_in, \
                dydzdMh, \
                flag_nu_in \
                )

        self.cosmo.struct_cleanup()
        return dydzdMh.value

