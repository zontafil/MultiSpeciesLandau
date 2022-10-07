import numpy as np

mp = 1.67262*10**(-27)
e = 1.602176565*10**(-19) * 3*10**9
qe_SI = 1.602176565*10**(-19)
qe_cgs = 4.8*10**(-10)
eps0 = 8.854*10**(-12)

# boltzmann constant, expressed in erg/eV
kb_cgs = 1.6*10**(-12)

def time_thermalization(t, ma, mb, n, l, Ta, Tb):
    # taken from NRL formulary
    ma = ma * 1000
    mb = mb * 1000
    n = n * 10**(-6)
    nu =  1.8*10**(-19) * np.sqrt(ma*mb) * n * l / (ma*Ta + mb*Tb)**(1.5)
    tau = 1. / (2*nu)
    print("tau_thermal {:.2E}".format(tau))
    return tau

def time_isotropization(t, m, n, l, Tx, Ty):
    # taken from NRL formulary
    Tpar = Tx
    Tperp = Ty
    # convert to CGS units
    m = m * 1000
    n = n * 10**(-6)
    A = Tperp/Tpar - 1
    if (A>0):
        f = np.arctan(np.sqrt(A))/np.sqrt(A)
    else:
        f = np.arctanh(np.sqrt(-A))/np.sqrt(-A)
    nu = (2*np.sqrt(np.pi)*e**4*n*l)/(np.sqrt(m)*(kb_cgs*Tpar)**(1.5)) \
        * A**(-2) * (-3 + (A+3)*(f))
    tau = 1 / nu
    ret = tau/4
    print("tau_isotrop {:.2E}".format(tau))
    return ret

def Tisothermal(t, ma, mb, n, Tax, Tay, Tbx, Tby, l00, l01, l11):
    # compute times
    Ta0 = (Tax+Tay)/2
    Tb0 = (Tbx+Tby)/2
    tau_a_iso = time_isotropization(t, ma, n, l00, Tax, Tay)
    tau_b_iso = time_isotropization(t, mb, n, l11, Tbx, Tby)
    tau_a_thermal = time_thermalization(t, ma, mb, n, l01, Ta0, Tb0)
    tau_b_thermal = time_thermalization(t, mb, ma, n, l01, Tb0, Ta0)
    
    ret = [0,0,0,0]
    ret[0] = (Ta0+Tb0)/2. + (Ta0-Tb0)/2.* np.exp(-t/tau_a_thermal) + (Tax-Tay)/2 * np.exp(-t/tau_a_iso)
    ret[1] = (Ta0+Tb0)/2. + (Ta0-Tb0)/2.* np.exp(-t/tau_a_thermal) + (Tay-Tax)/2 * np.exp(-t/tau_a_iso)
    ret[2] = (Ta0+Tb0)/2. + (Tb0-Ta0)/2.* np.exp(-t/tau_b_thermal) + (Tbx-Tby)/2 * np.exp(-t/tau_b_iso)
    ret[3] = (Ta0+Tb0)/2. + (Tb0-Ta0)/2.* np.exp(-t/tau_b_thermal) + (Tby-Tbx)/2 * np.exp(-t/tau_b_iso)

    return ret

def Tthermal(t, ma, mb, n, l, Tax, Tay, Tbx, Tby):
    # compute times
    Ta = np.sqrt(Tax*Tay)
    Ta = (Tax+Tay)/2
    Tb = np.sqrt(Tbx*Tby)
    Tb = (Tbx+Tby)/2
    tau_a_thermal = time_thermalization(t, ma, mb, n, l, Ta, Tb)
    tau_b_thermal = time_thermalization(t, mb, ma, n, l, Tb, Ta)
    
    ret = [0,0]
    ret[0] = (Ta+Tb)/2. + (Ta-Tb)/2.* np.exp(-t/tau_a_thermal)
    ret[1] = (Ta+Tb)/2. + (Tb-Ta)/2.* np.exp(-t/tau_b_thermal)

    return ret

def time_relaxation(n, ma, mb, l, Ta):
    # taken from i.e. Hinton 1976
    nu = n * 16 * np.sqrt(np.pi*ma) * qe_SI**4 * l / (3*(4*np.pi*eps0)**2 * (Ta*qe_SI)**(1.5)) / ma
    tau = 1 / nu
    return tau

def PslowDown(t, va, ma, mb, Ta, Tb, l, n):
    tau_relax = time_relaxation(n, ma, mb, l, Ta)
    print("tau relaxation: {:.2E}".format(tau_relax))
    ret = va * np.exp(-t/tau_relax)
    return ret