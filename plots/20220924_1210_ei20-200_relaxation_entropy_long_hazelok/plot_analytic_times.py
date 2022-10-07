import numpy as np

mp = 1.67262*10**(-27)
e = 1.602176565*10**(-19) * 3*10**9
qe_SI = 1.602176565*10**(-19)
qe_cgs = 4.8*10**(-10)
eps0 = 8.854*10**(-12)

# boltzmann constant, expressed in erg/eV
kb_cgs = 1.6*10**(-12)

def time_thermalization(t, ma, mb, n, l, Ta, Tb):
    ma = ma * 1000
    mb = mb * 1000
    n = n * 10**(-6)
    nu =  1.8*10**(-19) * np.sqrt(ma*mb) * n * l / (ma*Ta + mb*Tb)**(1.5)
    tau = 1. / (2*nu)
    return tau

def time_isotropization(t, m, n, l, Tx, Ty):
    Tpar = Tx
    Tperp = Ty
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
    ret = [tau/4, tau/4]
    # ret[0] =  (Tpar+Tperp)/2. + (Tpar-Tperp)/2.* np.exp(-4.*t/tau)
    # ret[1] =  (Tperp+Tpar)/2. + (Tperp-Tpar)/2.* np.exp(-2.*t/tau)
    return ret

def Tisothermal(t, ma, mb, n, l, Tax, Tay, Tbx, Tby):
    # compute times
    Ta = np.sqrt(Tax*Tay)
    Tb = np.sqrt(Tbx*Tby)
    Ta = (Tax+Tay)/2
    Tb = (Tbx+Tby)/2
    tau_a_iso = time_isotropization(t, ma, n, 14.55, Tax, Tay)
    tau_b_iso = time_isotropization(t, mb, n, 14.11, Tbx, Tby)
    tau_a_thermal = time_thermalization(t, ma, mb, n, 6.6, Ta, Tb)
    tau_b_thermal = time_thermalization(t, mb, ma, n, 6.6, Tb, Ta)
    
    ret = [0,0,0,0]
    ret[0] = (Ta+Tb)/2. + (Ta-Tb)/2.* np.exp(-t/tau_a_thermal) + (Tax-Tay)/2 * np.exp(-t/tau_a_iso[0])
    ret[1] = (Ta+Tb)/2. + (Ta-Tb)/2.* np.exp(-t/tau_a_thermal) + (Tay-Tax)/2 * np.exp(-t/tau_a_iso[1])
    ret[2] = (Ta+Tb)/2. + (Tb-Ta)/2.* np.exp(-t/tau_b_thermal) + (Tbx-Tby)/2 * np.exp(-t/tau_b_iso[0])
    ret[3] = (Ta+Tb)/2. + (Tb-Ta)/2.* np.exp(-t/tau_b_thermal) + (Tby-Tbx)/2 * np.exp(-t/tau_b_iso[1])

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

def time_relaxation_lingam(n, ma, mb, l, Ta):
    mu = mb / mp
    print("nu {} n {:.2E} ma {:.2E} l {:.2E} Ta {}".format(np.sqrt(mu), n, ma, l, Ta))
    nu = n * 4 * np.sqrt(2*np.pi*ma) * qe_SI**4 * l / (3*(4*np.pi*eps0)**2 * (Ta*qe_SI)**(1.5)) / ma
    # nu = nu * (1 + ma/mb)
    # nu = nu * np.sqrt(mu)
    tau = 1 / nu
    # tau = tau / np.sqrt(mu)

    return tau

def time_relaxation_hazel2(n, ma, mb, l, Ta):
    ncm = n * 10**-6
    mg = ma * 1000
    vt = 4.19*10**7 * np.sqrt(Ta)
    print("n {:.2E} ma {:.2E} l {:.2E} Ta {} vt {:.2E}".format(n, ma, l, Ta, vt))

    tau = 3 * np.sqrt(np.pi) * vt**3 * mg**2 / (4 * ncm * 4 * np.pi *qe_cgs**4 * l)

    print("tau {:.2E}".format(tau))

    return tau

def time_relaxation_hazel(n, ma, mb, l, Ta):
    ncm = n * 10**-6
    mg = ma * 1000
    vt = 4.19*10**7 * np.sqrt(Ta)
    print("n {:.2E} ma {:.2E} l {:.2E} Ta {} vt {:.2E}".format(n, ma, l, Ta, vt))
    vt = np.sqrt(qe_SI * Ta / ma) * 100
    vt = vt * 1.5
    print("n {:.2E} ma {:.2E} l {:.2E} Ta {} vt {:.2E}".format(n, ma, l, Ta, vt))

    nu = ncm * 4 * np.pi * qe_cgs**4 * l / (mg**2 * vt**3)

    print("tau {:.2E}".format(1/nu))

    return 1 / nu
    


def PslowDown(t, va, ma, mb, Ta, Tb, l, n):
    tau_relax = 4/np.sqrt(2)*time_relaxation_lingam(n, ma, mb, l, Ta)
    # tau_relax = time_relaxation_hazel2(n, ma, mb, l, Ta)
    # tau_relax = time_relaxation_hazel(n, ma, mb, l, Ta)
    print("tau relaxation: {:.2E}".format(tau_relax))
    ret = va * np.exp(-t/tau_relax)
    return ret