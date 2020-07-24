import numpy as np

#constants
coul = 332.0641

def degree_to_rad(degree):
    return (degree*(np.pi/180))

def calc_Uvdw(Ai,Aj,Bi,Bj,r):
    Uvdw = ((Ai * Aj)/r**12)  - ((Bi*Bj)/r**6)
    return Uvdw

def calc_Ucoul(qi,qj,r):
    Ucoul = coul*((qi*qj)/r)
    return(Ucoul)

def calc_harmonic(kb, b, b0):
    Uharmonic = (0.5*kb)*(b-b0)**2.
    return Uharmonic

def calc_torsion(parameters, radian):
    Utors=0
    for parameter in parameters:
        k = parameter[0]
        n = abs(parameter[1])
        d = degree_to_rad(parameter[2])
        Utors += ((0.5*k) * (1 + (np.cos((n * radian) - d ))))
    return(Utors)