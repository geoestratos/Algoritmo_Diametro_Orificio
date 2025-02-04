import numpy as np
from Lagrange_correlaciones import CNL_factor, NL_factor

def Beggs_Brill (t, p, drg, z):
    tpc = 168 + 325 * drg - 12.5 * np.power(drg, 2)
    ppc = 677 + 15 * drg - 37.5 * np.power(drg, 2)
    ppr = (p + 14.7) / ppc
    tpr = (t + 460) / tpc
    
    D = np.exp(0.715 - 1.128 * tpr + 0.42 * np.power(tpr, 2))
    C = 0.132 - 0.32 * np.log(tpr) / np.log(10)
    B = ppr * (0.62 - 0.23 * tpr) + np.power(ppr, 2) * (0.066 / (tpr - 0.86) - 0.037) + (0.32 * np.power(ppr, 6)) / np.exp(20.723 * (tpr - 1))
    a = 1.39 * np.power(tpr - 0.92, 0.5) - 0.36 * tpr - 0.101
    z = a + (1 - a) * np.exp(-B) + C * np.power(ppr, D)
    return z

def Lee (t, dengas, visgas):
    m = 16.04
    K = ((9.4 + 0.02 * m) * np.power(t + 460, 1.5)) / (209 + 19 * m + (t + 460))
    x = 3.5 + 986 / (t + 460) + 0.01 * m
    y = 2.4 + 0.2 * x
    
    visgas = K * 0.0001 * np.exp(x * np.power(dengas / 62.4, y))
    return visgas

def Beggs_Robinson (t, Rs, api, visoil):
    z = 3.0324 - 0.02023 * api
    y = np.power(10, z)
    x = y * np.power(t , -1.163)
    vod = np.power(10, x) - 1
    a = 10.715 * np.power((Rs + 100), -0.515)
    B = 5.44 * np.power((Rs + 150), -0.338)
    
    visoil = a * np.power(vod, B)
    
    P1 = 5.411
    P2 = 5.809934
    
    visoil = P1 * visoil + P2
    return visoil

def Lasater (p, t, drg, dro, api, Rs):
    tol = 0
    P1 = 1.18344
    P2 = 0
    Mo = (63.506 - api) / 0.0996 if api <= 40 else np.power((1048.33 / api), 1.6736)
    while (True):
        yg = (Rs / 379.3) / ((Rs / 379.3) + 350 * dro / Mo)
        pf = 5.043 * np.power(yg , 3) + 3.10526 * np.power(yg, 2) + 1.36226 * yg + 0.119118
        pb = pf * (t + 460) / drg
        
        a = 0.1
        tol = p - pb
        if tol < 0:
            a = a * -1
        Rs = Rs + a
        if np.abs(tol) < 1: 
            break
        
    Rs = P1 * Rs + P2
    return Rs
        
def Factor_z (t, p, drg, z):
    tol = 0
    continue_loop = True
    tolaux = 10
    cont = 0
    
    a1 = 0.31506
    a2 = -1.0467
    a3 = -0.5783
    a4 = 0.5353
    a5 = -0.6123
    a6 = -0.10489
    a7 = 0.68157
    a8 = 0.68446
    
    tpc = 238 + 210 * drg
    ppc = 740 - 100 * drg
    tpr = (t + 460) / tpc
    ppr = (p + 14.7) / ppc
    
    z = 0.5
    
    while continue_loop:
        dr = 0.27 * ppr / (z * tpr)
        aux = z
        z = 1 + dr * (a1 + a2 / tpr + a3 / (np.power(tpr, 3))) + np.power(dr, 2) * (a4 + a5 / tpr) + a5 * a6 * np.power(dr, 5) / tpr + (a7 * np.power(dr, 2) / np.power(tpr, 3)) * (1 + a8 * np.power(dr, 2)) * np.exp(-a8 * np.power(dr, 2))
        tol = np.abs(z - aux)
        if tol < 0.01:
            return z
        if tol < tolaux:
            tolaux = tol
            aux11 = z
        if cont > 10:
            z = aux11
            continue_loop = False
        cont =+ 1 
    return z

def Friccion(Nre, D, Rug, f):
    C = 0.000001
    f = 0.0056 + 0.5 * np.power(Nre, -0.32)
    
    if Nre < 2000:
        f = 64 / Nre
        return f
    
    while (True):
        a = np.power((1 / f), 0.5) 
        B = -2 * (np.log(Rug / (3.71 * D) + 2.51 / (Nre * np.power(f, 0.5))) / np.log(10))
        aux = a - B
        f = f + C
        if aux < 0:
            break
    return f

def Lagrange_Factor_Numero_Viscosidad(NL, CNL):
    x = []
    y = []
    polinomio = 0
    if NL >= 0.1:
        for i in range(0, 9):
            x.append(NL_factor[i + 18])
            y.append(CNL_factor[i + 18])
    if NL >= 0.1:
        for i in range(0, 9):
            x.append(NL_factor[i + 9])
            y.append(CNL_factor[i + 9])
    if NL >= 0.001:
        for i in range(0, 9):
            x.append(NL_factor[i])
            y.append(CNL_factor[i])
    
    if NL < 0.001:
        polinomio = 0.0019
        return polinomio
    if NL > 1:
        polinomio = 0.0133
        return polinomio
    
    for j in range(0, 10):
        c1 = 1
        for i in range(0, 10):
            cociente = 1 if i == j else (NL - x[i]) / (x[j] - x[i])
            c1 = cociente * c1
        termino = c1 * y[j]
        polinomio = termino + polinomio
    CNL = polinomio          
    return CNL