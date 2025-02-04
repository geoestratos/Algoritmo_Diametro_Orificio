import pandas as pd
import numpy as np
from FormulasVarias import Lasater, Beggs_Robinson, Factor_z, Lee, Friccion, Lagrange_Factor_Numero_Viscosidad
from scipy.optimize import root_scalar

def ecuacion(d_orificio, D_tuberia, p1, p2, CD, po, objetivo_qo):
    m = np.power((d_orificio / D_tuberia), 2)
    qo = (CD * (p2/p1) * ((1 - np.power(m, 2)) / (1 - np.power(((p2/p1)*m), 2)))) * (np.pi * np.power(d_orificio, 2) / 4) * np.power(((2 * 32.2 * po * (p1 - p2) * 14.223) / (1 - np.power(m, 2))), 0.5)
    qo = qo * (86400 / (144 * 5.615)) * (1 / po)
    return qo - objetivo_qo 

def gastoEnEstrangulador(alfa, m, d_orificio, p1, p2, po):
    qo = alfa * (np.pi * np.power(d_orificio, 2) / 4) * np.power(((2 * 32.2 * po * (p1 - p2) * 14.223) / (1 - np.power(m, 2))), 0.5)
    qo = qo * (86400 / (144 * 5.615)) * (1 / po)
    return qo

def clearDataFrame(df, first, last):
    if not isinstance(df, pd.DataFrame):
        raise ValueError("El argumento 'df' debe ser un DataFrame de Pandas.")
    if first < 0 or last >= len(df.columns):
        raise IndexError("Los índices de las columnas están fuera del rango del DataFrame.")
    df.iloc[:, first:last + 1] = np.nan
    return df

def Objetivo(df_perfil, df_HB, values, const):
    clearDataFrame(df_perfil, 1, 1)
    clearDataFrame(df_perfil, 4, 4)
    delta_p = values['pws'] - values['pwf']
    
    for row in range(0, int(values['orificios'])):
        lista = df_perfil.iloc[row].tolist()
        lista[1] = values['pws'] - delta_p * (values['orificios'] - row)
        lista[2] = values['pws'] - lista[1]
        lista[5] = lista[1] / values['pws']
        df_perfil.iloc[row] = lista
    
    const['pwh'] = df_perfil.iloc[0, 1]
    contador = 1
    while True:
        print(f'iteraciones realizadas {contador}')
        df_perfil = Perfil(df_perfil, df_HB, values, const) 
        df_HB = HagdeornAndBrown(const, df_HB, values['orificios'])
        
        for row in range(0, int(values['orificios'])):
            lista_perfil = df_perfil.iloc[row].tolist()
            lista_HB = df_HB.iloc[row].tolist()
            lista_perfil[1] = lista_HB[4]
            lista_perfil[2] = values['pws'] - lista_perfil[1]
            lista_perfil[5] = lista_perfil[1] / values['pws']
            df_perfil.iloc[row] = lista_perfil
            
        x = values['pwf'] - const['Pwf']
        if x < 0.0001:
            break
        const['pwh'] = const['pwh'] + x / 2
        contador += 1
            
    df_perfil = Perfil(df_perfil, df_HB, values, const)
    print(f'Finalizo el proceso con {contador} iteraciones')
    return df_perfil


def Perfil(df, df_HB, values, data):
    d_min, d_max = 0.0001, 3.0  
    for row in range(0, int(values['orificios'])):
        lista = df.iloc[row].tolist()
        sol = root_scalar(ecuacion, bracket=[d_min, d_max], method='brentq', args=(values['D'], lista[3], lista[1], values['CD'], values['po'], values['Qo_seccion']))
        lista[4] = sol.root
        lista[6] = np.power((lista[4] / values['D']), 2)
        lista[7] = values['CD'] * lista[5] * ((1 - np.power(lista[6], 2)) / (1 - np.power((lista[5] * lista[6]), 2)))
        lista[8] = lista[4] * 64
        lista[10] = gastoEnEstrangulador(lista[7], lista[6], lista[4], lista[3], lista[1], values['po'])
        df.iloc[row] = lista
    for row in range(0, int(values['orificios'])):
        lista = df.iloc[row].tolist()
        acum = df.iloc[row:, 10].sum()
        lista[9] = acum
        df.iloc[row] = lista
    return df

def HagdeornAndBrown(data, df, length_df):
    Rs = 0
    CNL = 0
    pws = data['pws'] * 14.223
    Tws = data['tws'] * 1.8 + 32
    pwh = data['pwh'] * 14.223
    
    clearDataFrame(df, 3, 41)
    dro = 141.5 / (131.5 + data['api'])
    gradiente = 0
    p = pws
    t = Tws
    
    aux = p
    if p <= pws:
        p = pws
    
    drg = data['og']
    Rs = Lasater(p, t, drg, dro, data['api'], Rs)
    if Rs > (38.38 * 5.615):
       Rs =  38.38 * 5.615
    RGA = Rs
    data['RGA'] = RGA / 5.615
    
    if aux <= pws:
        p = aux   
    
    p = pwh
    
    for i in range(1, int(length_df) + 1):
        lista = df.iloc[i - 1].tolist()
        ql = lista[42]
        D = lista[1] / 12
        Rug = lista[2] / 12 
        L = 0 if i == 1 else ((df.iloc[i - 1, 0] - df.iloc[i - 2, 0]) / 0.3048)
        
        lista[3] = L * 0.3048
        
        p = p + gradiente * L
        lista[4] = p / 14.223
        
        if p < 0:
            print('Se requiere mayor presion de fondo')
            return 0
    
        Twh = (Tws - 130)
        Twh = Twh * 1.8 + 32
        dt = 0
        if i > 1: 
            dt = (Tws - Twh) / ((df.iloc[int(length_df) - 1,0] - df.iloc[0,0]) / 0.3048)
            
        t = Twh if i == 1 else t + dt * L
        lista[5] = (t - 32) / 1.8 
        
        drg = data['og']
        Rs = Lasater(p, t, drg, dro, data['api'], Rs)
        if Rs > (38.38 * 5.615): Rs =  38.38 * 5.615
        if Rs > RGA: Rs = RGA 
        lista[6] = Rs / 5.615
        
        Bo = 0.972 + 0.00014 * np.power((Rs * np.power((drg / dro), 0.5) + 1.25 * t), 1.175)
        
        P1 = 0.908758
        P2 = 0.0617551
        Bo = Bo * P1 + P2
        lista[7] = Bo
        
        denoil = (62.4 * dro + 0.0764 * drg * Rs / 5.615) / Bo
        lista[8] = denoil
        
        visoil = 0
        visoil = Beggs_Robinson(t, Rs,  data['api'], visoil)
        lista[9] = visoil
        
        aux = 141.5 / (denoil / 62.4) - 131.5
        tensioninterfacial_g_o = (39 - 0.2571 * aux) - ((t - 68) * (1.5)) / 32
        if t < 68: tensioninterfacial_g_o = 39 - 0.2571 * aux
        if t > 100: tensioninterfacial_g_o = 37.5 - 0.2571 * aux 
        
        C = 1 - 0.02 * np.power((p + 14.7), 0.5)
        if C < 0: C = 1 - 0.02 * np.power((p + 14.7), 0.45)
        tensioninterfacial_g_o = C * tensioninterfacial_g_o
        lista[10] = tensioninterfacial_g_o
        
        z = 0
        z = Factor_z(t, p, drg, z)
        lista[11] = z
        
        
        Bg = 0.028269 * z * (t + 460) / (p + 14.7)
        lista[12] = Bg
        
        dengas = 0.0764 * drg / Bg
        lista[13] = dengas
        dengas = dengas * data['fact_corr_fricc']
        
        visgas = 0
        visgas = Lee(t, dengas, visgas)
        lista[14] = visgas
        
        Bw = (0.9911 + 0.0000635 * t + 0.00000085 * np.power(t, 2)) + (0.000001093 - 0.000000003497 * t + 0.00000000000457 * np.power(t, 2)) * p + (-0.00000000005 + 6.429E-13 * t - 1.43E-15 * np.power(t, 2)) * np.power(p, 2)
        Bw = Bw * (1 + (0.000000051 * p + (t - 60) * (0.00000574 - 0.000000000195 * p) + np.power((t - 60), 2) * (-0.0000000323 + 0.00000000000085 * p)) * data['sal_agua'] * 0.0001)
        lista[15] = Bw
        
        viswat = ((-0.04518 + 0.0000009313  * data['sal_agua'] - 0.00000000000393 * np.power(np.int64(data['sal_agua']), 2)) + (70.634 + 0.0000000009576 * np.power(np.int64(data['sal_agua']), 2)) / t) * (1 + 0.00000035 * np.power(p, 2) * (t - 40))
        lista[16] = viswat
    
        m = (1 / (1 + data['WOR'])) * dro * 5.615 * 62.4 + 0.0764 * drg  * Rs + (data['WOR'] / (1 + data['WOR'])) * data['ow'] * 5.615 * 62.4
        denliq = ((1 / (1 +  data['WOR'])) * dro * 62.4 + 0.0764 * drg * Rs / 5.615) / Bo + ( data['WOR'] / (1 +  data['WOR'])) * data['ow'] * 62.4 / Bw
        denliq = denliq * data['fact_corr_fricc']
        Ap = (np.pi * np.power(D, 2)) / 4
        Qo = ql * (1 / (1 + data['WOR'])) * Bo * (5.615 / 86400)
        Qw = ql * (data['WOR'] / (1 + data['WOR'])) * Bw * (5.615 / 86400)
        Qg = ql * (1 / (1 + data['WOR'])) * Bg * (RGA - Rs) / 86400
        
        if data['piny'] > lista[0] and data['Qgi'] > 0:
            drg = data['og_lift']
            z = Factor_z(t, p, drg, z)
            Bg = 0.028269 * z * (t + 460) / (p + 14.7) 
            dengas = (dengas * Qg + (0.0764 * drg / Bg) * (data['Qgi'] * (1000000 / 86400)))
            Qg = Qg + data['Qgi'] * (1000000 / 86400)
            dengas = dengas / Qg 
        
        vsl = (Qo + Qw) / Ap
        vsg = Qg / Ap
        visliq = visoil * (1 / (1 + data['WOR'])) + viswat * (data['WOR'] / (1 + data['WOR']))
        NLV = 1.938 * vsl * np.power((denliq / tensioninterfacial_g_o), 0.25)
        NGV = 1.938 * vsg * np.power((denliq / tensioninterfacial_g_o), 0.25)
        ND = 120.872 * D * np.power((denliq / tensioninterfacial_g_o), 0.5)
        NL = 0.15726 * visliq * np.power((1 / (denliq * np.power(tensioninterfacial_g_o, 3))), 0.25)
        vs = 0.8 * 0.3048
        vm = (vsl + vsg) * 0.3048
        a = 1.071 - 0.2218 * np.power(vm, 2) / (D * np.power(0.3048, 2))
        if a < 0.13: a = 0.13
        B = vsg / vm
        Hl = 0
        CNL = 0
        HLPSI = 0
        psi = 0
        if (B - a) < 0:
            Hl = 1 - 0.5 * (1 + (vm / vs) - np.power((np.power((1 + vm / vs), 2) - 4 * vsg / vs), 0.5))
        else:
            CNL = Lagrange_Factor_Numero_Viscosidad(NL, CNL)
            if vsg > 0:
                aux = (NLV / np.power(NGV, 0.575)) * np.power((p / 14.7), 0.1) * (CNL / ND)
                aux = 0.00326 * (np.power(ql, 0.425) * ((1 / (1 + data['WOR'])) * Bo + (data['WOR'] / (1 + data['WOR'])) * Bw) * np.power(p, 0.675) * np.power(tensioninterfacial_g_o, 0.39375) * CNL) / (np.power(D, 1.85) * np.power((t * z), 0.575) * np.power(denliq, 0.39375) * np.power((RGA - Rs * (1 / (1 +  data['WOR']))), 0.575))
            else:
                aux = 1
            if lista[0] < data['piny'] and data['Qgi'] > 0:
                aux = 0.00326 * (np.power(ql, 0.425) * ((1 / (1 + data['WOR'])) * Bo + (data['WOR'] / (1 + data['WOR'])) * Bw) * np.power(p, 0.675) * np.power(tensioninterfacial_g_o, 0.39375) * CNL) / (np.power(D, 1.85) * np.power((t * z), 0.575) * np.power(denliq, 0.39375) * np.power((RGA - Rs * (1 / (1 + data['WOR'])) + (data['Qgi'] * 1000000 * vsl / ql)), 0.575))
            HLPSI = np.power(((0.0047 + 1123.32 * aux + 729489.64 * np.power(aux, 2)) / (1 + 1097.1566 * aux + 722153.97 * np.power(aux, 2))), 0.5)
            aux = NGV * np.power(NL, 0.38) / np.power(ND, 2.14)
            psi = (1.0886 - 69.9473 * aux + 2334.3497 * np.power(aux, 2) - 12896.683 * np.power(aux, 3)) / (1 - 53.4401 * aux + 1517.9369 * np.power(aux, 2) - 8419.8115 * np.power(aux, 3))
            if NGV == 0: 
                psi = 1.00
            Hl = HLPSI * psi
        if Hl > 1:
            Hl = 1
        
        denmez = denliq * Hl + dengas * (1 - Hl)
        
        Nre = 0.022 * ql * m / (D * np.power(visliq, Hl) * np.power(visgas, (1 - Hl)))
        f = 0
        f = Friccion(Nre, D, Rug, f)
        gradiente = (1 / 144) * denmez
        lista[38] = gradiente

        gradiente = (1 / 144) * (f * np.power(ql, 2) * np.power(m, 2) / ((2.9652 * np.power(np.int64(10), 11)) * np.power(D, 5) * denmez))
        lista[39] = gradiente
        
        gradiente = (1 / 144) * (f * np.power(ql, 2) * np.power(m, 2) / ((2.9652 * np.power(np.int64(10), 11)) * np.power(D, 5) * denmez))
        lista[40] = gradiente
        
        lista[17] = m
        lista[18] = Ap
        lista[19] = Qo
        lista[20] = Qw
        lista[21] = Qg
        lista[22] = vsl
        lista[23] = vsg
        lista[24] = visliq
        lista[25] = NLV
        lista[26] = NGV
        lista[27] = ND
        lista[28] = NL
        lista[29] = CNL
        lista[30] = (NLV / np.power(NGV, 0.575)) * np.power((p / 14.7), 0.1) * (CNL / ND) if vsg > 0 else 1
        lista[31] = HLPSI
        lista[32] = NGV * np.power(NL, 0.38) / np.power(ND, 2.14)
        lista[33] = psi
        lista[34] = Hl
        lista[35] = denmez
        lista[36] = Nre
        lista[37] = f
        lista[41] = gradiente * (1 / (14.22 * 0.3048))
        df.iloc[i - 1] = lista
    for i in range(2, int(length_df) + 1):
        lista = df.iloc[i - 1].tolist()
        lista[4] = df.iloc[(i-2), 4] + lista[3] * lista[41]
        df.iloc[i - 1] = lista
    data['Pwf'] = df.iloc[:, 4].max()
    return df