import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from funciones import Perfil, Objetivo, HagdeornAndBrown
import time

data = {
       'api': 15,
       'og': 0.73,
       'ow': 1.026,
       'Ql': 0,
       'WOR': 0,
       'pws': 350,
       'tws': 117,
       'piny': 0,
       'Qgi': 0,
       'og_lift': 0.675,
       'pwh': 0,
       'fact_corr_fricc': 1,
       'fact_corr_colg': 1,
       'sal_agua': 120000,
       'RGA': 0,
       'Pwf': 0
}

init_val = {
        'CD': 0.8,
        'D': 3.958,
        'po': 51.168,
        'pws': 350,
        'IP': 175,
        'Qo': 2000,
        'orificios': 0,
        'pwf': 0,
        'Qo_seccion': 0
}

tuberia = {
    'longitud': 500,
    'separacion': 5
}

def graficarValores(x, y, title, x_label, y_label):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, y)
    ax.set(title = title, xlabel = x_label, ylabel = y_label)
    ax.grid()
    plt.show()


def init_df(df, values, tub):
    df['L'] = np.arange(0, tub['longitud'] + tub['separacion'], 5, dtype=np.float64)
    df['Pws'] = np.full(values['orificios'].as_integer_ratio(), values['pws'])
    pws = values['pws']
    pwf = values['pwf']
    delta_p = pws - pwf
    for row in range(0, int(values['orificios'])):
        lista = df.iloc[row].tolist()
        lista[1] = pws - delta_p * (values['orificios'] - row)
        lista[2] = pws - lista[1]
        lista[5] = lista[1] / pws
        df.iloc[row] = lista
    return 0



#Se calculan los valores iniciales de la distancia entre orificio, el gasto en cada seccion ideal y la presion de fondo con respecto al 
#indice de produccion del pozo
init_val['orificios'] = tuberia['longitud'] / tuberia['separacion'] + 1
init_val['Qo_seccion'] = init_val['Qo'] / init_val['orificios']
init_val['pwf'] = init_val['pws'] - init_val['Qo_seccion'] / init_val['IP'] / 14.223

#Generacion de los titulos de las columnas de los diferentes DataFrames a generar 
df_hagedorn_and_brown = pd.DataFrame(columns=['Prof', 'Diametro', 'Rug', 'L', 'Presion', 'Temp', 'Rs', 'Bo', 'po', 'uo',
                                              'ogo', 'z', 'Bg', 'pg', 'ug', 'Bw', 'uw', 'M', 'Ap', 'Qo', 
                                              'Qw', 'Qg', 'Vsl', 'Vsg', 'ul', 'NLV', 'NVG', 'ND', 'NL', 'CNL', 
                                              'NLv*P', 'HL/w', 'Ngv*NL', 'w', 'HL', 'pm', 'Nre', 'f', 'Deltap/Deltah_elevacion', 'Deltap/Deltah_friccion', 
                                              'Deltap/Deltah_psi', 'Deltap/Deltah_kg', 'Ql'])

df_perfil = pd.DataFrame(columns=['L', 'Pwf', 'DeltaP', 'Pws', 'd', 'X', 'M', 'a', 'd(64")', 'QoAcum', 'qo'])

init_df(df_perfil, init_val, tuberia)

df_hagedorn_and_brown['Prof'] = np.arange(0, tuberia['longitud'] + tuberia['separacion'], 5, dtype=np.float64)
df_hagedorn_and_brown['Diametro'] = np.full(init_val['orificios'].as_integer_ratio(), init_val['D'])
df_hagedorn_and_brown['Rug'] = np.full(init_val['orificios'].as_integer_ratio(), 0.0006)

#Se realizan los calculos necesarios para determinar el diametro del orificio en cada seccion para que se genere un gasto por seccion 
#previamente calculado en condiciones 
inicio = time.time()
df_perfil = Perfil(df_perfil, df_hagedorn_and_brown, init_val, data)
df_hagedorn_and_brown.iloc[:, 42] = df_perfil.iloc[:, 9].copy()
data['pwh'] = df_perfil.iloc[0, 1]

#df_perfil = Objetivo(df_perfil, df_hagedorn_and_brown, init_val, data)
df_hagedorn_and_brown = HagdeornAndBrown(data, df_hagedorn_and_brown, init_val['orificios'])
df_perfil = Objetivo(df_perfil, df_hagedorn_and_brown, init_val, data)
#print('saliendo de objetivo')
fin = time.time()
print(fin - inicio)
df_hagedorn_and_brown.to_csv('test-2-02-25.csv', index=False)
df_perfil.to_csv('testBPD-01.csv', index=False)

#Se grafican los resultados de los calculos realizados

gasto_redondeado = np.round(df_perfil.iloc[:, 10].tolist(), 4).tolist()
presion_redondeada = np.round(df_perfil.iloc[:, 1].tolist(), 4).tolist()
diametro_redondeado = np.round(df_perfil.iloc[:, 4].tolist(), 4).tolist()
graficarValores(df_perfil.iloc[:,0].tolist(), presion_redondeada, 'Presión', 'Distancia (m)', 'Presión (psi)')
graficarValores(df_perfil.iloc[:,0].tolist(), diametro_redondeado, 'Diametro Orificio', 'Distancia (m)', 'Diametro (plg)')
graficarValores(df_perfil.iloc[:,0].tolist(), gasto_redondeado, 'Gasto de Produccion', 'Distancia (m)', 'Gasto (bpd)')
