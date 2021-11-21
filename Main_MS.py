"""
Nome: Alessandro Melo de Oliveira
NºUSP: 10788662

Script principal para calcular as margens de segurança.
"""

import sympy as sp

from PropriedadesLaminado import InitLaminado
from Calcula_EpsilonSigma_TensaoDeformacao import Calcula_EpsilonSigma_TensaoDeformacao
from CalculaCriterios import Criterios


#%% Propriedades gerais do laminado
init_matrices = InitLaminado()
A_linha, C_linha, Q_barra_list, T_angle_list = init_matrices.get_values()

#%% Cálculo das margens de segurança
N1 = sp.Matrix([100000,0,0]) #Caso de carregamento 1
N2 = sp.Matrix([27000,0,0]) #Caso de carregamento 2
casos_de_carga = [N1, N2]


for (caso_carga, i) in zip(casos_de_carga, range(0, len(casos_de_carga))):
    print("Caso de carga: {}".format(i+1))
    
    EpsilonSigma = Calcula_EpsilonSigma_TensaoDeformacao(caso_carga, A_linha, C_linha, Q_barra_list, T_angle_list)
    tensao_deformacao_local_list = EpsilonSigma.get_values_tensao_deformacao()


    CalculaCriterio = Criterios(tensao_deformacao_local_list)
    margens_seguranca = CalculaCriterio.get_margens()
    CalculaCriterio.extract_min(display = True)
    print("##############################3")
