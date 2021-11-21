"""
Nome: Alessandro Melo de Oliveira
NºUSP: 10788662

Script para calcular as tensões e deformações, tanto no sistema global quanto no local
"""
class Calcula_EpsilonSigma_TensaoDeformacao:
    
    def __init__(self, cargas, A_linha, C_linha, Q_barra_list, T_angle_list):
        self.cargas = cargas
        
        h = 0.1/1000 # Espessura de cada lâmina [m]
        self.h_list = [-4*h, -3*h, -2*h, -h, 0, h, 2*h, 3*h, 4*h]
        
        epsilon0 = A_linha*cargas
        K = C_linha*cargas
        
        sigma_global_list = []
        self.epsilon_local_list = []
        self.sigma_local_list = []
        
        for i in range(0, len(Q_barra_list)):
            
            #Sigma global
            sigma_global_i = Q_barra_list[i]*epsilon0
            sigma_global_list.append(sigma_global_i)
            
            #Sigma local
            sigma_local_i = T_angle_list[i]*sigma_global_list[i]
            self.sigma_local_list.append(sigma_local_i)
            
            #Epsilon local
            epsilon_local_i = T_angle_list[i]*epsilon0
            self.epsilon_local_list.append(epsilon_local_i)
            
        self.tensao_deformacao_local_list = []

                        
        for (Q_barra_X, z, T_X) in zip(Q_barra_list,self.h_list, T_angle_list):
            
            dict_tensao_deformacao_local = {"tensao":{"botton":[], "mid":[], "top":[]}, "deformacao":{"botton":[], "mid":[], "top":[]}}
            dict_tensao_deformacao_global = {"tensao":{"botton":[], "mid":[], "top":[]}, "deformacao":{"botton":[], "mid":[], "top":[]}}
            
            #Posições botton, mid e top
            z_botton = z
            z_mid = z + h/2
            z_top = z + h
            
            # Tensão e deformacao na posição botton
            tensao_global_botton = Q_barra_X*(epsilon0 + z_botton*K)
            tensao_local_botton = T_X*tensao_global_botton
            deformacao_global_botton = epsilon0 + z_botton*K
            deformacao_local_botton = T_X*deformacao_global_botton
            
            # Tensão e deformacao na posição mid
            tensao_global_mid = Q_barra_X*(epsilon0 + z_mid*K)
            tensao_local_mid = T_X*tensao_global_mid
            deformacao_global_mid = epsilon0 + z_mid*K
            deformacao_local_mid = T_X*deformacao_global_mid
            
            # Tensão e deformacao na posição top
            tensao_global_top = Q_barra_X*(epsilon0 + z_top*K)
            tensao_local_top = T_X*tensao_global_top
            deformacao_global_top = epsilon0 + z_top*K
            deformacao_local_top = T_X*deformacao_global_top
            
            
            dict_tensao_deformacao_global['tensao']['botton'] = tensao_global_botton/1e6
            dict_tensao_deformacao_global['tensao']['mid'] = tensao_global_mid/1e6
            dict_tensao_deformacao_global['tensao']['top'] = tensao_global_top/1e6
            dict_tensao_deformacao_global['deformacao']['botton'] = deformacao_global_botton
            dict_tensao_deformacao_global['deformacao']['mid'] = deformacao_global_mid
            dict_tensao_deformacao_global['deformacao']['top'] = deformacao_global_top
            
            
            #Retorna a tensão em [MPa] e a deformação em [m]
            dict_tensao_deformacao_local['tensao']['botton'] = tensao_local_botton/1e6
            dict_tensao_deformacao_local['tensao']['mid'] = tensao_local_mid/1e6
            dict_tensao_deformacao_local['tensao']['top'] = tensao_local_top/1e6
            dict_tensao_deformacao_local['deformacao']['botton'] = deformacao_local_botton
            dict_tensao_deformacao_local['deformacao']['mid'] = deformacao_local_mid
            dict_tensao_deformacao_local['deformacao']['top'] = deformacao_local_top
            
            self.tensao_deformacao_local_list.append(dict_tensao_deformacao_local)
                
    def get_values_epsilon_sigma(self):
        return self.epsilon_local_list, self.sigma_local_list
    
    def get_values_tensao_deformacao(self):
        return self.tensao_deformacao_local_list
            