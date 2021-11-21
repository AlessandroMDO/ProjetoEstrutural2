"""
Nome: Alessandro Melo de Oliveira
NºUSP: 10788662

Script para calcular as margens de segurança.
"""

import sympy as sp

class Criterios:
    
    def __init__(self, dict_deformacao_tensao):
        
        self.dict_deformacao_tensao = dict_deformacao_tensao
        
        #[MPa]
        self.XT = 1400
        self.YT = 47
        self.XC = -930
        self.YC = -130
        self.S12 = 53

        #[m]
        self.XT_linha = 0.025
        self.YT_linha = 0.2/100
        self.XC_linha = -1/100
        self.YC_linha = -1.8/100
        self.S12_linha = 3/100
        
        
    #Critério da máxima tensão
    def max_tensao(self, lamina):
    
        sigma_1 = lamina[0] 
        sigma_2 = lamina[1]
        sigma_12 = lamina[2]
        
        if sigma_1 > 0:
            MS_1 = self.XT/sigma_1 - 1
        else:
            MS_1 = self.XC/sigma_1 - 1 
            
        if sigma_2 > 0:
            MS_2 = self.YT/sigma_2 - 1
        else:
            MS_2 = self.YC/sigma_2 - 1
            
        MS_12 = self.S12/abs(sigma_12) - 1
        
        return sp.Matrix([MS_1,MS_2,MS_12])
    
    #Critério da máxima deformacao
    def max_deformacao(self, lamina):
    
        epsilon_1 = lamina[0] 
        epsilon_2 = lamina[1]
        epsilon_12 = lamina[2]
                
        if epsilon_1 > 0:
            MS_1 = self.XT_linha/epsilon_1 - 1
        else:
            MS_1 = self.XC_linha/epsilon_1 - 1 
            
        if epsilon_2 > 0:
            MS_2 = self.YT_linha/epsilon_2 - 1
        else:
            MS_2 = self.YC_linha/epsilon_2 - 1
            
        MS_12 = self.S12_linha/abs(epsilon_12) - 1
        
        return sp.Matrix([MS_1,MS_2,MS_12])
    
    
    #Critério de Tsai-Wu
    def tsai_wu(self, lamina):
    
        sigma_1 = lamina[0] 
        sigma_2 = lamina[1]
        sigma_12 = lamina[2]
        
        
        F1 = 1/self.XT + 1/self.XC
        F11 = -1/(self.XT*self.XC)
        F2 = 1/self.YT + 1/self.YC
        F22 = -1/(self.YT*self.YC)
        F66 = (1/self.S12)**2
        F12 = -0.5*(F11*F22)**.5
        
        A = F11*sigma_1**2 + F22*sigma_2**2 + F66*sigma_12**2 + 2 *F12*sigma_1*sigma_2
        B = F1*sigma_1 + F2*sigma_2
        
        Sf_menos = abs((-B - (B**2 + 4*A)**.5)/(2*A))
        Sf_mais = (-B + (B**2 + 4*A)**.5)/(2*A)
        
        MS_menos =  Sf_menos - 1
        MS_mais = Sf_mais - 1
        
        return MS_mais
    
    #Critério de Tsai-Hill
    def tsai_hill(self, lamina):
    
        sigma_1 = lamina[0] 
        sigma_2 = lamina[1]
        sigma_12 = lamina[2]
        
        if sigma_1 >0 and sigma_2 > 0:
            self.f_sigma = (sigma_1**2)/self.XT**2 + (sigma_2**2)/self.YT**2 - (sigma_1*sigma_2)/self.XT**2 + (sigma_12**2)/self.S12**2
        elif sigma_1>0 and sigma_2<0:
            self.f_sigma = (sigma_1**2)/self.XT**2 + (sigma_2**2)/self.YC**2 - (sigma_1*sigma_2)/self.XT**2 + (sigma_12**2)/self.S12**2
        elif sigma_1<0 and sigma_2<0:
            self.f_sigma = (sigma_1**2)/self.XC**2 + (sigma_2**2)/self.YC**2 - (sigma_1*sigma_2)/self.XC**2 + (sigma_12**2)/self.S12**2
        elif sigma_1<0 and sigma_2>0:
            self.f_sigma = (sigma_1**2)/self.XC**2 + (sigma_2**2)/self.YT**2 - (sigma_1*sigma_2)/self.XC**2 + (sigma_12**2)/self.S12**2
        
        self.FS = self.f_sigma**.5
        self.MS = 1/self.FS - 1
        
        return self.MS
        
        
    #Função que calcula todas as margens de segurança, para todas as lâminas em cada uma das três posições
    def get_margens(self):
        self.margens_seguranca_list = []
        for lamina in self.dict_deformacao_tensao:
            
            dict_margens = {"botton" : {"max_tensao": {}, "max_deformacao": {}, "tsai_wu": {}, "tsai_hill":{}},
            "mid":{"max_tensao": {}, "max_deformacao": {}, "tsai_wu": {}, "tsai_hill":{}},
            "top":{"max_tensao": {}, "max_deformacao": {}, "tsai_wu": {}, "tsai_hill":{}}}
            
            for posicao in ['botton', 'mid', 'top']:
                dict_margens[posicao]['max_tensao'] = self.max_tensao(lamina['tensao'][posicao])
                dict_margens[posicao]['max_deformacao'] = self.max_deformacao(lamina['deformacao'][posicao])
                dict_margens[posicao]['tsai_wu'] = self.tsai_wu(lamina['tensao'][posicao])
                dict_margens[posicao]['tsai_hill'] = self.tsai_hill(lamina['tensao'][posicao])
            self.margens_seguranca_list.append(dict_margens)
            
        return self.margens_seguranca_list
    
    #Função que extrai as menores margens de segurança entre aquelas calculadas na função acima
    def extract_min(self, display = False):
        margens_seguranca_list = self.get_margens()
        
        
        
        #Max tensão
        dict_margens_minimas_max_tensao = {"sigma_1": {},
                                           "sigma_2": {},
                                           "tau_12": {}}
        sigma1_list_max_tensao = []
        sigma2_list_max_tensao = []
        tau12_list_max_tensao = []
        
        
        #Max deformacao
        dict_margens_minimas_max_deformacao = {"epsilon_1": {},
                                           "epsilon_2": {},
                                           "gamma_12": {}}
        sigma1_list_max_deformacao = []
        sigma2_list_max_deformacao = []
        tau12_list_max_deformacao = []
        
        #Tsai-Hill
        ms_list_tsai_hill = []
        
        
        #Tsai-Wu
        ms_list_tsai_wu = []
        
        
        for lamina_margens in margens_seguranca_list:
            for posicao in ['top', 'mid', 'botton']:
                sigma1_list_max_tensao.append(lamina_margens[posicao]['max_tensao'][0])
                sigma2_list_max_tensao.append(lamina_margens[posicao]['max_tensao'][1])
                tau12_list_max_tensao.append(lamina_margens[posicao]['max_tensao'][2])
                
                sigma1_list_max_deformacao.append(lamina_margens[posicao]['max_deformacao'][0])
                sigma2_list_max_deformacao.append(lamina_margens[posicao]['max_deformacao'][1])
                tau12_list_max_deformacao.append(lamina_margens[posicao]['max_deformacao'][2])
                
                ms_list_tsai_hill.append(lamina_margens[posicao]['tsai_hill'])
                ms_list_tsai_wu.append(lamina_margens[posicao]['tsai_wu'])
        
        dict_margens_minimas_max_tensao["sigma_1"] = round(min(sigma1_list_max_tensao),2) 
        dict_margens_minimas_max_tensao["sigma_2"] = round(min(sigma2_list_max_tensao),2)
        dict_margens_minimas_max_tensao["tau_12"] = round(min(tau12_list_max_tensao),2)
        
        dict_margens_minimas_max_deformacao["epsilon_1"] = round(min(sigma1_list_max_deformacao),2) 
        dict_margens_minimas_max_deformacao["epsilon_2"] = round(min(sigma2_list_max_deformacao),2)
        dict_margens_minimas_max_deformacao["gamma_12"] = round(min(tau12_list_max_deformacao),2)
        
        ms_tsai_wu = round(min(ms_list_tsai_wu),2)
        ms_tsai_hill = round(min(ms_list_tsai_hill),2)
        
        
        if display == True:
            
            print("Mínimas margens de segurança Max Tensão:")
            print(dict_margens_minimas_max_tensao)
            print("----------------------------------------")
            print("Mínimas margens de segurança Max Deformação:")
            print(dict_margens_minimas_max_deformacao)
            print("----------------------------------------")
            print("Mínima margen de segurança Tsai-Hill:")
            print(ms_tsai_hill)
            print("----------------------------------------")
            print("Mínima margen de segurança Tsai-Wu:")
            print(ms_tsai_wu)
        

        # return dict_margens_minimas_max_tensao, dict_margens_minimas_max_deformacao, ms_tsai_wu, ms_tsai_hill
        
            
        
        
