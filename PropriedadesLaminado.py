"""
Nome: Alessandro Melo de Oliveira
NºUSP: 10788662

Script para calcular as matrizes de referência do lâminado
ou seja, as matrizes A,B,C,D...
"""


from sympy import pi
import sympy as sp

m,n = sp.symbols('m n')

class InitLaminado:
    
    def __init__(self):
    
        E11 = 100e9 #[Pa]
        E22 = 10e9 #[Pa]
        v12 = 0.34
        G12 = 5.4e9 #[Pa]
        G23 = 3.05e9 #[Pa]
        
        h = 0.1/1000 # Espessura de cada lâmina [m]    
        self.h_list = [-4*h, -3*h, -2*h, -h, 0, h, 2*h, 3*h, 4*h]
        
        #Ângulo de cada lâmina
        list_angle = [0, pi/4, -pi/4, pi/2, pi/2, -pi/4, pi/4, 0]
        
        
        Q11_value = (E11**2)/(E11-v12**2 * E22)
        Q21_value =  (v12*E11*E22)/(E11-v12**2 * E22)
        Q22_value = (E11*E22)/(E11-v12**2 * E22)
        Q66_value = G12
        
        Q = sp.zeros(3,3)
        Q[0,0] = Q11_value
        Q[0,1] = Q21_value
        Q[0,2] = 0
        Q[1,0] = Q21_value
        Q[1,1] = Q22_value
        Q[1,2] = 0
        Q[2,0] = 0
        Q[2,1] = 0
        Q[2,2] = Q66_value
        
        Q11, Q12, Q16, Q22, Q26, Q16, Q66 = sp.symbols('Q11 Q12 Q16 Q22 Q26 Q16 Q66')
        Q_barra = sp.zeros(3,3)
        Q_barra[0,0] =  Q11*m**4 + 2*m**2 * n**2 * (Q12 + 2*Q66) + Q22*n**4
        Q_barra[0,1] = (Q11 + Q22 - 4*Q66)*n**2 * m**2 + Q12*(n**4 + m**4)
        Q_barra[0,2] = (Q11 - Q12)*n*(m**3) + (Q12 - Q22)*(n**3)*m - 2*m*n*(m**2 - n**2)*Q66
        Q_barra[1,0] = (Q11 + Q22 - 4*Q66)*n**2 * m**2 + Q12*(n**4 + m**4)
        Q_barra[1,1] = Q11*n**4 + 2*(Q12 + 2*Q66)*n**2 * m**2 + Q22*m**4
        Q_barra[1,2] = (Q11 - Q12)*n**3*m + (Q12 - Q22)*n*m**3 + 2*m*n*(m**2 - n**2)*Q66
        Q_barra[2,0] = (Q11 - Q12)*n*m**3 + (Q12 - Q22)*n**3*m - 2*m*n*(m**2 - n**2)*Q66
        Q_barra[2,1] = (Q11 - Q12)*(n**3)*m + (Q12 - Q22)*n*(m**3) + 2*m*n*(m**2 - n**2)*Q66
        Q_barra[2,2] = (Q11 + Q22 - 2*Q12 - 2*Q66)*n**2 * m**2 + Q66*(n**4 + m**4)
        
        
        T = sp.zeros(3)
        T[0,0] = m**2 
        T[0,1] = n**2
        T[0,2] = 2*m*n
        T[1,0] = n**2
        T[1,1] = m**2
        T[1,2] = -2*m*n
        T[2,0] = -m*n
        T[2,1] = m*n
        T[2,2] = (m**2 - n**2)
        
    
        self.Q_barra_list = []
        self.T_angle_list = []
        
        for angle in list_angle:
            Q_barra_angle = Q_barra.subs({
            m:sp.cos(angle),
            n:sp.sin(angle),
            Q11 :Q11_value,
            Q12 : Q21_value,
            Q66 : Q66_value,
            Q22 : Q22_value
            })
            
            T_angle = T.subs({
                m:sp.cos(angle),
                n:sp.sin(angle)
            })
            
            self.Q_barra_list.append(Q_barra_angle)
            self.T_angle_list.append(T_angle)
            
            
        self.A = sp.zeros(3)
        self.B = sp.zeros(3)
        self.D = sp.zeros(3)
        
        for i in range(0,len(self.Q_barra_list)):
            
            self.A += self.Q_barra_list[i]*(self.h_list[i+1] - self.h_list[i])
            self.B += (self.Q_barra_list[i]*(self.h_list[i+1]**2 - self.h_list[i]**2))/2
            self.D += (self.Q_barra_list[i]*(self.h_list[i+1]**3 - self.h_list[i]**3))/3
            
        A_ast= self.A.inv()
        B_ast= -self.A.inv() * self.B
        C_ast= -B_ast.T
        D_ast= self.D - self.B*self.A.inv()*self.B
        
        self.A_linha = A_ast + (B_ast)*(D_ast.inv())*B_ast.T
        self.B_linha = (B_ast)*(D_ast.inv())
        self.C_linha = self.B_linha
        self.D_linha = D_ast.inv()
        
        
    #Retorna as matrizes A,B,D    
    def get_matrices(self):
        return self.A, self.B, self.D
    
    #Retorna as matrizes A', C', assim como Q_barra e T_list
    def get_values(self):
        return self.A_linha, self.C_linha, self.Q_barra_list, self.T_angle_list
        

        
        
        
    
