# -*- coding: utf-8 -*-
"""
Created on Sun Dec  5 20:28:33 2021

@author: User
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


data = pd.read_csv('Convergencia_Malha.csv',skiprows = 2)

x_range = np.linspace(0,11000, 10)

f_tensao = interp1d(data.Total, data.Mpa, fill_value = "extrapolate")
f_error = interp1d(data.Total, data.Error, fill_value = "extrapolate")

interp_tensao = f_tensao(x_range)
interp_error = f_error(x_range)

#%% Tensão

fig_tensao = plt.figure(figsize= (10,5))
# plt.scatter(data.Total, data.Mpa)
plt.plot(x_range, interp_tensao, color = 'k')
plt.grid()
plt.xlabel("Total de elementos (CQUAD4 + CBAR)", fontsize = 12)
plt.ylabel("Max Shear [MPa]", fontsize = 12)
# plt.savefig("convergencia_malha_tensao.pdf")
plt.show()

#%% Error
fig_error = plt.figure(figsize=(10,5))
# plt.scatter(data.Total, data.Error)
plt.plot(x_range, interp_error, color = 'k')
plt.grid()
plt.xlabel("Total de elementos (CQUAD4 + CBAR)", fontsize = 12)
plt.ylabel("Error [%]", fontsize = 12)
# plt.savefig("convergencia_malha_error.pdf")
plt.show()
