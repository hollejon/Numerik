import numpy as np
import matplotlib.pyplot as plt



f =     10**5       # Frequenz in HZ
R1 =    800         # Primärwiderstand in Ohm
Lp =    50*10**-6   # Primärinduktivität in H    
Lps =   150*10**-6  # Streuinduktivität in H
Lsp =   150*10**-6  # Streuinduktivität in H
R2 =    2           # Sekundärwiderstand in Ohm
Ls =    500*10**-6  # Sekundärinduktivität in H
U0 =    4           # Leerlaufspannung in V



# Anfangsbedingungen




# DGL System Ax =b
A = np.array([[(R1*1/Lp)/(1-Lps*Lsp*1/Ls),(R2*1/Ls*Lps)/(1-Lps*Lsp*1/Ls)],[(R2*(1/Ls)*Lsp*Lps)/(1+U0*(1/Lp)*Lsp*Lps),(R2(1/Ls)*Lsp)/(1+U0*(1/Lp)*Lsp*Lps*(1/Ls))]])
b = np.array([(U0*(1/Ls))/(1-Lps*Lsp*(1/Ls)),0])


