import numpy as np
from scipy.linalg import solve_triangular
import matplotlib.pyplot as plt
from numpy.linalg import cholesky


data = np.genfromtxt('data.txt')

T0 = 0.01
f0 = 1/T0
t = data[:,0]
y_m = data[:,1]
y_a = np.sin(2*np.pi*f0*t)

# Erstellen Sie eine Figur
fig, ax = plt.subplots(figsize=(10, 5))
# Messwerte Plot
ax.plot(t, y_m, 'o', label='Messwerte', markersize=2)
# Approx Plot
ax.plot(t, y_a, label='Approximierende Funktion')
# Gitter und legende hinzufügen
ax.grid(visible=True, which='major', axis='both')
ax.set_title('Messwerte und Approximierende Funktion')
# Fügen Sie eine Legende hinzu
ax.legend()
# Zeigen Sie den Plot an
plt.show()


# L = cholesky(A)


# class myFit:
#   def __init__(self, data):
#       self.setData(data)
#       self.c = None       # Klassen Variable für die Modell Koeffizienten
#       self.omega = omega

#   def setData(self, data):
#       self.ti = data[:,0] # Zeitstempel
#       self.ui = data[:,1] # Messwerte

#   def computeCoefs(self, n=5):
#       self.n = n

#     #   <<snipp>>

#     #   self.c = <<snipp>>

#   def compute(self,t):
#       if not type(mf.c) == np.ndarray:
#           self.computeCoefs()

#     #   y = <<snipp>>

#       return y
  
# # Objekt instanzieren
# mf = myFit(data)
# # Koeffizienten berechnen
# mf.computeCoefs(5)
# # Visualisieren
# plt.plot(t,mf.compute(t))
# plt.show()