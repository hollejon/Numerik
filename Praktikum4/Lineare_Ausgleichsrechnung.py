import numpy as np
from scipy.linalg import solve_triangular
import matplotlib.pyplot as plt
from numpy.linalg import cholesky


data = np.genfromtxt('data.txt')

L = cholesky(A)


class myFit:
  def __init__(self, data):
      self.setData(data)
      self.c = None       # Klassen Variable für die Modell Koeffizienten
      self.omega = omega

  def setData(self, data):
      self.ti = data[:,0] # Zeitstempel
      self.ui = data[:,1] # Messwerte

  def computeCoefs(self, n=5):
      self.n = n

    #   <<snipp>>

    #   self.c = <<snipp>>

  def compute(self,t):
      if not type(mf.c) == np.ndarray:
          self.computeCoefs()

    #   y = <<snipp>>

      return y
  
# Objekt instanzieren
mf = myFit(data)
# Koeffizienten berechnen
mf.computeCoefs(5)
# Visualisieren
plt.plot(t,mf.compute(t))
plt.show()