import sympy as sp



# Definiere die Symbole
x1, x2, x3, x4, x5, x6, x7, t, y, phi1, phi2 = sp.symbols('x1 x2 x3 x4 x5 x6 x7 t y phi1 phi2')

# Definiere die Funktionen
#f1 = sp.sin(phi1) * x1 + sp.sin(phi1+phi2) * y
#f2 = sp.cos(phi1) * x1 + sp.cos(phi1+phi2) * y
f3 = x1 * sp.exp(-x2 * t) * sp.sin(x3 * t + x4)

g1 = x1 / (1 + ((t - x5) / (x7*0.5))**2)
g2 = x2
g3 = x3*t
g4 = x4*t**2
g5 = x5*t**3


# Erstelle eine Liste der Funktionen
functions = [g1, g2, g3, g4, g5]

def jacobi(f, x):
    # Berechne die Jacobimatrix
    J = sp.Matrix([[sp.diff(f, var) for var in x] for f in functions])
    return J

print("Jacobi:\n",jacobi(functions, (x1, x2, x3, x4, x5, x6, x7)))
