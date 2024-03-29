import sympy as sp



# Definiere die Symbole
x, y, phi1, phi2 = sp.symbols('x y phi1 phi2')

# Definiere die Funktionen
f1 = sp.sin(phi1) * x + sp.sin(phi1+phi2) * y
f2 = sp.cos(phi1) * x + sp.cos(phi1+phi2) * y

# Erstelle eine Liste der Funktionen
functions = [f1, f2]

def jacobi(f, x):
    # Berechne die Jacobimatrix
    J = sp.Matrix([[sp.diff(f, var) for var in x] for f in functions])
    return J

print("Jacobi 2",jacobi(functions, (phi1, phi2)))
