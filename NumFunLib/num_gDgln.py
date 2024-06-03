"""
                                        _____            _                                                                  
                                       |  __ \          | |                                                                 
  _ __    _   _   _ __ ___       __ _  | |  | |   __ _  | |  _ __                                                           
 | '_ \  | | | | | '_ ` _ \     / _` | | |  | |  / _` | | | | '_ \                                                          
 | | | | | |_| | | | | | | |   | (_| | | |__| | | (_| | | | | | | |                                                         
 |_| |_|  \__,_| |_| |_| |_|    \__, | |_____/   \__, | |_| |_| |_|                                                         
                                 __/ |            __/ |                                                                     
                                |___/            |___/                    
"""
"""
Autor: Raphael Romann
Datum: 30.05.2024   


"""
print("Load: num_gDgln.py")
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve,lu_factor,lu_solve

"""
Uebersicht:

def euler_forward(x0, X, N, f):    
    x0: Startwert x, y
    X: Endwert
    N: Anzahl Schritte
    f: Funktion f(x, y) = y'

    
def euler_backward(x0, X, N, f, df, tol):
    x0: Startwert x, y
    X: Endwert
    N: Anzahl Schritte
    f: Funktion f(x, y) = y'
    df: Jacobi-Matrix der Funktion f
    tol: Toleranz

    
def calc_error(x_data, y_data, model_analytisch):
    x_data: x-Werte - Numerisch
    y_data: y-Werte - Numerisch
    model_analytisch: analytische Funktion

    
# Verfahren von Runge (RK2)                                         Nicht getestet
def runge_explizit_RK2(x0, X, N, f):
    x0: Startpunkt
    X: Endpunkt
    N: Anzahl der Schritte
    f: Funktion

    
# Verfahren von Heun (RK2)                                              Nicht getestet
def heun_explizit(x0, X, N, f):
    x0: Startpunkt
    X: Endpunkt
    N: Anzahl der Schritte
    f: Funktion

    
# Verfahren Runge-Kutta-Verfahren (RK4)                                 Nicht getestet
def runga_kutta_RK4(x0, X, N, f):
    x0: Startpunkt
    X: Endpunkt
    N: Anzahl der Schritte
    f: Funktion


    
# implizite Mittelpunktsregel (nicht für Systeme)                       Nicht getestet
# 0.5 | 0.5 
#     | 1.0
def implizit_mittelpunkt(x0, X, f, df, N, tol=1e-3, A=0.5, B=1, C=0.5):
    # x0: Startpunkt
    # X: Endpunkt
    # f: Funktion
    # df: Jacobi Matrix der Funktion f
    # N: Anzahl der Schritte
    # tol: Toleranz

    

# Implizite Trapezregel (nicht für Systeme)                             Nicht getestet
# 0 | 0
# 1 | 0.5   0.5
#   | 0.5   0.5
def implizit_Trapez(x0, X, f, df, N, tol=1e-3):                     
    # x0: Startpunkt
    # X: Endpunkt
    # f: Funktion
    # df: Jacobi Matrix der Funktion f
    # N: Anzahl der Schritte
    # tol: Toleranz

    

# Verfahren Runge-Kutta-Verfahren (RK4) (für ein autonomes System)      Getestet
    # x0: Anfangswerte [x1, x2, ...]
    # t0: Anfangswert
    # h: Schrittweite
    # tend: Zeit bis zum Ende
    # f: Funktion mit n Gleichungen
    return x, t

    

# Implizite Mittelpunktsregel (für 2 Systeme)                           Nicht getestet
# C | A  =  1/2 | 1/2
#   | B  =      | 1
def implizite_mittelpunkt_2Systeme(x0, h, tend, f, df, tol):
    # x0: Startpunkt
    # X: Endpunkt
    # f: Funktion
    # df: Jacobi Matrix der Funktion f
    # N: Anzahl der Schritte
    # tol: Toleranz

    
# Partielle Differentialgleichungen sind auf Zeile 611 beschrieben

"""   


# ----------------------------------------- Praktikum 8 -----------------------------------------
# Euler-vorwaerts / RK explizit
def euler_forward(x0, X, N, f):
    h = (X-x0[0])/N
    N = N+1
    #print("EXP: Schrittweite: ", h)

    x = np.zeros((N))
    y = np.zeros((N))
    x[0] = x0[0]
    y[0] = x0[1]

    for i in range(1, N):
        x[i] = x[i-1] + h
        k = f(x[i-1], y[i-1])
        y[i] = y[i-1] + h*k

    return x, y


"""
# ========== Example Euler-Vorwaerts / RK-Explizit ==========
def model_analytisch(x):
    return np.exp(-4 * x)

def model1(x, y):
    return -4 * y

x0 = [0, 1] # Startwert x, y
res_x, res_y = euler_forward(x0, 2, 20, model1)

x_an = np.linspace(0,2,100)
plt.plot(x_an, model_analytisch(x_an), '-', label='y(x) Analytisch')
plt.plot(res_x, res_y, 'o-', label='Euler Vorwaerts')
plt.ylabel('y-Achse')
plt.xlabel('x-Achse')
plt.legend()
plt.grid(True)
plt.show()

"""

# Euler-rueckwaerts / RK implizit
def euler_backward(x0, X, N, f, df, tol):
    # x0: Startpunkt
    # X: Endpunkt
    # f: Funktion
    # df: Jacobi Matrix der Funktion f
    # N: Anzahl der Schritte
    # tol: Toleranz
    max_iter=20

    h= (X-x0[0])/N
    N = N+1
    #print( "IMP: Schrittweite: ", h)
    x = np.zeros(N)
    y = np.zeros(N)
    x[0] = x0[0]
    y[0] = x0[1]

    #Anfangspunkt berechnen
    k = f(x[0],y[0])

    for i in range(1,N):
        step = 0
        x[i] = x[i-1] + h
        r = k - f(x[i-1] + h, y[i-1] + (h*k))

        while np.abs(r) > tol and step < max_iter:
            j = df(x[i-1] + h, y[i-1] + (h*k))
            delta_k = -r / (1 - (h*j))
            k = k + delta_k
            r = k - f(x[i-1] + h, y[i-1] + (h*k))
            step += 1
        y[i] = y[i-1] + (h*k)
    return x,y

"""
# ========== Example: Euler-Rueckwaerts / RK-Implizit ==========
def model2(x,y):
    return -4*y

def model2_df(x,y):
    return -4

def model_analytisch(x):
    return np.exp(-4 * x)

x0 = (0, 1)
x_imp,y_imp = euler_backward(x0, 2, 10, model2, model2_df ,1e-3)

x_an = np.linspace(0,2,100)
plt.plot(x_imp, y_imp, 'o-', label='Implizites Euler-Verfahren')
plt.plot(x_an, model_analytisch(x_an), '-', label='y(x) Analytisch')
plt.ylabel('y-Achse')
plt.xlabel('x-Achse')
plt.title('Schrittweite: 0.2')
plt.legend()
plt.grid(True)
plt.show()
"""

# Berechnet den Fehler zwischen numerischer und analytischer Lösung

def calc_error(x_data, y_data, model_analytisch):
    r = []
    for i in range(len(x_data)):
        print(x_data[i])
        x_anal = model_analytisch(x_data[i])
        r.append(abs(y_data[i] - x_anal))
    return r

#plt.semilogy(x_num, calc_error(x_num, y_num), 'o-', label='Fehler')





# ----------------------------------------- Praktikum 9 -----------------------------------------
# Verfahren von Runge (RK2)
def runge_explizit_RK2(x0, X, N, f):
    # x0: Startpunkt
    # X: Endpunkt
    # N: Anzahl der Schritte
    # f: Funktion
    h = (X-x0[0])/N
    N = N+1
    #print("EXP: Schrittweite: ", h)

    x = np.zeros((N))
    y = np.zeros((N))
    x[0] = x0[0]
    y[0] = x0[1]

    for i in range(1, N):
        x[i] = x[i-1] + h
        k1 = f(x[i-1], y[i-1])
        k2 = f(x[i-1] + h*0.5, y[i-1] + h*0.5*k1)
        y[i] = y[i-1] + h*k2

    return x, y


# Verfahren von Heun (RK2)
def heun_explizit(x0, X, N, f):
    # x0: Startpunkt
    # X: Endpunkt
    # N: Anzahl der Schritte
    # f: Funktion
    h = (X-x0[0])/N
    N = N+1
    #print("EXP: Schrittweite: ", h)

    x = np.zeros((N))
    y = np.zeros((N))
    x[0] = x0[0]
    y[0] = x0[1]

    for i in range(1, N):
        x[i] = x[i-1] + h
        k1 = f(x[i-1], y[i-1])
        k2 = f(x[i-1] + h, y[i-1] + h*k1)
        y[i] = y[i-1] + h*(0.5*k1 + 0.5*k2)

    return x, y


# Verfahren Runge-Kutta-Verfahren (RK4)
def runga_kutta_RK4(x0, X, N, f):
    # x0: Startpunkt
    # X: Endpunkt
    # N: Anzahl der Schritte
    # f: Funktion
    h = (X-x0[0])/N
    N = N+1
    #print("EXP: Schrittweite: ", h)

    x = np.zeros((N))
    y = np.zeros((N))
    x[0] = x0[0]
    y[0] = x0[1]

    for i in range(1, N):
        x[i] = x[i-1] + h
        k1 = f(x[i-1], y[i-1])
        k2 = f(x[i-1] + 0.5*h, y[i-1] + 0.5*h*k1)
        k3 = f(x[i-1] + 0.5*h, y[i-1] + 0.5*h*k2)
        k4 = f(x[i-1] + h, y[i-1] + h*k3)
        y[i] = y[i-1] + h*((1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4)

    return x, y


# ========== Example: RK2 / Heun / RK4 ==========
# def model1(x, y):
#     return -4 * y
# 
# def model_analytisch(x):
#     return np.exp(-4 * x)
# 
# x_an = np.linspace(0,2,100)
# 
# x0 = [0, 1] # Startwert x, y
# R_RK2_x, R_RK2_y = runge_explizit_RK2(x0, 2, 10, model1)
# H_RK1_x, H_RK1_y = heun_explizit(x0, 2, 10, model1)
# RK_RK4_x, RK_RK4_y = runga_kutta_RK4(x0, 2, 10, model1)
# 
# plt.plot(R_RK2_x, R_RK2_y, 'o-', label='Runge (p=2)')
# plt.plot(H_RK1_x, H_RK1_y, 'o-', label='Heun (p=2)')
# plt.plot(RK_RK4_x, RK_RK4_y, 'o-', label='Runke-Kutta (p=4)')
# plt.plot(x_an, model_analytisch(x_an), '-', label='Analytisch')
# plt.ylabel('y-Achse')
# plt.xlabel('x-Achse')
# plt.legend()
# plt.grid(True)
# plt.show()

# ----------------------------------------- Praktikum 10 -----------------------------------------
# implizite Mittelpunktsregel (nicht für Systeme)
# 0.5 | 0.5
#     | 1.0
def implizit_mittelpunkt(x0, X, f, df, N, tol=1e-3, A=0.5, B=1, C=0.5):
    # x0: Startpunkt
    # X: Endpunkt
    # f: Funktion
    # df: Jacobi Matrix der Funktion f
    # N: Anzahl der Schritte
    # tol: Toleranz
    max_iter=20

    h= (X-x0[0])/(N-1)
    print( "IMP: Schrittweite: ", h)
    x = np.zeros(N)
    y = np.zeros(N)
    x[0] = x0[0]
    y[0] = x0[1]
    s = 0

    #Anfangspunkt berechnen Für die Berechung des ersten K wertes werden die Randwerte in die Funktion gross F eingesetzt
    k = f(x[0],y[0])

    for i in range(1,N):
        print("K: ", k)
        step = 0
        
        x[i] = x[i-1] + h
        r = k - f(x[i-1] + h, y[i-1] + (h*k))

        while np.abs(r) > tol and step < max_iter:
            j = -df(x[i-1] + h*C, y[i-1] + (h*k*A))
            delta_k = -r / (1 - (h*j))
            k = k + delta_k
            r = k - f(x[i-1] + h*C, y[i-1] + (h*k*A))
            step += 1

        y[i] = y[i-1] + (h*k*B)
    return x, y


# ========== Example: implizite Mittelpunktsregel ==========
# def model1(x, y):
#     return -4*y
# 
# def model1_df(x, y):
#     return 4
# 
# def model1_analytisch(x):
#     return np.exp(-4 * x)
# 
# x0 = [0, 1]
# x_imp_mittelpunkt,y_imp_mittelpunkt = implizit_mittelpunkt(x0, 1, model1, model1_df, 10, 1e-3)
# 
# x_an = np.linspace(0,1,100)
# plt.figure(1)
# plt.plot(x_imp_mittelpunkt, y_imp_mittelpunkt, 'o-', label='Implizites Mittelpunktregel')
# plt.plot(x_an, model1_analytisch(x_an), '-', label='y(x) Analytisch')
# plt.ylabel('y-Achse')
# plt.xlabel('x-Achse')
# plt.title('Implizites Mittelpunktregel')
# plt.legend()
# plt.grid(True)
# plt.show()


# Implizite Trapezregel (nicht für Systeme)
# 0 | 0
# 1 | 0.5   0.5
#   | 0.5   0.5
def implizit_Trapez(x0, X, f, df, N, tol=1e-3):
    # x0: Startpunkt
    # X: Endpunkt
    # f: Funktion
    # df: Jacobi Matrix der Funktion f
    # N: Anzahl der Schritte
    # tol: Toleranz
    max_iter=20
    # c | a
    #   | b
    A = [[0.0, 0.0], [0.5, 0.5]]
    B = [0.5, 0.5]
    C = [0.0, 1.0]

    h = (X-x0[0])/(N-1)
    print( "IMP: Schrittweite: ", h)
    x = np.zeros(N)
    y = np.zeros(N)
    x[0] = x0[0]
    y[0] = x0[1]
    s = 0

    for i in range(1,N):
        #print("K: ", k)
        step = 0
        
        x[i] = x[i-1] + h
        k_1 = f(x[i-1], y[i-1])
        k_2 = k_1
        r = k_2 - f(x[i-1]+h, y[i-1] + h*(k_1*A[s][0] + k_2*A[s][1]))
        
        while np.abs(r) > tol and step < max_iter:
            j = -df(x[i-1] + h*C[1], y[i-1] + h*(k_1*A[s][0] + k_2*A[s][1]))
            delta_k = -r / (1-h*j)
            k_2 += delta_k
            r = k_2 - f(x[i-1] + h*C[1], y[i-1] + h*(k_1*A[s][0] + k_2*A[s][1]))
            step += 1
        
        #print("k's: ", kl[0], kl[1])
        y[i] = y[i-1] + h*(k_1*B[0] + k_2*B[1])
    return x, y


# ========== Example: implizite Trapezregel ==========
# def model2(x,y):
#     return (-x**2)/(y)
#     
# def model2_df(x,y):
#     return (x**2)/(y**2)
#     
# def model2_analytisch(x):
#     return -np.sqrt(2/3) * np.sqrt(24 - x**3)
# x0 = (0, -4)
# x_imp_dirk2, y_imp_dirk2 = implizit_Trapez(x0, 1, model2, model2_df, 100, 1e-3)
# 
# 
# x_an = np.linspace(0,1,100)
# plt.figure(2)
# plt.plot(x_imp_dirk2, y_imp_dirk2, 'o-', label='Implizites Trapez')
# plt.plot(x_an, model2_analytisch(x_an), '-', label='y(x) Analytisch')
# plt.ylabel('y-Achse')
# plt.xlabel('x-Achse')
# plt.legend()
# plt.grid(True)
# plt.show()
# ----------------------------------------- Praktikum 11 -----------------------------------------
#  Runge-Kutta-Verfahren (RK4) (für ein autonomes System)
def runga_kutta_RK4_2Systeme(x0, t0, h, tend, f):
    # x0: Anfangswerte [x1, x2, ...]
    # t0: Anfangswert
    # h: Schrittweite
    # tend: Zeit bis zum Ende
    # f: Funktion mit n Gleichungen

    N = int((tend - t0) / h)
    t = np.zeros(N+1)
    x = np.zeros((N+1,np.shape(x0)[0]))
    
    t[0] = t0
    x[0,:] = x0
    for i in range(1, N+1):
        t[i] = t[i-1] + h

        k1 = f(t[i-1], x[i-1,:])
        k2 = f(t[i-1] + 0.5*h, x[i-1,:] + 0.5*h*k1)
        k3 = f(t[i-1] + 0.5*h, x[i-1,:] + 0.5*h*k2)
        k4 = f(t[i-1] + h, x[i-1,:] + h*k3)
        x[i,:] = x[i-1,:] + h*(k1 + 2*k2 + 2*k3 + k4)/6
    return x, t


# ========== Example: Runge-Kutta-Verfahren (RK4) (für ein autonomes System) ==========
# Parameter
# p       = 1.184  # kg/m^3
# Ca      = 0.45  # constant
# Rkugel  = 0.1   # m
# m       = 0.05  # kg
# g       = 9.81  # m/s^2
# ks      = 0.05  # N/m
# v0      = 0.0     # m/s
# k       = 0.5 * Ca * p * np.pi * Rkugel**2 
# print(k)
# # Hier Modelfunktion definieren
# def model1(t, x):
#     x0 = x[0] # dx
#     x1 = x[1] # x
#     return np.array([x1,  -(k/m)*np.sign(x1+v0)* (x1+v0)**2 - ((ks*x0)/m) + g])
# 
# # Hier Anfangswerte, Endzeitpunktun und schrittweite definieren
# x0 = [0, 0]
# t0 = 0
# tend = 100
# h = 0.25
# 
# v0 = 0.0
# xdata, t = runga_kutta_RK4_2Systeme(x0, t0, h, tend, model1)
# 
# plt.figure(1)
# plt.plot(t, xdata,  '.-', label='v0=0')
# plt.ylabel('dx[m]')
# plt.xlabel('t')
# plt.title('RK4 mit Schrittweite h=0.25')
# plt.legend()
# plt.grid(True)
# plt.show()

# Implizite Mittelpunktsregel (für 2 Systeme)
# C | A  =  1/2 | 1/2
#   | B  =      | 1
def implizite_mittelpunkt_2Systeme(x0, h, tend, f, df, tol):
    # x0: Startpunkt
    # X: Endpunkt
    # f: Funktion
    # df: Jacobi Matrix der Funktion f
    # N: Anzahl der Schritte
    # tol: Toleranz
    max_iter=20

    N = int(tend / h)+1
    t = np.zeros(N)
    x = np.zeros((N,np.shape(x0)[0]))
   
    x[0,:] = x0
    k = f(t[0], x[0,:])

    for i in range(1,N):
        step = 0
        t[i] = t[i-1] + h
        r = k - f(x[i-1], x[i-1] + (h*0.5*k))
        J = -df(t[i-1] + h*0.5, x[i-1] + (h*k*0.5))
        # einheitsmatrix:
        I2 = np.eye(2)
        M = I2 - h * 0.5 * J

        while (np.abs(r[0]) > tol or np.abs(r[1]) > tol) and step < max_iter :
            # Solve linear system with LR
            b = -r
            ATA = np.dot(M.T, M)
            ATb = np.dot(M.T, b)
            delta_k = np.linalg.solve(ATA, ATb)
            k = k + delta_k
            r = k - f(t[i-1] + h*0.5, x[i-1] + (h*k*0.5))
            step += 1

        x[i,:] = x[i-1,:] + h*k
       
    return x, t

# ========== Example: Runge-Kutta-Verfahren (RK4) (für ein autonomes System) ==========
# Modelfunktion definieren
# def model1(t, x):
#     x1 = x[0] # x
#     x0 = x[1] # dx
#     return np.array([x0,  -(k/m)*(x0+v0)*np.abs(x0+v0) - ((ks*x1)/m) + g])
# 
# # Jacobi-Matrix der Funktion f
# def model1_df(t, x):
#     x0 = x[0] # dx
#     x1 = x[1] # x
#     return np.array([[0, 1], 
#                     [-ks/m, -(k/m)*2*np.abs(x0)]])
# 
# # Anfangswerte, Endzeitpunkt und Schrittweite einstellen
# x0 = (0, 0)
# tend = 100
# h = 0.25
# 
# v0 = 0.0
# v0_x, t = implizite_mittelpunkt_2Systeme(x0, h, tend, model1, model1_df, 1e-3)
# 
# plt.figure(1)
# plt.plot(t, v0_x[:,0],  '.-', label='v0=0')
# plt.plot(t, v0_x[:,1],  '.-', label='v0=0')
# plt.ylabel('x[m]')
# plt.xlabel('t')
# plt.title('Implizite Mittelpunktregel mit Schrittweite h=0.25')
# plt.legend()
# plt.grid(True)
# plt.show()


# ----------------------------------------- Praktikum 12 -----------------------------------------
# Partielle Differentialgleichungen

# Schritt 1: Allgemeine Daten definieren
# Parameter festlegen
Tl = 100         # °C
Tr = 0           # °C
L = 0.1          # m
a = 3.8 * 10**-6 # m^2/s

N = 30                      # Anzahl der Schritte
Nmiddle = (int)(N/2+1)      # Mitte zwischen den Platten bestimmen
h = L/(N+1)                 # Schrittweite bestimmen
dt = 1                      # delta t bestimmen

# Schritt 2: Matrix A erstellen (auf ZF schauen wie die Matrix aussieht)
A = np.zeros((N,N))
for i in range(0,N-1):
    A[i][i] = -2
    A[i+1][i] = 1
    A[i][i+1] =1
A[N-1][N-1] = -2
A = a/h**2 * A


# Schritt 3: Vektor v erstellen (auf ZF schauen wie)
v = np.zeros(N)
v[0] = a*Tl/h**2

# Schritt 4: Zeitvektor erstellen
K = 1000
#t = np.zeros(K)
x = np.linspace(0, 1, N+2)
#norm = np.zeros(K)

# Schritt 5: Anfangsbedingungen festlegen (Vektor u in ZF schauen)
u = np.zeros((N+2,K))
u[:Nmiddle, 0] = Tl
u[(Nmiddle+1):, 0] = Tr
u[0,:] = 100
u[-1,:] = 0
print(u[:,0])

# Schritt 6: System mit FTCS-Methode lösen
# 0 | 1
#   | 1

for k in range(1, K):
    u[1:N+1, k] = (np.identity(N)+dt*A)@u[1:N+1,k-1] + dt*v
    #norm[k] = np.linalg.norm(u[:, k] - u[:, k-1])
    #t[k] = k*dt

plt.figure(1)
plt.plot(x, u[:,0], '.-', label='t = 0s')
plt.plot(x, u[:,10-1], '.-', label='t = 10s')
plt.plot(x, u[:,100-1], '.-', label='t = 100s')
plt.plot(x, u[:,K-1], '.-', label='t = 1000s')
plt.xlabel("x [m]")
plt.ylabel("T [°]")
plt.legend()
plt.grid()

#plt.figure(2)
#plt.loglog(t[1:], norm[1:], '.')
#plt.axhline(0.1)
#plt.grid()
#plt.show()



# Alternative mit Crank-Nicolson Berechnen

LU,piv = lu_factor(np.identity(N)-dt/2*A)

for k in range(1, K):
    u[1:N+1, k] = lu_solve((LU,piv),(np.identity(N)+dt/2*A)@u[1:N+1,k-1] + dt*v)
    #norm[k] = np.linalg.norm(u[:, k] - u[:, k-1])
    #t[k] = k*dt

#plt.figure()
#plt.plot(x, u[:,0], '.-', label='t = 0s')
#plt.plot(x, u[:,10-1], '.-', label='t = 100s')
#plt.plot(x, u[:,K-1], '.-', label='t = 1000s')
#plt.xlabel("x [m]")
#plt.ylabel("T [°]")
#plt.legend()
#plt.grid()
#plt.show()




print("Load done: num_gDgln.py")








