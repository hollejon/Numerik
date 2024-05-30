"""
 _   _                           _           _   _____                 _   _                 
| \ | |_   _ _ __ ___   ___ _ __(_) ___ __ _| | |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
|  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ _` | | | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
| |\  | |_| | | | | | |  __/ |  | | (_| (_| | | |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
|_| \_|\__,_|_| |_| |_|\___|_|  |_|\___\__,_|_| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
"""
"""
LUT(M)..........................................LU decomposition for tridiagonal matrix
LU(A............................................LU decomposition for general matrix A
fbSubsT(LU, b)..................................Solve Ax = b for tridiagonal matrix LU
linsolve(A, b)..................................Solve Ax = b for general matrix A
fbSubs(LR, b)...................................Solve Ax = b for general matrix A
mysign(x).......................................Signum Funktion
e(n)............................................Einheitsvektor
HouseholderTransformation(w)....................Householder transformation
Newton(x0, F, df, tol, K).......................Newton Iteration
rndOrtho(n).....................................random orthogonal matrix
gauss_newton(function, jacobi, x0, xdata, ydata, tol, max_iter,
              damped, maxDampingIter)...........Gauss-Newton Algorithmus für nicht lineare Gleichungen
"""


print("Load Library: numeric_functions.py")
import matplotlib.pyplot as plt
from matplotlib.pylab import norm
import numpy as np
import scipy as sp
from scipy.linalg import solve_triangular,lu,lu_solve
import sympy as sym


"""
Uebersicht:

-------------------- Praktiukum 2 --------------------
LU-Zerlegung
def LU(M):
    M: Quadratische MAtrix   

Vorwaerts / Rueckwaertseinsetzen
def fbSubs(LR, b):
    LR: LR-Zerlegte Matrix
    b: rechte seite von Ax=b 

lineares Gleichungssystem A*x = b loesen.
def linsolve(A, b):
    A: MAtrix
    b: rechte seite von A*x = b


-------------------- Praktiukum 3 --------------------
# LU-Zerlegung für tridiagonale Matrix
# in: tridiagonale Matrix M
def LUT(M):
    M: tridiagonale Matrix M     

# Solve Ax = b for tridiagonal matrix LU
def fbSubsT(LU, b):
    LU: LU zerlegte tridiagonal Matrix
    b: rechte seite von A*x=b

    
-------------------- Praktiukum 4 --------------------
cholesky Zerlegung Beispiel: Z 239


-------------------- Praktiukum 5 --------------------
QR-Zerlegung Beispiel: Z: 270


-------------------- Praktiukum 6 --------------------
# NichtLineare Gleichungssysteme und Ausgleichsrechnung
# Newton-Verfahren

def Newton(x0, F, df, tol = 1e-5, K = 1000):
    x0: Anfangswerte
    F: Modellfunktion
    df: Jacobifunktion
    tol: toleranz
    K: Anzahl Iterationen

-------------------- Praktiukum 7 --------------------
# NichtLineare Gleichungssysteme und Ausgleichsrechnung
# Gauss-Newton-Verfahren

def gauss_newton(function, jacobi, x0, xdata, ydata, tol=1e-6, max_iter=100, damped=False, maxDampingIter=100):
    function:   Modellfunktion
    jacobi:     Jacobi-Matrix
    x0:         Anfangswerte
    xdata:      Vektor mit Zeitdaten
    ydata:      Vektor mit Messdaten
    tol:        Toleranz
    max_iter:   Maximale iterationen



def rndOrtho(n):                                #Beispielcode
def rndCond(n, cond):                           #Beispielcode

"""


# ----------------------------------------- Praktikum 1 -----------------------------------------
# Differenzenquotienten, log-log-Plots für Fehler gegen Schrittweite
def f(x):
    return np.cos(x)
x0 = 1

# exakter Wert der Ableitung
y = -np.sin(1)

DeltaRechts = []
DeltaLinks = []
DeltaZentral = []
Hs = 10.**np.arange(-20,-1) # logarithmische Schrittweite

for h in Hs:
    # Fehler des rechtsseitigen Differenzenquotient
    DeltaRechts.append(np.abs(y-(f(x0+h)-f(x0))/h))
    
    # Fehler des linksseitigen Differenzenquotient
    DeltaLinks.append(np.abs(y-(f(x0)-f(x0-h))/h))
    
    # Fehler des zentralen Differenzenquotient
    DeltaZentral.append(np.abs(y-(f(x0+h)-f(x0-h))/(2*h)))


#plt.loglog(Hs,DeltaRechts,'o-',label='Vorwaertsdifferenzenquotient')
# die beiden folgenden Zeilen können Sie nach dem Implementieren
# der weiteren Differenzenquotienten aktivieren
# plt.loglog(Hs,DeltaLinks,'.-',label='Rueckwaertsdifferenzenquotient')
# plt.loglog(Hs,DeltaZentral,'.-',label='zentraler DiffQuotient')
# plt.xlabel('Schrittweite')
# plt.ylabel('absoluterFehler')
# plt.title('Fehlerentwicklung der Differenzenquotienten fuer h->0')
# plt.legend()
# plt.grid()
# plt.show()


# ----------------------------------------- Praktikum 2 -----------------------------------------
# LR-Zerlegung ohne und mit Spaltenpivotisierung, Vorwaerts- und Rückwaertseinsetzen

# LR-Zerlegung der quadratischen Matrix A    
# in: quadratische Matrix A
#out: 
# - A wird überschrieben, unteres Dreieck = L (ohne Diagonale), oberes Dreieck = R
# - idx: Indexvektor der Zeilenvertauschungen
def LU(A):
    m = A.shape[0]
    idx = np.array(range(m))    
    for k in range(0,m):
        for i in range(k+1,m):
            if(A[k][k] == 0):
                #Pivot ist gleich 0 -> Vertausche Zeilen
                # Aktuelle zeile Speichern
                row_act = A[k].copy()
                
                # Zeile überschreiben
                A[k] = A[k+1]
                A[k+1] = row_act

                # P Vektor anpassen
                idx[k] = k+1
                idx[k+1] = k

            A[i][k]= A[i][k]/A[k][k]
 
            for j in range(k+1,m):
                A[i][j] = A[i][j] - A[i][k]*A[k][j]

    return A, idx


# Vorwawrts- und Rückwaertseinsetzen

# Solve Ax = b for general matrix A
def fbSubs(LR, b):
    n = len(b)   
    # Vorwaertseinsetzen
    y = np.zeros(n)
    y[0]=b[0]
    for i in range(1,n):
        y[i] = b[i] - np.dot(LR[i][:i], y[:i])
   
    # Rueckwaertseinsetzen
    x = np.zeros(n)
    x[n-1]= y[n-1]/LR[n-1][n-1]
    for j in range(n-2, -1, -1):  
        x[j] = (y[j] - np.dot(LR[j][j+1:], x[j+1:]))/LR[j][j]
 
    return x

# lineares Gleichungssystem A*x = b loesen.
def linsolve(A, b):
    M_LU, idx = LU(A)
    #Fuer L*R*x = P*b muessen wir P*b berechnen
    rows_b = len(b)
    P_b = np.zeros(rows_b)
    
    for i, val in enumerate(idx):
        # i := 0,1,2,...
        # val := index der Zeile von P
        P_b[i] = b[val]

    res = fbSubs(M_LU, P_b)
    return res

# ========== Example: Solve linear system with LR ==========
 
#b=y_m
#ATA = np.dot(A.T, A)
#ATb = np.dot(A.T, b)
#alphalr = np.linalg.solve(ATA, ATb)

# Berechnete alphas auf einen Zeitvektor anwenden

# t = np.linspace(0, 10, 100)
# for i in range(0, len(t)):
#     y_mLR[i] = alphalr[0]*A1(t[i]) + alphalr[1] + alphalr[2]*t[i] + alphalr[3]*t[i]**2 + alphalr[4]*t[i]**3


# plt.plot(t,y_mLR,'-', label='Approximation mit LR')
# plt.xlabel('Frequenz')
# plt.ylabel('Signal')
# plt.legend()
# plt.grid()
# plt.show()


# ----------------------------------------- Praktikum 3 -----------------------------------------
# auch für Tridiagonalmatrizen, Anwendung auf lineare Randwertprobleme

# LU-Zerlegung für tridiagonale Matrix
# in: tridiagonale Matrix M
def LUT(M): 
    n = M.shape[1] 
    for k in range(1,n): 
        for i in range(0,2): 
            if i == 0: 
                M[i][k] = M[i][k]/M[i+1][k-1] 
            else: 
                M[i][k] = M[i][k]-M[i-1][k]*M[i+1][k-1]
    return M


# Vorwärts- und Rückwärtseinsetzen für tridiagonale Matrix
#in: LU (output from LUT), vector b
#out: vector x s.t. L@U@x == b
  
# Solve Ax = b for tridiagonal matrix LU
def fbSubsT(LU, b):
    n = len(b)          #calculate length of b
    y = np.zeros(n)     #array with size of n
    x = np.zeros(n)     #array with size of n
    
    #vorwaerts
    y[0] = b[0]         #The first element of y is set equal to the first element of vector b
    for i in range(1,n):
        y[i] = b[i] - LU[0, i] * y[i-1]   #For each element of y (starting from the second one), it subtracts the dot prodt of the corresponding row of LU with the elements of y from b[i]
                    #This line computes the last element of x by dividing the last element of y by the last diagonal element of LU    
    #rueckwaerts
    x[-1]= y[-1] / LU[1,-1]
    for i in range(n-2,-1,-1):
        x[i] = (y[i] - LU[2,i] * x[i+1]) / LU[1,i]  #This loop computes the remaining elements of x iteratively starting from the second last one. It subtracts the dot product of the corresponding row of LU with the elements of x from the corresponding element of y and then divides it by the diagonal element of LU
    
    return x

# Lösen des linearen Systems
#y = solve_triangular(L, A @ y_m, lower=True)
#c = solve_triangular(L.T, y, lower=False)

# ----------------------------------------- Praktikum 4 -----------------------------------------
# Lineare Ausgleichsrechnung
# Cholesky zerlegung
y_m = [5, 6, 7, 8, 9]
A = np.array([[1, 2, 3, 6, 7], [4, 5, 6, 7, 7], [7, 8, 9, 7, 7], [10, 11, 12, 8, 8], [13, 14, 15, 8, 8]])

b=y_m
ATA = np.dot(A.T, A)
ATb = np.dot(A.T, b)

#    L = np.linalg.cholesky(ATA)
#    y = np.linalg.solve(L, ATb)
#    alphacy = np.linalg.solve(L.T, y)


# Berechnete alphas auf einen Zeitvektor anwenden
# for i in range(0, len(t)):
#     y_mcy[i] = alphacy[0]*A1(t[i]) + alphacy[1] + alphacy[2]*t[i] + alphacy[3]*t[i]**2 + alphacy[4]*t[i]**3
# 
# 
# plt.plot(t,y_mcy,'-', label='Approximation mit Cholesky')
# plt.xlabel('Frequenz')
# plt.ylabel('Signal')
# plt.legend()
# plt.grid()
# plt.show()


# ----------------------------------------- Praktikum 5 -----------------------------------------
# Lineare Ausgleichsrechnung
# QR-zerlegung
Q,R = np.linalg.qr(A)

# Solve liner system with QR
y_num =  np.dot(Q.T, y_m)
alpha = solve_triangular(R, y_num, trans=0, lower=False, unit_diagonal=False, overwrite_b=False, check_finite=True)

# Berechnete alphas auf einen Zeitvektor anwenden
# for i in range(0, len(t)):
#     y_m1[i] = alpha[0]*A1(t[i]) + alpha[1] + alpha[2]*t[i] + alpha[3]*t[i]**2 + alpha[4]*t[i]**3
# 
# 
# plt.plot(t,y_m1,'-', label='Approximation mit QR')
# plt.xlabel('Frequenz')
# plt.ylabel('Signal')
# plt.legend()
# plt.grid()
# plt.show()


# ----------------------------------------- Praktikum 6 -----------------------------------------
# NichtLineare Gleichungssysteme und Ausgleichsrechnung
# Newton-Verfahren
def Newton(x0, F, df, tol = 1e-5, K = 1000):
    k=0
    x = x0
    r = 1
    res_k = np.array([])

    while(r > tol and k < K):
        k = k+1
        J = df(x)
        y = F(x)

        # Gleichungssystem lösen
        ATA = np.dot(J.T,J)
        ATB = np.dot(J.T,y) 

        dx = np.linalg.solve(ATA,ATB)
        #print("Resultat: ", x, dx)
        x = x - dx
        r = np.linalg.norm(F(x))  
        res_k = np.append(res_k,[r],axis= 0)

    return x, res_k


# ========== Example: Newton-Verfahren ==========
# p_l1 = 2.0 # Länge des 1. Roboter Arms
# p_l2 = 1.0 # Länge des 2. Roboter Arms
# 
# P = [-1.0, 2.0]
# 
# # Nullstellenform
# def modell(t):
#     p = np.array([p_l1*np.cos(t[0]) + p_l2 *np.cos(t[0]+t[1]) - P[0], 
#                   p_l1*np.sin(t[0]) + p_l2 *np.sin(t[0]+t[1]) - P[1]])
#     return p
# 
# #Jacobi-Matrix berechnen
# def modell_df(x):
#     J = np.array([[-p_l1*np.sin(x[0])-p_l2*np.sin(x[0]+x[1]),-p_l2*np.sin(x[0]+x[1])],
#                   [ p_l1*np.cos(x[0])+p_l2*np.cos(x[0]+x[1]), p_l2*np.cos(x[0]+x[1])]])
#     return J
# 
# x0 = np.array([1, 1])
# a_k, res_k= Newton(x0, modell, modell_df)
# 
# plt.plot(res_k, 'o', label='Messwerte')
# plt.xlabel('Step')
# plt.ylabel('Abweichung')
# plt.legend()
# plt.grid()
# plt.yscale("log")
# plt.show()


# ----------------------------------------- Praktikum 7 -----------------------------------------
# NichtLineare Gleichungssysteme und Ausgleichsrechnung
# Gauss-Newton-Verfahren
def gauss_newton(function, jacobi, x0, xdata, ydata, tol=1e-6, max_iter=100, damped=False, maxDampingIter=100):
    """
    Gauss-Newton Algorithmus für nicht lineare Gleichungen

    Eingabe
    function:   Funktion welche gefitet werden soll
    jacobi:     Jakobi-Matrix der Funktion welche gefittet werden soll: Jacobi = Df(x)
    x0:         Erster Schätzwert der Parameter. Müssen in der nähe der Lösung liegen damit der Algorithmus konvergiert.
    ydata:      Vektor mit Messdaten
    xdata:      Vektor mit Zeitdaten
    max_iter:   Maximale Interationsschritte (default: 100)
    tol:        Toleranz der Kovnergenz (default: 1e-6)

    Rückgabe:
    parameter:  Lösungsvektor, Zahlen der Parameter, damit die Lösung den kleinsten fehler Quadrate der Messdaten entspricht
    """

    k = 0 #Anzahl Interationen
    parameter = x0.copy() #muss in der nähe der Lösung liegen damit es konvergiert
    y = function(xdata, parameter) - ydata
    j = jacobi(xdata, parameter)

    r = np.array([np.linalg.norm(j.T@y)])
    print("k", "x1\t", "x2\t", "x3\t", "x4\t", "r\t\t", "s")
    
    while r > tol and k < max_iter:
        j = jacobi(xdata, parameter)
        y = function(xdata, parameter) - ydata
       
        # Normalengleichung lösen
        #deltax = np.linalg.solve(j.T@j, -j.T@y)

        # Mit QR-Zerlegung
        Q,R = np.linalg.qr(j)
        deltax = solve_triangular(R, -Q.T@y, lower=False)
        # Falls Dämpfung aktiviert ist
        if damped == True:
            delta_k = 1
            n = 1
            while  np.linalg.norm(parameter + deltax*delta_k) > np.linalg.norm(parameter) and n < maxDampingIter:
                print("Dämpfung", n)
                delta_k = delta_k / 2
                n = n + 1
            deltax = deltax*delta_k

        parameter = parameter + deltax
        k += 1

        r = np.linalg.norm(j.T@y)
        s = np.sum((y-function(xdata, parameter))**2)
        print(k, round(parameter[0], 3), 
               round(parameter[1], 3), "\t",
               round(parameter[2], 3), "\t",
               round(parameter[3], 3),"\t",
               "{:.3e}\t".format(r),
               "{:.3e}".format(s))

    return parameter




# ============ Beispiel für Gauss-Newton Verfahren: ============

# Parameter for Gauss-Newton test
# tk = np.array([0.1, 0.3, 0.7, 1.2, 1.6, 2.2, 2.7, 3.1, 3.5, 3.9])
# yk = np.array([0.558, 0.569, 0.176, -0.207, -0.133, 0.132, 0.055, -0.090, -0.069, 0.027])
# 
# # Modellfunktion definieren
# def model(t, x):
#     return x[0]*np.exp(-x[1]*t)*np.sin(x[2]*t+x[3])
# 
# # Jakobi-Matrix der Modellfunktion definieren
# def model_df(t, x):
#     print(t)
#     J =  np.array([        np.exp(-x[1]*t)*np.sin(x[2]*t+x[3]),
#                    -x[0]*t*np.exp(-x[1]*t)*np.sin(x[2]*t+x[3]),
#                     x[0]*t*np.exp(-x[1]*t)*np.cos(x[2]*t+x[3]),
#                      x[0]*np.exp(-x[1]*t)*np.cos(x[2]*t+x[3])])
#     return J.T
# 
# 
# # Startwerte
# startwerte = np.array([0.5, 1, 3, 1]) # = Anzahl Parameter in der Jacobi-Funktion
# alpha = gauss_newton(model, model_df, startwerte, tk, yk, tol=1e-6, max_iter=100)
# print(alpha)
# 
# # Grössere Anzahl an Messpunkten
# time = np.linspace(0,4,501)
# 
# # Messpunkte und gefittete Funktion plotten
# plt.plot(tk, yk, 'o', label='Messpunkte')
# plt.plot(time, model(time, alpha),'-', label='Gauss-Newton')
# plt.ylabel('Signal')
# plt.xlabel('Zeit')
# plt.legend()
# plt.grid(True)
# plt.show()



# Erzeugen von speziellen Matrizen

# random orthogonal matrix
def rndOrtho(n):
    S = np.random.rand(n,n)
    S = S - S.T
    O = np.linalg.solve(S - np.identity(n), S + np.identity(n))
    return O


# random matrix with specified condition number 
# n =: Zeilenzahl 
# cond =: geforderte Kondition
def rndCond(n, cond):    
    d = np.logspace(-np.log10(cond)/2, np.log10(cond)/2,n);
    A = np.diag(d)
    U,V = rndOrtho(n), rndOrtho(n)
    return U@A@V.T



print("Load Library done: numeric_functions.py")



"""
https://docs.scipy.org/doc/scipy/reference/integrate.html
def DGLs(y_vec,t):
    y_1, y_2 = y_vec
    dy_1_dt = 3*y_1 + 2*y_2 - (2*t**2 + 1)*np.exp(2*t)
    dy_2_dt = 4*y_1 + y_2 + (t**2 +2*t - 4)*np.exp(2*t)
    return np.array([dy_1_dt,dy_2_dt])


def dgl2(y_vec,t):
    y_1, y_2 = y_vec
    dy_1_dt = (t**2)/y_1
    dy_2_dt = 0
    return np.array([dy_1_dt,dy_2_dt])

N=100
tval = np.linspace(0, 1, N+1)
y_1_0 = 0
y_2_0 = -4
initialval = np.array([y_1_0,y_2_0])
Loes = sp.integrate.odeint(dgl2, initialval, tval)



plt.xlabel(r"$\rm t$")
plt.ylabel(r"$\rm y_1,y_2$")
plt.plot(tval, Loes[:, 0],c="blue", label=r"$\rm y_1(t)$")
plt.plot(tval, Loes[:, 1],c="red", label=r"$\rm y_2(t)$")
plt.legend(loc='upper center',fontsize=16)
plt.grid(True)
plt.show()

"""