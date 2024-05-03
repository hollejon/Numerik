import numpy as np
import matplotlib.pyplot as plt



f =     10**5       # Frequenz in HZ
R1 =    800.0         # Primärwiderstand in Ohm
Lp =    0.00005   # Primärinduktivität in H    
Lps =   150.0*10**-6  # Streuinduktivität in H
Lsp =   150.0*10**-6  # Streuinduktivität in H
R2 =    6.0           # Sekundärwiderstand in Ohm
Ls =    500*10**-6  # Sekundärinduktivität in H
U0 =    4.0           # Leerlaufspannung in V

print("f:", f, " R1:", R1, " Lp:", Lp, " Lps:", Lps, " Lsp:", Lsp, " R2:", R2, " Ls:", Ls, " U0:", U0)

# Anfangsbedingungen




# DGL System Ax =b
A = np.array([[(R1*1/Lp)/(1-Lps*Lsp*1/Ls),(R2*1/Ls*Lps)/(1-Lps*Lsp*1/Ls)],
              [(R2*(1/Ls)*Lsp*Lps)/(1+U0*(1/Lp)*Lsp*Lps),(R2*(1/Ls)*Lsp)/(1+U0*(1/Lp)*Lsp*Lps*(1/Ls))]])

b = np.array([(U0*(1/Ls))/(1-Lps*Lsp*(1/Ls)),0])


def runga_kutta_RK4_2Systeme(x0, tend, h, f):
    # x0: Anfangswerte
    # y0: Anfangswerte
    # tend: Zeit bis zum Ende
    # h: Schrittweite
    # f: Funktion mit n Gleichungen
    N = int(tend / h) +1
    # Vektor mit dx-daten und x-daten
    x = np.array([np.zeros((N)), np.zeros((N))])
    # Zeitvektor
    t = np.zeros((N))
    
    x[0][0] = x0[0]
    x[1][0] = x0[1]

    for i in range(1, N):
        t[i] = t[i-1] + h

        k1 = f([t[i], x[0][i-1],               x[1][i-1]])
        k2 = f([t[i], x[0][i-1] + 0.5*h*k1[1], x[1][i-1] + 0.5*h*k1[0]])
        k3 = f([t[i], x[0][i-1] + 0.5*h*k2[1], x[1][i-1] + 0.5*h*k2[0]])
        k4 = f([t[i], x[0][i-1] + h*k3[1],     x[1][i-1] + h*k3[0]])

        x[0][i] = x[0][i-1] + h*((1/6)*k1[1] + (1/3)*k2[1] + (1/3)*k3[1] + (1/6)*k4[1]) #dx
        x[1][i] = x[1][i-1] + h*((1/6)*k1[0] + (1/3)*k2[0] + (1/3)*k3[0] + (1/6)*k4[0]) #x

    return x, t


def RK4(t0,x0,h,tend,f):
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
    return t, x

#Versuch 1
def model1(t, x):
    #t = x[0]
    #x1 = x[1] # I1
    #x0 = x[2] # I2
    t = t
    x0 = x[0]
    x1 = x[1]

    U0 = 4*np.sin(2*np.pi*f*t)
    return np.array([(Lsp*Lps*(1/Ls)*R1*x0*(1/Lp) + R2*x1*(1/Ls)*Lsp) / (1+U0*(1/Lp)*Lsp*Lps*(1/Ls)),
                      (U0*(1/Ls) + R2*x1*(1/Ls)*Lps - R1*x0*(1/Lp))/(1-Lps*Lsp*(1/Ls))])

#Versuch 2
def model2(x):
    t = x[0]
    x0 = x[1] # I1
    x1 = x[2] # I2

    #t = t
    #x1 = x[0]
    #x0 = x[1]


    U0 = 4*np.sin(2*np.pi*f*t)
    return np.array([(R2*x0*Lps*Ls - R1*x1*(1/Lp)-U0) / (1 - (Lsp*Lps*Ls)),
                      ((R1*x1)-U0)/((Lp/Lsp)-Lps)])

x0 = (0.0, 0.0)
tend = 0.00005
h =    0.000000001


I_x, t = runga_kutta_RK4_2Systeme(x0, tend, h, model2)
#t2, I_x2 = RK4(0,x0,h,tend, model2)



# =================== Plot erstellen ===================
plt.figure(2)
plt.plot(t, I_x[0],  '.-', label='v0=0')
#plt.plot(t2, I_x2[:,0],  '.-', label='I1')
plt.ylabel('I1[A]')
plt.xlabel('t')
plt.title('RK4 mit Schrittweite h=0.1ns')
plt.legend()
plt.grid(True)


plt.figure(3)
plt.plot(t, I_x[1],  '.-', label='v0=0')
#plt.plot(t2, I_x2[:,1],  '.-', label='v0=0')
plt.ylabel('I2[A]')
plt.xlabel('t')
plt.title('RK4 mit Schrittweite h=0.1ns')
plt.legend()
plt.grid(True)
plt.show()
