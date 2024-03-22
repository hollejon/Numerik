
import numpy as np










"""

LU decomposition for tridiagonal matrix
in: a  =  [[0,      a_{21}, ..., a_{n-1,n-2}, a_{n,n-1}],
           [a_{11}, a_{22}, ..., a_{n-1,n-1}, a_{nn}],
           [a_{12}, a_{23}, ..., a_{n-1,n},   0]]

out: LU = [[0,      l_{21}, ..., l_{n-1,n-2}, l_{n,n-1}],
           [r_{11}, r_{22}, ..., r_{n-1,n-1}, r_{nn}],
           [r_{12}, r_{23}, ..., r_{n-1,n},   0]]
           
"""

def LUT(M): 
    n = M.shape[1] 

    for k in range(1,n): 
        for i in range(0,2): 
            if i == 0: 
                M[i][k] = M[i][k]/M[i+1][k-1] 
            else: 
                M[i][k] = M[i][k]-M[i-1][k]*M[i+1][k-1]
    return M



# LU decomposition for general matrix A
def LU(A):
    m = A.shape[0]
    idx = np.array(range(m))    
    #print("M", m, "idx: ", idx, "A: ", A)

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



"""
in: LU (output from LUT), vector b
out: vector x s.t. L@U@x == b
"""  
# Solve Ax = b for tridiagonal matrix LU
def fbSubsT(LU, b):

    n = len(b)          #calculate length of b
    y = np.zeros(n)     #array with size of n
    x = np.zeros(n)     #array with size of n
    
    #vorwärts
    y[0] = b[0]         #The first element of y is set equal to the first element of vector b
    for i in range(1,n):
        y[i] = b[i] - LU[0, i] * y[i-1]   #For each element of y (starting from the second one), it subtracts the dot prodt of the corresponding row of LU with the elements of y from b[i]
                    #This line computes the last element of x by dividing the last element of y by the last diagonal element of LU    
    #rückwärts
    x[-1]= y[-1] / LU[1,-1]
    for i in range(n-2,-1,-1):
        x[i] = (y[i] - LU[2,i] * x[i+1]) / LU[1,i]  #This loop computes the remaining elements of x iteratively starting from the second last one. It subtracts the dot product of the corresponding row of LU with the elements of x from the corresponding element of y and then divides it by the diagonal element of LU
    
    return x



# Solve Ax = b for general matrix A
def fbSubs(LR, b):
    # code
    n = len(b)   
   
     # Vorwärtseinsetzen
    y = np.zeros(n)
    y[0]=b[0]
    for i in range(1,n):
        y[i] = b[i] - np.dot(LR[i][:i], y[:i])
   
 
    # Rückwärtseinsetzen
    x = np.zeros(n)
    x[n-1]= y[n-1]/LR[n-1][n-1]
    for j in range(n-2, -1, -1):  
        x[j] = (y[j] - np.dot(LR[j][j+1:], x[j+1:]))/LR[j][j]
 
    return x



















# Erzeugen von speziellen Matrizen












# random orthogonal matrix
def rndOrtho(n):
    S = rnd.rand(n,n)
    S = S - S.T
    O = lin.solve(S - np.identity(n), S + np.identity(n))
    return O


# random matrix with specified condition number 
# n =: Zeilenzahl 
# cond =: geforderte Kondition
def rndCond(n, cond):    
    d = np.logspace(-np.log10(cond)/2, np.log10(cond)/2,n);
    A = np.diag(d)
    U,V = rndOrtho(n), rndOrtho(n)
    return U@A@V.T
