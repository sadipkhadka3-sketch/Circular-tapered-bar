import numpy as np

# INPUT DATA
P  = 10000.0      # N
L  = 1000.0       # mm
d1 = 50.0         # mm
d2 = 25.0         # mm
E  = 200000.0     # N/mm^2

# Exact displacement
u_exact = 4*P*L/(np.pi*E*d1*d2)

print("Exact displacement =",u_exact)

def tapered_bar_FEM(n):

    le = L/n

    K = np.zeros((n+1,n+1))

    for e in range(n):

        xm = (e+0.5)*le

        de = d1 + ((d2-d1)/L)*xm

        Ae = np.pi*de**2/4

        ke = (E*Ae/le)*np.array([[1,-1],
                                 [-1,1]])

        K[e:e+2,e:e+2] += ke

    F = np.zeros(n+1)

    F[-1] = P

    Kred = K[1:,1:]
    Fred = F[1:]

    Ured = np.linalg.solve(Kred,Fred)

    U = np.zeros(n+1)
    U[1:] = Ured

    return U[-1]

for n in [1,2,3]:

    ufem = tapered_bar_FEM(n)

    error = abs((ufem-u_exact)/u_exact)*100

    print("n =",n)
    print("Displacement =",ufem)
    print("Error =",error,"%")
