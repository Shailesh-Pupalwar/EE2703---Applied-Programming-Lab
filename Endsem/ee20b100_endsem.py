

#             ENDSEM Applied Programming Lab (EE2703)
#             Name: SHAILESH PUPALWAR
#             Roll No.: EE20B100  
#             Finding the Antenna currents in a half-wave dipole antenna


#Question-1
# Defining all the given variables 

import pylab
import numpy as np

pi = np.pi 

N = 4  # Number of sections in each half section of the antenna
Im = 1.0  # Current injected into the antenna
len = 0.5  # Quarter wavelength
w_no = pi  # Wave number = (2*pi)/(lambda)
dz = len/N  # Spacing of current samples

z = np.linspace(-len,len,2*N+1) # Points at which we determine the currents.
u = np.delete(z,[0,N,2*N]) #2*(N-1) locations of unknown currents
a = 0.01 # Radius of wire
mu_0 = 4e-7*pi # Permeability of free space



#Question-2
# Defining matrix M

def M(): # Function to determine and return the matrix M.
  M = (np.identity(2*(N-1)))*(1/(2*pi*a))
  return M
M = M()



#Question-3
# Computing vectors R_z, R_u and matrices PB, P

Z = np.meshgrid(z,z)
z_i, z_j = Z[0], Z[1]

#R_z determines distances from source and observer where source is point of wire.
R_z = np.sqrt((z_i-z_j)**2 + np.ones([2*N+1,2*N+1],dtype=complex)*(a**2))

U = np.meshgrid(u, u)
u_i, u_j = U[0], U[1]

#R_u determines distances from source and observer where observer is point where we want the field.
R_u = np.sqrt((u_i-u_j)**2 + np.ones([2*N-2,2*N-2],dtype=complex)*(a**2))

# Distances with respect to z = 0
R_in = np.delete(R_z[:][N],[N*2,0,N])  # Removing the three elements (first, middle, last) 


# P is the contribution to the vector potential due to unknown currents
P = (mu_0/(4*pi))*(np.cos(w_no*R_u)-(np.sin(w_no*R_u))*1j)*(1/R_u)*(dz)

# PB is the contribution to the vector potential due to current "In"
PB = (mu_0/(4*pi))*(np.cos(w_no*R_in)-(np.sin(w_no*R_in))*1j)*(R_in)*(dz)


#Question-4
# Computing Q and QB.

# Matrix corresponding to unknown currents
Q = (a/mu_0)*(P)*((1j)*(w_no)+(1/R_u))*(1/R_u)

# Matrix corresponding to the boundary current
QB = (a/mu_0)*(PB)*((1j)*(w_no)+(1/R_in))*(1/R_in)



#Question-5
# Here we are computing Estimated currents and Assumed currents.

J = np.linalg.inv(M-Q)@QB*Im
#estimated currents
I_esm = np.concatenate(([0],J[:N-1],[Im],J[N-1:],[0]))
#assumed currents
I_asm = Im*np.sin(w_no*(len-abs(z)))

# Printing out all the vectors and matrices asked in the question

print("R_z")
print((R_z).round(2))

print("R_u")
print((R_u).round(2))

print("R_in")
print((R_in).round(2))

print("P*1e8")
print((P*1e8).round(2))

print("PB*1e8")
print((PB*1e8).round(2))

print("Q")
print((Q).round(2))

print("QB")
print((QB).round(2))

print("I_esm")
print((abs(I_esm)).round(2))

print("I_asm")
print((I_asm).round(2))

print("z")
print((z).round(2))

print("u")
print((u).round(2))



# Plotting the required graph for both cases.

pylab.figure()
pylab.plot(z,abs(I_esm),label = "Estimated Current")
pylab.plot(z,I_asm,label = "Assumed Current")
pylab.xlabel(r"z")
pylab.ylabel(r"I")
pylab.legend()
pylab.title("Antenna currents at N=100")
pylab.grid()
pylab.savefig('Figure_1.png')

pylab.show()