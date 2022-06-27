"""
            EE2703 Applied Programming Lab 
            Assignment 3: Fitting Data To Models
            NAME: SHAILESH PUPALWAR
            ROLLL NO.: EE20B100
"""

# Importing required libraries
import scipy.special as sp 
from pylab import * 
import numpy as np
from scipy.linalg import lstsq #module used to return least-squares solution to the matrix equation

# Read the "fitting.dat" file to extract the data from it.


try:
    data = np.loadtxt("fitting.dat")       #Loading the text into the variable "data"
except IOError:
    print('Error No file found.Please keep the "fitting.dat" file in same folder') #Return a error message
    exit()

deviation = np.logspace(-1, -3, 9) # 9 numbers between 10^(-3) and 10^(-1) 

rn, cln = data.shape # Finding the dimensions of the data, i.e., no. of rows and columns using "shape"(a tuple)

t = np.linspace(0, 10, rn) # rn elements equally found between 0 and 10


# Defining our function excluding the noise
def g(t, A, B): 
    return A*sp.jv(2, t) + B*t # jv(n, t) --> Bessel Function

import matplotlib.pyplot as plt 

# Code for Q-4
for i in range(cln - 1):
    # Define legend()
    plot(data[:, 0], data[:, i + 1], label = '$\sigma$' + "=" + str(np.around(deviation[i], 3))) 
# x-axis = time(t) , y-axis = g(t, A, B) i.e., function value 
plt.legend()
# Plot the above legend on our graph

# Labelling the x and y axis
xlabel(r'$t$', size = 15)
ylabel(r'$f(t)+n$', size = 15)
# Titling the graph
title(r'Q4: Data to be fitted to theory')
grid(True)
# Plotting the original graph with true value
plot(t, g(t, 1.05, -0.105), label = r"True value") 
plt.legend()
# Drawing the graph
show() 

# Code for Q-5
errorbar(t[::5], data[::5, 1], deviation[1], fmt='ro', label = r"Error bar") 
xlabel(r"$t$", size = 15)
title(r"Q5: Data points for $\sigma$ = 0.1 along with exact function") 
plt.legend() 
plot(t, g(t, 1.05, -0.105), label = r"True value") # for y-axis
plt.legend() 
show()

# Define matrix M , N 
M = np.zeros((rn, 2)) # Order of matrix rn x 2 and is initialised
N = np.zeros((rn, 1)) # Initialising the matrix N, here N = MP , where P = Transpose(1.05  -0.105)  
for i in range(rn):
    M[i, 0] = sp.jv(2, data[i, 0]) 
    M[i, 1] = data[i, 0] 
    N[i] = M[i, 0]*1.05 + M[i, 1]*(-0.105) 


# Define matrix A, B
A = linspace(0, 2, 20) # A = 0, 0.1, ....., 2 
B = linspace(-0.2, 0, 20) # B = -0.2, -0.19, ....., 0

fk = g(data[:, 0], 1.05, -0.105) # fk = g(t, A0, B0) is a vector

# Epsilon calculation
epsilon = np.zeros((len(A), len(B))) #epsilon --> "mean squared error"

# Using dual for loops , as epsilon[i][j] = (1/101)* sigma(k = 0 to 101) {fk - g(tk, Ai, Bj)}
for i in range(len(A)):
    for j in range((len(B))):
        epsilon[i,j] = np.mean(np.square(fk - g(t, A[i], B[j]))) 


# Code for Q-8
cntrplt = plt.contour(A, B, epsilon, 15) 
plot(1.05, -0.105, "ro")
annotate(r"$Exact\ location$", xy=(1.05, -0.105))
plt.clabel(cntrplt, inline = True)
plt.xlabel(r"$A$", size = 15)
plt.ylabel(r"$B$", size = 15)
plt.title(r"Q8: Countour plot for $\epsilon_{ij}$")
show() # Plotting the contour plots      

# Calculation of errors(A, B) 
pred = [] 
Aerror = []
Berror = []
y_true = g(t, 1.05, -0.105) # True graph 
for i in range(cln - 1):
    p, resid, rank, sig = lstsq(M, data[:, i + 1])
    aerr = np.square(p[0] - 1.05) 
    ber = np.square(p[1] + 0.105)   
    # Updating the errors
    Aerror.append(aerr) 
    Berror.append(ber) 


# Code for Q10
plot(deviation, Aerror, "ro", linestyle = "--", linewidth = 1, label = r"$Aerr$") # Legend for Aerr in graph
plt.legend()
plot(deviation, Berror, "go", linestyle = "--", linewidth = 1, label = r"Berr") # Legend for Berr in graph
plt.legend()
grid(True)
plt.xlabel(r"Noise standard deviation") 
plt.ylabel(r"$MS Error$", size = 15)
plt.title("$Q10:Variation\ of\  error\  with\  noise$")
show()


# Code for Q11
plt.loglog(deviation, Aerror, "ro")
plt.errorbar(deviation, Aerror, np.std(Aerror), fmt="ro", label=r"$Aerr$") # Legend for Aerr in graph
plt.legend()
plt.legend()
plt.loglog(deviation, Berror,"go")
plt.errorbar(deviation, Berror, np.std(Berror), fmt="go", label=r"$Berr$") # Legend for Berr in graph
plt.legend()
grid(True)
plt.ylabel(r"$MS Error$", size = 15)
plt.title(r"$Q11: Variation\ of\ error\ with\ noise$")
plt.xlabel(r"$\sigma_{n}$", size = 15)
show()

# Function to compare two matrices
def matrixequal(P, Q):
    count = 0 
    
    for i in range(0, rn):
        if P[i] == Q[i]: 
            count += 1 

    if count == rn: # If all elements matched
        return True 
    else: # Else, return false
        return False

# Printing the output whether N is equal to fk or not, where N = MP and fk = g(t, A0, B0)
if matrixequal(N, fk):
    print("Both the matrices(N = MP and fk) are equal.")
else: 
    print("Both the matrices(N = MP and fk) are not equal.")            
