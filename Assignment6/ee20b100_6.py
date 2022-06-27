"""
            EE2703 Applied Programming Lab - 2022
            Assignment 5: The Laplace Transform 
            NAME: Shailesh Pupalwar
            ROLLL NO.: EE20B100
            
"""

# Importing the necessary modules.
from pylab import *            #Importing pylab
import scipy.signal as sp      

# The python code snippet for Q.1
ply11 = poly1d([1,0.5])
ply21 = polymul([1,1,2.5],[1,0,2.25]) 
X1 = sp.lti(ply11,ply21) # X(s)
t1,x1 = sp.impulse(X1,None,linspace(0,50,500)) # Calculate x(t) using impulse function

# The python code snippet for Q.2
# Repeating the above for 0.05 decay
ply12 = poly1d([1,0.05]) # numerator
ply22 = polymul([1,0.1,2.2525],[1,0,2.25]) # denominator
X2 = sp.lti(ply12,ply22)
t2,x2 = sp.impulse(X2,None,linspace(0,50,500)) # x(t) for deacy = 0.05

# The python code snippet for Q.3
# Calculate x(t) for different frequencies and plot them
H = sp.lti([1],[1,0,2.25])
for w in arange(1.4,1.6,0.05): # w = frequency and range is [1.4, 1.6]
	t = linspace(0,50,500)
	f = cos(w*t)*exp(-0.05*t)
	t,x,svec = sp.lsim(H,f,t) # Calculate y(t) from x(t) and H(s)

    # The plot of x(t) for various frequencies vs time.
	figure(2)
	plot(t,x,label='w = ' + str(w))
	title("x(t) for different frequencies")
	xlabel(r'$t\rightarrow$')
	ylabel(r'$x(t)\rightarrow$')
	legend(loc = 'upper left')
	grid(True)

# The python code snippet for Q.4: Coupled spring problem
t4 = linspace(0,20,500)
X4 = sp.lti([1,0,2],[1,0,3,0]) # X(s)
Y4 = sp.lti([2],[1,0,3,0])	# Y(s)
t4,x4 = sp.impulse(X4,None,t4) # Calculate x(t)
t4,y4 = sp.impulse(Y4,None,t4) # Calculate y(t)

# The python code snippet for Q.5: RLC circuit
temp = poly1d([1e-12,1e-4,1])
H5 = sp.lti([1],temp) # H(s)
w,S,phi = H5.bode() # Bode plot of transfer function

# The python code snippet for Q.6: Output for given input in RLC circuit
t6 = arange(0,25e-3,1e-7)
vi = cos(1e3*t6) - cos(1e6*t6) # Input function f(t)
t6,vo,svec = sp.lsim(H5,vi,t6) # Output voltage Vo(t)

# Plotting all the required plots that needs to be plotted.

# The plot x(t) vs t for Q.1
figure(0)
plot(t1,x1)
title("The solution x(t) for Q.1")
xlabel(r'$t\rightarrow$')
ylabel(r'$x(t)\rightarrow$')
grid(True)

# The plot of x(t) vs t for Q.2
figure(1)
plot(t2,x2)
title("The solution x(t) for Q.2")
xlabel(r'$t\rightarrow$')
ylabel(r'$x(t)\rightarrow$')
grid(True)

# The plot of x(t) and y(t) vs t for Q.4 
figure(3)
plot(t4,x4,label='x(t)')
plot(t4,y4,label='y(t)')
title("x(t) and y(t)")
xlabel(r'$t\rightarrow$')
ylabel(r'$functions\rightarrow$')
legend(loc = 'upper right')
grid(True)

# The magnitude bode plot for Q.5 
figure(4)
semilogx(w,S)
title("Magnitude Bode plot")
xlabel(r'$\omega\rightarrow$')
ylabel(r'$20\log|H(j\omega)|\rightarrow$')
grid(True)

# The phase bode plot for Q.5 
figure(5)
semilogx(w,phi)
title("Phase Bode plot")
xlabel(r'$\omega\rightarrow$')
ylabel(r'$\angle H(j\omega)\rightarrow$')
grid(True)

# The plot of Vo(t) vs t for large time interval.
figure(6)
plot(t6,vo)
title("The Output Voltage for large time interval")
xlabel(r'$t\rightarrow$')
ylabel(r'$V_o(t)\rightarrow$')
grid(True)

# The plot of Vo(t) vs t for small time interval.
figure(7)
plot(t6[0:300],vo[0:300])
title("The Output Voltage for small time interval")
xlabel(r'$t\rightarrow$')
ylabel(r'$V_o(t)\rightarrow$')
grid(True)

show()