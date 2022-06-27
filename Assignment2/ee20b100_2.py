"""
EE2703
Applied Programming Lab 2022- ASSIGNMENT 2
Name : Shailesh Pupalwar
Roll No: EE20B100

"""

#Importing required modules
from sys import argv, exit
import math
import numpy as np

#counting the number of voltage sources   
def volcount(data):
    vol_count=0
    for line in data:
        list = line.split('#')[0].split()

        if(list[0][0] == "V"):
            vol_count += 1
    return vol_count

#counting the number of nodes in the circuit
def nodecount(data):
    cnt=0
    for line in data:
        list = line.split('#')[0].split()
        
        if(list[1] != 'GND'):
            if(int(list[1]) > cnt):
                cnt = int(list[1])
        if(list[2] != 'GND'):
            if(int(list[2]) > cnt):
                cnt = int(list[2])        
    return cnt   

#Defining the class resistor
class RLC:
    
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = value

    # General MNA entries for resistor,inductor and capacitor
    def MatrixEntries(self,name, n1, n2, value, G):
        if(n1 != 'GND' and n2 != 'GND'):
            x = int(n1)-1
            y = int(n2)-1
            z = value      
        
            G[x][x] += 1.0 / z
            G[x][y] += -1.0 / z
            G[y][x] += -1.0 / z
            G[y][y] += 1.0 / z
        
        if(n1 == 'GND'):
            y = int(n2)-1
            z = value
        
            G[y][y] += 1.0 / z

        if(n2 == 'GND'):
            x = int(n1)-1
            z = value
        
            G[x][x] += 1.0 / z 


#Defining class for independent voltage source
class volt_Src:
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = value
     

    #MNA entries for a independent voltage source
    def MatrixEntries(self,name, n1, n2, value, G,I,count,nodecnt):
        count += 1
        vkl = nodecnt + count -1
        if n1 != 'GND' and n2 != 'GND' :
            x = int(n1)-1
            y = int(n2)-1
            G[x][vkl] = 1
            G[y][vkl] = -1
            G[vkl][x] = 1
            G[vkl][y] = -1

        elif n1 == 'GND' :
            y = int(n2)-1
            G[vkl][y] = -1
            G[y][vkl] = -1
        else :
            x = int(n1)-1
            G[vkl][x] = -1
            G[x][vkl] = -1
        
        I[vkl][0] = value



#Declaring constant variables for better readability
CIRC = '.circuit'
END = '.end'
AC = '.ac'

#initialising a list called 'list' with 0 entries
list = [] 
count =0 # To count the index of voltage source


# checking if no of arguments are 2 or not
if len(argv) != 2:
    print('Your input arguments are not two,please make sure it should be two for valid command')
    exit()

"""
Making sure to drop an error message if wrong file was given as an input using try-catch
"""
try:
    with open(argv[1]) as f:
        data = f.readlines()
        f.close() # closing the file after reading
        start = -1; end = -2

        for line in data:              # extracting circuit definition start and end data
            if  CIRC== line[:len(CIRC)]:
                start = data.index(line)
               
            elif END== line[:len(END)]:
                end = data.index(line)

            elif AC == line[:len(AC)]:
                ac = data.index(line)  
                ac_words = line.split()
                #angular frequency of the ac circuit
                Ang_Freq=2*3.1415*float(ac_words[-1]) 

#validating circuit block, which is:'.circuit' should always ahead of '.end'
        if start >= end:
            print("Invalid circuit definition")
        else :
            for i in range(start+1, end):
                L = ((data[i].split('#')[0].split()))
                list.append(L)
               

        m = nodecount(data[start+1 : end]) #m is the number of nodes in circuit
        n = volcount(data[start+1 : end])  #n is the number of voltage sources in the circuit
        order = m+n # m+n gives the order of the G matrix
        print(order)
        
        G = np.zeros((order,order), dtype="complex")  # Conductance matrix 
        V = np.zeros((order,1), dtype="complex")  # Variable vector
        I = np.zeros((order,1), dtype="complex")  # Vector of independent sources
#for matrix entries we have to find name of object and its speciefic values for three classes which are resistor, capacitor and inductor
        for line in list:
            #  for Resistor 
            if(line[0][0] == 'R'):
                Value = float(line[3])
                object = RLC(line[0],line[1],line[2],Value)
                object.MatrixEntries(line[0],line[1],line[2],Value, G)
            # for Capacitor    
            if(line[0][0] == 'C'):
                Value = complex(0,-1/(Ang_Freq*float(line[3])))
                object = RLC(line[0],line[1],line[2],Value)
                object.MatrixEntries(line[0],line[1],line[2],Value, G)
            # for Inductor
            if(line[0][0] == 'L'):
                Value = complex(0,(Ang_Freq*float(line[3])))
                object = RLC(line[0],line[1],line[2],Value)
                object.MatrixEntries(line[0],line[1],line[2],Value, G)
            # for Independent Voltage source
            if(line[0][0] == 'V'):
                if(line[3] == 'ac'):
                    Value = complex((float(line[-2])*math.cos(float(line[-1])))/2,(float(line[-2])*math.sin(float(line[-1])))/2)
                    object = volt_Src(line[0],line[1],line[2],Value)
                    object.MatrixEntries(line[0],line[1],line[2],Value,G,I,count, m)
                else:
                    Value = float(line[3])
                    object = volt_Src(line[0],line[1],line[2],Value)
                    object.MatrixEntries(line[0],line[1],line[2],Value,G,I,count, m)

            
                    
        #Solving [G][V] = [I] for V matrix; V = (inverse of[G])[I]
        V = np.linalg.solve(G, I)
        print("Required Voltage matrix is \n")
        print(V)

        print('\n Required Voltages are :\n')
        for i in range(0,m+n):
            if ac == 0:
                if i<m:print("Voltage at node %d is :" %(i+1))
                if i>=m:print("Current through V%d source is :" %(i-m+1))
                print(float("%.2f"%V[i]))
            else:       # for dc sources the node wise voltage and current through the voltage source wil be printed by folloing code.
                if i<m:print("Voltage at node %d is :"%(i+1))
                if i>=m:print("Current through V%d source is"%(i-m+1))
                print(V[i])
except IOError:
    print('Invalid file')
    exit()    