#!/usr/bin/env python3
"""

Name:           BSM_main.py 
Authors:        McAleavey,McCarthy,O'Brian 
Studnet ids:    15204643,14200512,98139495
Created:        6/11/2016    
Assignment:     Numerical Analytics and Software - Programming Assignment 2: 
                Linear equations  
            
Purpose:        Builds set of linear equations for solving the Black-scholes equation 
                and outputs the matrix to a file. It then converts the dense matrix to a TD 
                matrix and solves the system using the SOR emthod.

"""

###############################################################################

import os
import sys
sys.path.append('functions/')
#set woking dir 
#os.chdir('/Users/shanemccarthy/Dropbox/Masters/Semester 3/Numerical Analytics and Software/Assignment3/Code_final')

###############################################################################
# import librarys needed to execute code
#import matplotlib.pyplot as plt
import numpy 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import BuildSparse
import MIF

###############################################################################
# Black-Scholes Parameters
print('\n Building the Black-Scholes Matrix Ax=b and outputting to file\n');
print('See doc for usage');
print('Make sure to have cleared all variables before proceding by typing reset');

N = 80;                 # size of the stock price dimension (x-dim in ADE equation)
M = 100;                # size of the temporal array
T = 1.5 ;               # maturity date of option
R = 0.1;                # interest rate
sigma2 = 0.16;          # diffusivity
X = 1;                  # stock price t time T
Smax = 2;               # max stock price 
filename = 'nas_Sor.in.txt';   # filename of BSM matrix written to working directory
omega = 1.1
tol=1e-9
limit =100

'''
N = int(input("Please enter size of BSM Matrix: "));
M = int(input("Please enter number of time steps for solution: "));
T = float(input("Please enter maturity date of the option: "));
R = float(input("Please enter interest rate: "));
sigma2 = float(input("Please enter the: "));
X = float(input("Please enter number of time steps for solution: "));
Smax = float(input("Please enter number of time steps for solution: "));
filename = (input("Please enter number of time steps for solution: "));
'''

'''
N = int(sys.argv[1]);                 # size of the stock price dimension (x-dim in ADE equation)
M = int(sys.argv[2]);                 # size of the stock price dimension (x-dim in ADE equation)
T = float(sys.argv[3]);                 # size of the stock price dimension (x-dim in ADE equation)
R = float(sys.argv[4]);                 # size of the stock price dimension (x-dim in ADE equation)
sigma2 = float(sys.argv[5]);                 # size of the stock price dimension (x-dim in ADE equation)
X = float(sys.argv[6]);                 # size of the stock price dimension (x-dim in ADE equation)
Smax = float(sys.argv[7]);                 # size of the stock price dimension (x-dim in ADE equation)
filename = (sys.argv[8]);                 # size of the stock price dimension (x-dim in ADE equation)
omega = 1.1
tol=1e-9
limit =100
'''

omega = 1.1
tol=1e-9
limit =100

print('\n\t size of BSM Matrix =',N)
print('\t number of time steps for solution = ',M)
print('\t maturity date of the option = ',T)
print('\t interest rate = ',R)
print('\t diffusivity = ',sigma2)
print('\t stock price = ',X)
print('\t Max stock price = ',10*X)
print('\t filename = ',filename)

##############################################################################
##############################################################################
dS = Smax/N;            # price/spatial grid step 
dt = T/M;               # dt < 1/simga2/N/N
# set arrays needed to run BSM
BSM = numpy.zeros((N,N))
BSMTD = numpy.zeros((N,3))
#fprice = numpy.zeros((M,N));
B = numpy.zeros((N,1));

###############################################################################
# Populate Dense BSM Matrix

for i in range(0,N):
    for j in range(0,N):
        BSM[i][j]=0;
        if(i==j):
            BSM[i][j]=(1+dt*R+dt*sigma2*i*i);
        if(j==i+1):
            BSM[i][j]=(-0.5*i*dt*(i*sigma2+R));
        if(j==i-1):
            BSM[i][j]=(-0.5*i*dt*(i*sigma2-R));


###############################################################################
# Triadiagonal BSM Matrix
# convert full matrix to TD using function
BSMTD = BuildSparse.BTD(BSM);

#print('Matrix=\n',BSM);
#print('Matrix=\n',BSMTD);

###############################################################################
#  boundary conditions for European Put option

# Final condition
B[0] = max(X-0*dS,0) + 0.0*dt*(sigma2 - R)*X;
for i in range(1,N):
    B[i]=max(X-i*dS,0);

###############################################################################
# Writes out Matrix A and vectro b to file in the form 
# n 
# [A]
# [b]

save_path = '/output/'
wd=os.getcwd()
pathFilename= wd+save_path+filename

f = open(pathFilename, 'w');           
f.write("%d\n" %N);
for i in range(0,N):
    for j in range(0,N):
        f.write("%f " %BSM[i][j]); 
    f.write("\n" );
for i in range(0,N):
    f.write("%f " %B[i]); 
f.write("\n" );
f.close();
       
###############################################################################
###############################################################################
# option to test matrix here and
test_on = int(input("Type 1 to test or any other key to continue to solution: "));
if test_on==1: 
        locaction = (input("Please enter the path to the input file containing A and b: "));


###############################################################################
###############################################################################

fprice = numpy.zeros((M,N));

# Final condition
for i in range(0,N):
    fprice[0][i]=B[i];

# BC on S=0 and S=Smax
for i in range(0,M):
    fprice[i][0]=X;         # BC on S=0;
    fprice[i][N-1]=0;       # BC on S=Smax
 
# Main time loop for BSM Equation
for t in range(1,M):
    print(t)
    fprice[t][:] = MIF.SOR_Tridiagonal(N, omega, BSMTD, fprice[t-1][:], tol , limit);
    fprice[t][1] = fprice[t][1] + 0.5*dt*(sigma2 - R)*X;   
    fprice[t][0]=X;
    fprice[t][N-1]=0;

###############################################################################
# plotting the output
# generate 2 2d grids for the x & y bounds
TT, SS = numpy.mgrid[slice(0, T, dt),slice(0, Smax, dS)]

fig=plt.figure(figsize=(12,10))
ll=round(M/2);
plt.subplot(2,1,1);
plt.plot(SS[0][:],fprice[0][:],'r' );
plt.hold(True)
plt.plot(SS[ll][:],fprice[ll][:],'b' );
plt.plot(SS[M-1][:],fprice[M-1][:],'k' );
plt.grid(True)
plt.ylabel('Option Value ($)');
plt.xlabel('Stock Price ($)');
plt.title('BSM Solution')
    
plt.subplot(2,1,2);
plt.pcolor(SS, TT, fprice, cmap='rainbow_r');#, vmin=z_min, vmax=z_max)
plt.title('BSM Solution')
plt.ylabel('Time (yrs)');
plt.xlabel('Stock Price ($)');
#plt.axis([x.min(), x.max(), y.min(), y.max()])
plt.colorbar()

fig=plt.figure(figsize=(12,10))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(SS,TT,fprice, cmap='coolwarm',antialiased=True)
#ax.set_zlim(-1.01, 1.01)
#fig.colorbar(surf, shrink=0.5, aspect=5)  
ax.view_init(elev=30,azim=70) 
plt.ylabel('Time (yrs)');
plt.xlabel('Stock Price ($)');
#plt.zlabel('Option Value ($)');

###############################################################################




