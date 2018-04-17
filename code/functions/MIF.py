"""
Gareth, Shane, Richard - NAS 2016 Assignment 3 27/10/2016

"""
###############################################################################
# import librarys needed to execute code
import matplotlib.pyplot as plt
import time
import numpy #as np
import math
import decimal
import random

###############################################################################
# Gauss Seidel from lecture notes
# used in construction of the code but not in the solutions
def GaussSeidel(A, b, tol, limit):
    
    x = numpy.zeros_like(b)
    for iteration in range(limit):
        next = numpy.zeros_like(x)
        for i in range(A.shape[0]):
            s1 = numpy.dot(A[i, :i], next[:i])
            s2 = numpy.dot(A[i, i + 1:], x[i + 1:])
            next[i] = (b[i] - s1 - s2) / A[i, i]
        if numpy.allclose(x, next, rtol=tol): break
        x = next
    
    return x  

###############################################################################
# SOR for full Matrix
# benched against a matlab solution
def SOR_full(omega, A, b, tol, limit ):
   
    iterations = 0; 
    x = numpy.zeros_like(b);
    error = 1e3;

    while(iterations < limit and error > tol):
        error = 0;
        for i in range(A.shape[0]):
            sum1 = 0; 
            for j in range(A.shape[1]):
                sum1 = sum1 + A[i,j] * x[j];
            x[i] = x[i] + omega*(b[i]-sum1)/A[i,i]; 
            error = error +  abs(omega*(b[i]-sum1)/A[i,i]); 
        iterations=iterations+1;
        print ('Full SOR iterations=',iterations,' limit =',limit,'Error=',error)

    return x
      
###############################################################################
# SOR for Tridiagonal matrix inversion - using a priori knowledge of the
# finite-differnce stencil to calcualte solution with A3 := [N][3]
def SOR_Tridiagonal(N,omega, A3, b, tol , limit):
   
    iterations = 0; 
    x = numpy.zeros_like(b);
    error = 1e3;
    dia = 1; 
    #print ('dim',A3.shape[0])

    # while loop with tolerance and iteration thresholds
    while(iterations < limit and error > tol):
        error = 0; ## measure the sum of the error on iterations
        
        # for first row of the loop to avoid i-1 index address
        i=0;
        sum1 = A3[i,1]*x[i] + A3[i,2]*x[i+1];
        x[i] = x[i] + omega*(b[i]-sum1)/A3[i,dia]; 
        error = error +  abs(omega*(b[i]-sum1)/A3[i,dia]); 
                                    
        # full loop over all rows
        for i in range(1,N-1):
            sum1 = A3[i,0] * x[i-1] + A3[i,1] * x[i] + A3[i,2]*x[i+1];
            x[i] = x[i] + omega*(b[i]-sum1)/A3[i,dia]; 
            error = error +  abs(omega*(b[i]-sum1)/A3[i,dia]); 

        # for last row of the loop to avoid i+1 index address
        i=N-1;
        sum1 = A3[i,0]*x[i-1] + A3[i,1]*x[i];
        x[i] = x[i] + omega*(b[i]-sum1)/A3[i,dia]; 
        error = error +  abs(omega*(b[i]-sum1)/A3[i,dia]); 
                                    
        iterations=iterations+1;
        print ('TD SOR iterations=',iterations,' limit =',limit,'Error=',error)

    return x

###############################################################################
###############################################################################
###############################################################################
# SOR for Sparse  matrix inversion - using a priori knowledge of the
def SOR_Sparse(N,omega, AS, b, tol , limit):
   
    iterations = 0; 
    x = numpy.zeros_like(b);
    error = 1e3;
    #print ('Sparse=\n',AS)

    while(iterations < limit and error > tol):
        error = 0;
        for i in range(0,N):
            sum1 = 0 ;
            s1 = int(AS[i,2]); s2= int(AS[i+1,2]);
            for j in range(s1,s2):
                jj = int(AS[j,1]);
                sum1 = sum1 + AS[j,0] * x[jj];
                if jj == i:
                    dd = AS[j,0];
                #print ('i=',i,' j=',j,' jj=',jj,'s1=',s1,'s2=',s2,'val=',AS[j,0],'vec=',x[jj])
            x[i] = x[i] + omega*(b[i]-sum1)/dd; 
            error = error +  abs(omega*(b[i]-sum1)/dd); 
        iterations=iterations+1;
        print ('Sparse SOR iterations=',iterations,' limit =',limit,'Error=',error)

    return x
        
###############################################################################
###############################################################################
###############################################################################
