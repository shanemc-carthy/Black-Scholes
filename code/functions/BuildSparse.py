"""

Gareth, Shane, Richard - NAS 2016 Assignment 3 27/10/2016

"""
import numpy #as np
###############################################################################
###############################################################################

def BS(A):

# scan the matrix for the size of the array
    Msize = 0;
    N=A.shape[0];
    for i in range(0,N):
        for j in range(0,N):
            if A[i][j]!=0:
                Msize = Msize + 1;           
            
#define the sparse matrix arrays with #elements
    SpM = numpy.zeros((Msize,3))

#loop over array and store sparse representation
    k = 0 ; r = 1; io=0;
    for i in range(0,N):
        for j in range(0,N):
            if A[i][j]!=0:
                SpM[k][0] = A[i][j];
                SpM[k][1] = j;
                k = k + 1;
                if io!=i:
                    SpM[r][2] = k-1;
                    r = r + 1;
                    io=i;
                    
# fill end of rowstart array              
    for i in range(1,Msize):
        if(SpM[i][2]==0):
            SpM[i][2] = Msize;

    return SpM
    
###############################################################################
###############################################################################
# build tridiagional matrix
def BTD(A):
    
    N=A.shape[0];
    TM = numpy.zeros((N,3))

    i=0; 
    TM[0][0] = 0;
    TM[0][1] = A[0][0];
    TM[0][2] = A[0][1];
    for i in range(1,N-1):
        TM[i][0] = A[i][i-1];
        TM[i][1] = A[i][i];
        TM[i][2] = A[i][i+1];
    i=N-1;
    TM[i][0] = A[N-1][N-2];
    TM[i][1] = A[N-1][N-1];
    TM[i][2] = 0.0;

    return TM
    #print('Tridiagonal matrix = \n',M3)

  
    