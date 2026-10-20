# cython: language_level=3, boundscheck=False, wraparound=False
"""Modular evaluation kernel; coefficients are certified separately over QQ."""
import numpy as np
cimport numpy as np

def base_rows(points, weights, roots, extra, pairings, degrees, unsigned long long prime):
    """Characters and integral scaled power sums, evaluated modulo prime."""
    cdef Py_ssize_t i,j,k,a,b,n=len(points),nw=len(weights),nr=len(roots),ne=len(extra)
    cdef unsigned long long v,t,p1,p2,p3,p4,a1,a2,s,e2,e3,e4,ev,q,z
    cdef unsigned long long inv2=pow(2,-1,int(prime)),inv6=pow(6,-1,int(prime)),inv24=pow(24,-1,int(prime))
    cdef np.ndarray[np.int64_t,ndim=3] P=np.empty((n,7,9),dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=2] X=np.empty((n,7),dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=2] Y=np.empty((n,len(degrees)),dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=1] V=np.empty(nw+nr,dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=2] W=np.asarray(list(weights)+list(roots),dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=2] E=np.asarray([w for w,m in extra],dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=1] M=np.asarray([m for w,m in extra],dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=2] D=np.asarray(pairings,dtype=np.int64)%int(prime)
    for i in range(n):
        for j in range(7):
            t=int(points[i][j])%int(prime)
            P[i,j,4]=1
            for k in range(1,5):
                P[i,j,4+k]=pow(int(t),int(k),int(prime))
                P[i,j,4-k]=pow(int(t),-int(k),int(prime))
        p1=p2=p3=p4=0
        a1=a2=7
        for j in range(nw+nr):
            v=1
            for k in range(7):
                v=(v*<unsigned long long>P[i,k,4+W[j,k]])%prime
            V[j]=v
            if j<nw:
                q=v*v%prime
                p1=(p1+v)%prime
                p2=(p2+q)%prime
                p3=(p3+q*v)%prime
                p4=(p4+q*q)%prime
            else:
                a1=(a1+v)%prime
                a2=(a2+v*v)%prime
        s=0
        for j in range(ne):
            v=<unsigned long long>M[j]
            for k in range(7):
                v=v*<unsigned long long>P[i,k,4+E[j,k]]%prime
            s=(s+v)%prime
        # Add multiples of prime before subtraction: unsigned arithmetic.
        e2=(p1*p1%prime+prime-p2)%prime*inv2%prime
        e3=(p1*p1%prime*p1%prime+3*prime-3*(p1*p2%prime)+2*p3)%prime*inv6%prime
        e4=(p1*p1%prime*p1%prime*p1%prime+12*prime-6*(p1*p1%prime*p2%prime)+3*(p2*p2%prime)+8*(p1*p3%prime)-6*p4)%prime*inv24%prime
        X[i,0]=a1
        X[i,1]=s
        X[i,2]=((a1*a1%prime+prime-a2)*inv2%prime+prime-a1)%prime
        X[i,3]=(e4+prime-e2)%prime
        X[i,4]=(e3+prime-p1)%prime
        X[i,5]=(e2+prime-1)%prime
        X[i,6]=p1
        for j in range(len(degrees)):
            Y[i,j]=0
        for j in range(nw):
            ev=0
            for k in range(D.shape[1]):
                a=k if D.shape[1]==nw else nw+k
                ev=(ev+<unsigned long long>D[j,k]*<unsigned long long>V[a])%prime
            for k in range(len(degrees)):
                z=pow(int(ev),int(degrees[k]),int(prime))
                Y[i,k]=(Y[i,k]+z)%prime
    return X,Y

def evaluation_rows(chars, parents, factors, unsigned long long prime):
    cdef np.ndarray[np.int64_t, ndim=2] X=np.asarray(chars,dtype=np.int64)
    cdef np.ndarray[np.int64_t, ndim=1] P=np.asarray(parents,dtype=np.int64)
    cdef np.ndarray[np.int64_t, ndim=2] F=np.asarray(factors,dtype=np.int64)
    cdef Py_ssize_t m=X.shape[0], n=P.shape[0],i,j,a,b,p
    cdef unsigned long long v
    cdef np.ndarray[np.int64_t, ndim=2] A=np.empty((m,n),dtype=np.int64)
    for i in range(m):
        for j in range(n):
            p=P[j]
            if p<0:
                a=F[j,0]
                A[i,j]=1 if a<0 else X[i,a]
            else:
                a=F[j,0]
                b=F[j,1]
                v=(<unsigned long long>A[i,p]*<unsigned long long>X[i,a])%prime
                if b>=0:
                    v=(v*<unsigned long long>X[i,b])%prime
                A[i,j]=v
    return A
