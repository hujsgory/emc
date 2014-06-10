import numpy
import scipy.linalg as la

def block_lu(a, nc):
    la.inv(a[:nc,:nc],overwrite_a=True)
    a[nc:,:nc]  = numpy.dot(a[nc:,:nc], a[:nc,:nc])
    a[nc:,nc:] -= numpy.dot(a[nc:,:nc], a[:nc,nc:])
    la.inv(a[nc:,nc:], overwrite_a=True)

    
def block_solve(a, b, nc):
    x = b.copy()
    x[nc:]-= numpy.dot(a[nc:,:nc],x[:nc])
    x[nc:] = numpy.dot(a[nc:,nc:],x[nc:])
    x[:nc] = numpy.dot(a[:nc,:nc],x[:nc]-numpy.dot(a[:nc,nc:],x[nc:]))
    return x


def iterative(A, b, M, max_iter = 50, tol = 1e-30, block=False, nc=None):
    return bicgstab(A, b, M, max_iter, tol, block, nc)

    ## \function bicgstab
    #  \brief Stabilized bi-conjugate gradient method with preconditioning (BiCGStab)
    #  \param M Matrix object with factorized matrixes S
    #  \param b vector of right-hand members
def bicgstab(A, b, M, max_iter, tol, block, nc):
    n_cond = b.shape[1]
    n = A.shape[0]
    X = numpy.ones((n, n_cond))
    V = numpy.zeros((n, n_cond))
    P = numpy.zeros((n, n_cond))
    R = b - numpy.dot(A, X)
    Rt = R.copy()
    S = numpy.zeros((n, n_cond))
    T = numpy.zeros((n, n_cond))
    norm_R0 = la.norm(R)
    alpha = numpy.ones(n_cond)
    beta = numpy.zeros(n_cond)
    rho = numpy.zeros(n_cond)
    rho_old = numpy.ones(n_cond)
    omega = numpy.ones(n_cond)
    for niter in xrange(max_iter):
        for i in xrange(n_cond):
            rho[i] = numpy.dot(Rt[:,i], R[:,i])
            if rho[i] == 0.0:
                raise ValueError
            beta[i] = (rho[i]/rho_old[i])*(alpha[i]/omega[i])
        P = R + beta*(P - omega*V)
        if block:
            Pt = block_solve(M, P, nc)
        else:
            Pt = la.lu_solve(M, P)
        V = numpy.dot(A, Pt)
        for i in xrange(n_cond):
            alpha[i] = rho[i]/numpy.dot(Rt[:,i], V[:,i])
            S[:,i] = R[:,i] - alpha[i]*V[:,i]
            X[:,i] += alpha[i]*Pt[:,i]
        if la.norm(S)/norm_R0 <= tol:
            break
        if block:
            St = block_solve(M, S, nc)
        else:
            St = la.lu_solve(M, S)
        T = numpy.dot(A, St)
        for i in xrange(n_cond):
            omega[i] = numpy.dot(T[:,i], S[:,i])/numpy.dot(T[:,i], T[:,i])
            X[:,i] += omega[i]*St[:,i]
            R[:,i] = S[:,i] - omega[i]*T[:,i]
        if la.norm(R)/norm_R0 <= tol:
             break
        rho_old[:] = rho[:]
    return X


def bicg(A, b, M, max_iter, tol, block, nc):
    n_cond = b.shape[1]
    n = A.shape[0]
    Mt = (numpy.transpose(M[0]), M[1])
    At = numpy.transpose(A)
    X = numpy.ones((n, n_cond))
    R = b - numpy.dot(A, X)
    Rt = R.copy()
    P = numpy.zeros((n, n_cond))
    Pt = numpy.zeros((n, n_cond))
    rho = numpy.zeros(n_cond)
    rho_old = numpy.ones(n_cond)
    norm_R0 = la.norm(R)
    alpha = numpy.zeros(n_cond)
    beta = numpy.zeros(n_cond)
    for niter in xrange(max_iter):
        if block:
            Z = block_solve(M, R, nc)
            Zt = block_solve(Mt, R, nc)
        else:
            Z = la.lu_solve(M, R)
            Zt = la.lu_solve(Mt, R)
        for i in xrange(n_cond):
            rho[i] = numpy.dot(Zt[:,i], Rt[:,i])
            if rho[i] == 0.0:
                raise ValueError
            beta[i] = rho[i]/rho_old[i]
            P[:,i] = Z[:,i] + beta[i]*P[:,i]
            Pt[:,i] = Zt[:,i] + beta[i]*Pt[:,i]
        Q = numpy.dot(A, P)
        Qt = numpy.dot(At, Pt)
        for i in xrange(n_cond):
           alpha[i] = rho[i]/numpy.dot(Pt[:,i], Q[:,i])
           X[:,i] += alpha[i]*P[:,i]
           R[:,i] -= alpha[i]*Q[:,i]
           Rt[:,i] -= alpha[i]*Qt[:,i]
        res = la.norm(R)/norm_R0
        print res
        if res <= tol:
             break
        rho_old[:] = rho[:]
    return X

def tfqmr(A, b, M, max_iter, tol, block, nc):
    pass


def qmr(A, b, M, max_iter, tol, block, nc):
    n_cond = b.shape[1]
    n = A.shape[0]
    X = numpy.ones((n, n_cond))
    R = b - numpy.dot(A, X)
    norm_R0 = la.norm(R)
    for i in xrange(n_cond):
        xi[i] = rho[i] = norm_R0
    if block:
        Z = block_solve(M, R, nc)
    else:
        Z = la.lu_solve(M, R)
    for niter in xrange(max_iter):
        for i in xrange(n_cond):
            if rho[i]==0 or xi[i]==0:
                raise ValueError
            V[:,i]  = Vt[:,i]/rho[i]
            Y[:,i] /= rho[i]
            W[:,i]  = Wt[:,i]/xi[i]
            Z[:,i] /= xi[i]
            delta[i] = numpy.dot(Z[:,i],Y[:,i])
        if block:
            Yt = block_solve(M, Y, nc)
            Zt = block_solve(M, Z, nc)
        else:
            Yt = la.lu_solve(M, Y)
            Zt = la.lu_solve(M, Z)  
          
    return X

def cgs(A, b, M, max_iter, tol, block, nc):
    n_cond = b.shape[1]
    n = A.shape[0]
    X = numpy.ones((n, n_cond))
    Q = numpy.zeros((n, n_cond))
    R = b - numpy.dot(A, X)
    Rt = R.copy()
    U = numpy.zeros((n, n_cond))
    P = numpy.zeros((n, n_cond))
    norm_R0 = la.norm(R)
    alpha = numpy.ones(n_cond)
    beta = numpy.zeros(n_cond)
    rho = numpy.zeros(n_cond)
    rho_old = numpy.ones(n_cond)
    for niter in xrange(max_iter):
        for i in xrange(n_cond):
            rho[i] = numpy.dot(Rt[:,i], R[:,i])
            if rho[i] == 0.0:
                raise ValueError
            beta[i] = rho[i]/rho_old[i]
            U[:,i] = R[:,i] + beta[i]*Q[:,i]
            P[:,i] = U[:,i] + beta[i]*(Q[:,i] + beta[i]*P[:,i])
        if block:
            Pt = block_solve(M, P, nc)
        else:
            Pt = la.lu_solve(M, P)
        Vt = numpy.dot(A, Pt)
        for i in xrange(n_cond):
            alpha[i] = rho[i]/numpy.dot(Rt[:,i], Vt[:,i])
            Q[:,i] = U[:,i] - alpha[i]*Vt[:,i]
        if block:
            Ut = block_solve(M, U+Q, nc)
        else:
            Ut = la.lu_solve(M, U+Q)
        Qt = numpy.dot(A, Ut)
        for i in xrange(n_cond):
            X[:,i] += alpha[i]*Ut[:,i]
            R[:,i] -= alpha[i]*Qt[:,i]
        res = la.norm(R)/norm_R0
        #print res
        if res <= tol:
             break
        rho_old = rho
    return X

def gmres(A, b, M, max_iter, tol, block, nc):
    pass
    
def fom(A, b, M, max_iter, tol, block, nc):
    pass