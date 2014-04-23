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


    ## \function iterative
    #  \brief Stabilized bi-conjugate gradient method with preconditioning (BiCGStab)
    #  \param M Matrix object with factorized matrixes S
    #  \param b vector of right-hand members
def bicgstab(A, b, M, max_iter = 50, tol = 1e-30, block=False, nc=None):
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
                    break
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
            norm_R = la.norm(R)/norm_R0
            if norm_R <= tol:
                 break
            rho_old = rho
        return X