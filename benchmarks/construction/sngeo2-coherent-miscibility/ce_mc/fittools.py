import numpy as np
def ols(Phi,y): return np.linalg.lstsq(Phi,y,rcond=None)[0]
def loocv(Phi,y):
    H=Phi@np.linalg.pinv(Phi); r=y-H@y; return np.sqrt(np.mean((r/(1-np.diag(H)))**2))
def lasso(Phi,y,lam,iters=5000,tol=1e-10):
    # unpenalized intercept column 0 and point column 1? penalize all but col0
    n,p=Phi.shape; w=np.zeros(p); norm=(Phi**2).sum(0)/n; r=y-Phi@w
    for it in range(iters):
        wmax=0
        for j in range(p):
            rho=Phi[:,j]@r/n+norm[j]*w[j]
            if j==0: new=rho/norm[j]
            else: new=np.sign(rho)*max(abs(rho)-lam,0)/norm[j]
            d=new-w[j]
            if d!=0: r-=Phi[:,j]*d; w[j]=new; wmax=max(wmax,abs(d))
        if wmax<tol: break
    return w
def kfold(Phi,y,fitfun,k=5,seed=0):
    rng=np.random.default_rng(seed); idx=rng.permutation(len(y)); folds=np.array_split(idx,k); err=[]
    for f in folds:
        tr=np.setdiff1d(idx,f); w=fitfun(Phi[tr],y[tr]); err+=list(y[f]-Phi[f]@w)
    return np.sqrt(np.mean(np.square(err)))
def random_corr(orbs,x):
    return np.array([1.0]+[(2*x-1)**o['n'] for o in orbs])
