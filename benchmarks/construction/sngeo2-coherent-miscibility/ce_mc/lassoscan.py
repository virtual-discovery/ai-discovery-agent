import numpy as np, fittools as F, pickle
xs,ef,Phi,orbs=pickle.load(open('big.pkl','rb'))
y=ef*1000
for lam in [0.01,0.03,0.1,0.3,1,2,3,5]:
    w=F.lasso(Phi,y,lam)
    cv=np.mean([F.kfold(Phi,y,lambda P,t:F.lasso(P,t,lam),seed=s) for s in range(5)])
    rm=np.sqrt(np.mean((y-Phi@w)**2))
    print(lam,'nz',int((abs(w)>1e-8).sum()),'rmse',round(rm,2),'cv',round(cv,2),'Hrand(0.5)',round(F.random_corr(orbs,0.5)@w,1),'H(0.25)',round(F.random_corr(orbs,0.25)@w,1),'H(.75)',round(F.random_corr(orbs,0.75)@w,1))
