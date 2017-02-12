X = 2.*rand(11,2)-1;
Y = sign(sin(X(:,1))+X(:,2));

mex BRI.cpp -larmadillo -llapack -lblas
type = 'c';
L_fold = 5; 
kernel = 'RBF_kernel';

[gam,sig2] = tunelssvm({X,Y,type,[],[],kernel},'gridsearch','crossvalidatelssvm',{L_fold,'misclass'});
[alpha,b] = trainlssvm({X,Y,type,gam,sig2,kernel});
plotlssvm({X,Y,type,gam,sig2,kernel},{alpha,b});
