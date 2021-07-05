clear;clc;close all 
rng(0) % for reproducibility

load ionosphere
Ybool = strcmp(Y,'g');
X = X(:,3:end);

xtr = X(1:200,:);
xte = X(201:end,:);
ytr = Ybool(1:200);
yte = Ybool(201:end,:);

[B,FitInfo] = lassoglm(xtr, ytr, 'binomial','NumLambda',25,'CV',10);

indx = FitInfo.Index1SE;
B0 = B(:,indx);
cnst = FitInfo.Intercept(indx);
B1 = [cnst;B0];
preds = glmval(B1,xte,'logit');
predlog = (preds>0.5);
acc1 = sum(predlog==yte)/length(yte)*100;


yp0 = 1 - Ybool;
yp1 = Ybool;

yptr0 = yp0(1:200); ypte0 = yp0(201:end);
yptr1 = yp1(1:200); ypte1 = yp1(201:end); 

[bp0, fp0] = lassoglm(xtr, yptr0, 'poisson', 'NumLambda',25,'CV',10);
[bp1, fp1] = lassoglm(xtr, yptr1, 'poisson', 'NumLambda',25,'CV',10);

k = xte*(bp1(:,fp1.IndexMinDeviance) - bp0(:,fp0.IndexMinDeviance)) + fp1.Intercept(fp1.IndexMinDeviance) - fp0.Intercept(fp0.IndexMinDeviance);

pr = exp(k)./(exp(k)+1);
prlog = (pr>0.5);
acc2 = sum(prlog==yte)/length(yte)*100;
figure()
corr(pr,preds)
scatter(preds, pr)
