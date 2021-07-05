clear;clc;close all 

load('cc_realamp.mat')

figure()
hold on
bar(0:24,[cc_lasso cc_ols cc_eof cc_cca])
plot([-1 25],[1 1]*0.389,'k--')
plot([-1 25],-1*[1 1]*0.389,'k--')
legend({'regularized','ols','eof','cca'},'location','northwest')
axis([-1 25 -0.6 0.7])
grid
ylabel('Correlation Coef')
xlabel('Lead Time (mo)')
