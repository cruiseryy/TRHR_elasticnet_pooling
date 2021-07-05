clear;clc;close all 
rng(0)

load('../sst.mat')
load('../TRHR_precip.mat')

precip(:,:,457:458) = [];
sprecp = mean(mean(precip));
TRHR = zeros(456/12,1);
for k = 5:10
    TRHR = TRHR + reshape(sprecp(k:12:456),456/12,1);
end

TRHR = (TRHR - mean(TRHR))./std(TRHR);

sst(isnan(sst)==1) = -9.9;
mm_sst = zeros(360,180,12);
std_sst = zeros(360,180,12);
for k = 1:12
    mm_sst(:,:,k) = mean(sst(:,:,k:12:12*39),3);
    std_sst(:,:,k) = std(sst(:,:,k:12:12*39),[],3);
end 
ssta = (sst - repmat(mm_sst,1,1,39))./repmat(std_sst,1,1,39);
sst_vector = zeros(468,2);
mark = 1;
Loc = zeros(2,1);
for i = 1:360
    for j = 1:180
        if ocean_mask(i,j) == 1
            sst_vector(:,mark) = reshape(ssta(i,j,:),468,1);
            Loc(:,mark) = [i; j];
            mark = mark + 1;
        end
    end
end
pause = 1;
% ssta(:,:,468) = [];
% ssta(:,:,1:11) = [];


pr_vector = sst_vector;

tmp = pr_vector';
R = tmp'*tmp;
% R = pr_vector*pr_vector';
[C, L] = eig(R);
% plot(diag(L)/trace(L))
eigenvector = pr_vector'*C;
% TTT = eigenvector'*eigenvector;
% coef = diag(TTT);
% coef = coef';
% coef = repmat(coef,size(eigenvector,1),1);
% eigenvector = eigenvector./sqrt(coef);
PC = C;
eigenvalue = diag(L)/trace(L);

score1 = zeros(13,1);

txx = C(:,end-29:end);
for lag = 0:12
    tmp_vector = txx(12+5-lag:12:468-lag,:);
    tmp_vector = [tmp_vector ones(38,1)];
    tfu = zeros(38,1);
    
    for kk = 1:38
        tmpxx = tmp_vector;
        tmpxx(kk,:) = [];
        tmpyy = TRHR;
        tmpyy(kk) = [];
%         tic
        b = regress(tmpyy, tmpxx);
%         toc
        tmpfit = tmp_vector*b;
        tfu(kk) = tmpfit(kk);
        pause = 1;

   
    end
    score1(lag+1) = sum(sign(tfu).*sign(TRHR) > 0)/38;
end

var = sum(eigenvalue(end-29:end))





