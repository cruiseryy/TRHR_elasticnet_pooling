clear;clc;close all
rng(0)

load('../data/sst/sst_1979to2019.mat')
% load('../data/sst/lsmask.mat')
load('Indian_Pacific_mask.mat')
load('../data/precip/RS/trhr_precip_2021.mat')

% the wet season is defined as May to September (consistent with the wet
% season of TB from other studies and the growing season in Chen's work
% ~87% of the annual precip fall during the wet season based on 1981-2019
% (CHIRPS precip)
trhr_prcp = zeros(39,1);
for k = 1:39
    trhr_prcp(k) = sum(precip((k-1)*12+5:(k-1)*12+9));
end
% standardization of precipitation to remove the systematic shift around
% 2000.5
prcp2 = trhr_prcp;
trhr_prcp(1:20) = (prcp2(1:20) - mean(prcp2(1:20)))/std(prcp2(1:20));
trhr_prcp(21:end) = (prcp2(21:end) - mean(prcp2(21:end)))/std(prcp2(21:end));

% figure(1)
% hold on
% plot(1981:2019, trhr_prcp,'ro-')
% plot([1 1]*2000.5, 2.5*[-1 1], 'k--')
% grid
% ylabel('wet season precip anomaly')
% close(1)

wsz = 4; % non-overlapping pooling window size

% max pooling or average pooling is tested
ni = 360/wsz;
nj = 180/wsz;
nt = 492;
rsst = zeros(ni, nj, nt);
rmask = zeros(ni, nj);

for i = 1:ni
    for j = 1:nj
        if sum(sum(ocean_mask((i-1)*wsz+1:i*wsz, (j-1)*wsz+1:j*wsz))) == wsz^2
            
%             rsst(i,j,:) = reshape(mean(mean(sst((i-1)*wsz+1:i*wsz, (j-1)*wsz+1:j*wsz, :))),nt,1); % average pooling
            rsst(i,j,:) = reshape(max(max(sst((i-1)*wsz+1:i*wsz, (j-1)*wsz+1:j*wsz, :))),nt,1); % max pooling
%             pause = 1;
            rmask(i, j) = 1; % 'pooled' mask
        end
    end
end

% standardization
msst1 = zeros(ni, nj, 12);
msst2 = zeros(ni, nj, 12);
std_sst1 = zeros(ni, nj, 12);
std_sst2 = zeros(ni, nj, 12);
for k = 1:12
    msst1(:,:,k) = mean(rsst(:,:,k:12:12*22),3);
    std_sst1(:,:,k) = std(rsst(:,:,k:12:12*22),[],3);
    msst2(:,:,k) = mean(rsst(:,:,12*22+k:12:end),3);
    std_sst2(:,:,k) = std(rsst(:,:,12*22+k:12:end),[],3);
end 
ssta = zeros(ni, nj, nt);
ssta(:,:,1:12*22) = (rsst(:,:,1:12*22) - repmat(msst1,1,1,22))./repmat(std_sst1,1,1,22);
ssta(:,:,12*22+1:end) = (rsst(:,:,1+12*22:end) - repmat(msst2,1,1,19))./repmat(std_sst2,1,1,19);

sstv = zeros(492,2);
loc = zeros(2,1);
tmark = 1;

for i = 1:ni;
    for j = 1:nj
        if rmask(i,j) == 1
            sstv(:,tmark) = reshape(ssta(i,j,:),nt,1);
            loc(:,tmark) = [i;j];
            tmark = tmark + 1;
        end
    end
end

%----------------------------------------------

F = sstv'; 
F = detrend(F,0);
R = F'*F;
[C, L] = eig(R);
PC = F*C;

eigva = diag(L)/trace(L);
plot(eigva)

ll = 50;
sum(eigva(end-ll+1:end))
xx0 = C(:,end-ll+1:end);

cc = [];
for lag = 0:24
    
    xx = xx0(24+5-lag:12:end-lag,:);
    b = regress(trhr_prcp(1:20,:),[xx(1:20,:) ones(20,1) ]);
    tmpfit = xx*b(2:end) + b(1);
%     cc(lag+1) = corr(tmpfit(1:20), trhr_prcp(1:20));
    cc(lag+1) = corr(tmpfit(21:end), trhr_prcp(21:end));
%     lag
    
end

    
%----------------------------------------------
cc2 = cc;
load('elastic_net_cc_max_pooling.mat')
figure()
hold on
bar(0:24, [cc cc2'])
plot([-1 25],[1 1]*0.389,'k--')
plot([-1 25],-1*[1 1]*0.389,'k--')
legend('regularized','EOF OLS','location','northwest')
axis([-1 25 -0.4 0.75])
grid














