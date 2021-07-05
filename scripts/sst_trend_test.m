clear;clc;close all 

load('../data/sst/sst_1979to2019.mat')
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

figure()
subplot(2,1,1)
prcp3 = (prcp2 - mean(prcp2))/std(prcp2);
bar(1981:2019, prcp3)

subplot(2,1,2)
bar(1981:2019, trhr_prcp)

pause = 1;

wsz = 4;
ni = 360/wsz;
nj = 180/wsz;
nt = 492;
rsst = zeros(ni, nj, nt);
rmask = zeros(ni, nj);

nlon = 2.5:4:360;
nlat = -87.5:4:90;

for i = 1:ni
    for j = 1:nj
        if sum(sum(ocean_mask((i-1)*wsz+1:i*wsz, (j-1)*wsz+1:j*wsz))) == wsz^2
            
%             rsst(i,j,:) = reshape(min(min(sst((i-1)*wsz+1:i*wsz, (j-1)*wsz+1:j*wsz, :))),nt,1); % average pooling
            rsst(i,j,:) = reshape(max(max(sst((i-1)*wsz+1:i*wsz, (j-1)*wsz+1:j*wsz, :))),nt,1); % max pooling
%             pause = 1;
            rmask(i, j) = 1; % 'pooled' mask
        end
    end
end

% standardization
msst1 = zeros(ni, nj, 12);
msst2 = zeros(ni, nj, 12);
msst3 = zeros(ni, nj, 12);
std_sst1 = zeros(ni, nj, 12);
std_sst2 = zeros(ni, nj, 12);
std_sst3 = zeros(ni, nj, 12);
for k = 1:12
    msst1(:,:,k) = mean(rsst(:,:,k:12:12*22),3);
    std_sst1(:,:,k) = std(rsst(:,:,k:12:12*22),[],3);
    msst2(:,:,k) = mean(rsst(:,:,12*22+k:12:end),3);
    std_sst2(:,:,k) = std(rsst(:,:,12*22+k:12:end),[],3);
    msst3(:,:,k) = mean(rsst(:,:,k:12:end),3);
    std_sst3(:,:,k) = std(rsst(:,:,k:12:end),[],3); 
end 
ssta = zeros(ni, nj, nt);
ssta(:,:,1:12*22) = (rsst(:,:,1:12*22) - repmat(msst1,1,1,22))./repmat(std_sst1,1,1,22);
ssta(:,:,12*22+1:end) = (rsst(:,:,1+12*22:end) - repmat(msst2,1,1,19))./repmat(std_sst2,1,1,19);

ssta2 = (rsst - repmat(msst3,1,1,41))./repmat(std_sst3,1,1,41);

ssta_diff = mean(ssta2(:,:,12*22+1:end),3) - mean(ssta2(:,:,1:12*22),3);
ssta_diff(isnan(ssta_diff)==1) = 0;

figure()
test_pic(ssta_diff)
tlat = 0*ones(1,360);
tlon = 1:360; 
poslat = []; poslon = [];
neglat = []; neglon = [];
% p = 0.1 two-tailed (90% CI) = [-1.645*sigma, 1.645*sigma] = [-0.5152 0.5152] 
% assuming independent random variables for SSTA at each grid,
% Var(ssta_diff) = 1/22 + 1/19 
sig_threshold = 0.5152;

for i = 1:90
    for j = 1:45
        if ssta_diff(i,j) >= sig_threshold
            poslat = [poslat nlat(j)]; 
            poslon = [poslon nlon(i)];
        end
        if ssta_diff(i,j) <= -1*sig_threshold
            neglat = [neglat nlat(j)];
            neglon = [neglon nlon(i)];
        end
    end
end
scatterm(poslat, poslon, 'k.')
scatterm(neglat, neglon, 'k.')
setm(gca,'grid','on')
h = colorbar;
ylabel(h, 'SSTA')

% total # of SSTA predictors = 1343
% significantly positive # = 511
% significantly negative # = 38
% overall a consistent positive shift is observed in SSTAs around 2000

sstv = zeros(492,2);
loc = zeros(2,1);
tmark = 1;

for i = 1:ni
    for j = 1:nj
        if rmask(i,j) == 1
            sstv(:,tmark) = reshape(ssta(i,j,:),nt,1);
            loc(:,tmark) = [i;j];
            tmark = tmark + 1;
        end
    end
end



