clear;clc;close all 

rng(100)

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

B = ones(4,4); B = B/sum(sum(B)); 
% B = [1 3 3 1]'*[1 3 3 1]; B = B/sum(sum(B));
% B = 1;

% C = conv2(A,B,'same')

msst1 = zeros(360,180,12); 
msst2 = zeros(360,180,12);
std_sst1 = zeros(360,180,12);
std_sst2 = zeros(360,180,12);
for k = 1:12 
    msst1(:,:,k) = mean(sst(:,:,k:12:12*22),3);
    std_sst1(:,:,k) = std(sst(:,:,k:12:12*22),[],3);
    msst2(:,:,k) = mean(sst(:,:,12*22+k:12:end),3);
    std_sst2(:,:,k) = std(sst(:,:,12*22+k:12:end),[],3);
end

ssta = zeros(360,180, 41*12);
ssta(:,:,1:12*22) = (sst(:,:,1:12*22) - repmat(msst1,1,1,22))./repmat(std_sst1,1,1,22);
ssta(:,:,12*22+1:end) = (sst(:,:,1+12*22:end) - repmat(msst2,1,1,19))./repmat(std_sst2,1,1,19);

for i = 1:41
    ssta(:,:,i) = conv2(ssta(:,:,i), B, 'same'); 
end

mask0 = zeros(360,180);
mask0(isnan(mean(ssta,3))~=1) = 1;
omask = mask0.*ocean_mask;

sstv = zeros(492,2);
loc = zeros(2,1);
tmark = 1;

for i = 1:360
    for j = 1:180
        if omask(i,j) == 1
            sstv(:,tmark) = reshape(ssta(i,j,:),41*12,1);
            loc(:,tmark) = [i;j];
            tmark = tmark + 1;
        end
    end
end

% lag = 0;
% xx = sstv(24+5-lag:12:end-lag,:);
% testlambda = exp(linspace(log(1), log(200), 20));
% [Beta, FitInfo] = lasso(xx(1:20,:), trhr_prcp(1:20),'Lambda',testlambda, 'CV', 5, 'Alpha', 0.01);
% % [Beta, FitInfo] = lassoglm(xx(1:end-1,:), trhr_log(1:end-1), 'binomial', 'Lambda', testlambda, 'CV',5 , 'Alpha', 0.01);
% 
% lassoPlot(Beta,FitInfo,'PlotType','CV');
% lambdaMinMSE = 49.6/37.5 for conv2d layer 

cc = zeros(25,1);
tmap = zeros(360,180);
tmpidx = 1:39;
for lag = 0:24
    xx = sstv(24+5-lag:12:end-lag,:);

    [Beta, FitInfo] = lasso(xx(tmpidx(1:20),:), trhr_prcp(tmpidx(1:20)), 'Lambda', 37.5 , 'Alpha', 0.01);
    tmpfit = FitInfo.Intercept + xx*Beta;
    
    cc(lag+1) = corr(tmpfit(tmpidx(21:end)), trhr_prcp(tmpidx(21:end)));
    
    for k = 1:length(Beta)
        tmap(loc(1,k), loc(2,k)) = Beta(k);
    end
%     tmap(tmap==0) = NaN;
%     figure()
%     test_pic(tmap)
%     
    lag
end

col = zeros(360,180);

% colinearily calculation window size: cwz = 3/5

% cwz = 3; 
% 
% for i = 1:360
%     for j = 1:180
%         if omask(i,j) == 1
%             
%             tmpc = -1; 
%             tmpsst = reshape(ssta(i,j,:), 492,1);
%             for ti =  (i-(cwz-1)/2): (i+(cwz-1)/2)  
%                 for tj = (j-(cwz-1)/2): (j+(cwz-1)/2)  
%                     if ti ~= i || tj ~= j 
%                         
%                        tmpsst1 = reshape(ssta(ti,tj,:), 492,1);
%                         if corr(tmpsst, tmpsst1) > tmpc
%                             tmpc = corr(tmpsst, tmpsst1); 
%                             col(i,j) = tmpc;
%                         end
%                     end
%                     
%                 end
%             end
%         end
%     end
% end




