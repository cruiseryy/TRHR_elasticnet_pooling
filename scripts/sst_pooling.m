clear;clc;close all

rng(100)

load('../data/sst/sst_1979to2019.mat')
% load('../data/sst/lsmask.mat')
load('Indian_Pacific_mask.mat')
load('../data/precip/RS/trhr_precip_2021.mat')
% -------------------------------------------------------------------------
% the wet season is defined as May to September (consistent with the wet
% season of TB from other studies and the growing season in Chen's work
% ~87% of the annual precip fall during the wet season based on 1981-2019
% (CHIRPS precip)
% -------------------------------------------------------------------------
trhr_prcp = zeros(39,1);
for k = 1:39
    trhr_prcp(k) = sum(precip((k-1)*12+5:(k-1)*12+9));
end
% -------------------------------------------------------------------------
% standardization to remove the systematic shift around 2000
% -------------------------------------------------------------------------
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
% ---------------------------------------
% some pooling types are tested
% ---------------------------------------
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
% -------------------------------
% standardization
% ----------------------------
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


sstv = zeros(492,2); % reshape to vectors of SST anomalies
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

trhr_log = (trhr_prcp>=0);
pause = 1;

% ---------------------------------------------
% to determine lambda value (at 0 lag)
% --------------------------------------------
% lag = 0;
% xx = sstv(24+5-lag:12:end-lag,:);
% testlambda = exp(linspace(log(1), log(200), 20));
% [Beta, FitInfo] = lasso(xx(1:20,:), trhr_prcp(1:20),'Lambda',testlambda, 'CV', 5, 'Alpha', 0.01);
% % [Beta, FitInfo] = lassoglm(xx(1:end-1,:), trhr_log(1:end-1), 'binomial', 'Lambda', testlambda, 'CV',5 , 'Alpha', 0.01);
% 
% lassoPlot(Beta,FitInfo,'PlotType','CV');
% % minMSE lambda = 28/21 for real amplitude predictand 
% % minMSE lambda = 12/9 for binary predictand
% pause = 1;

% ----------------------------------------------

cc = zeros(25,1);
tmap = zeros(90,45);
tmpidx = 1:39;
for lag = 0:24
    xx = sstv(24+5-lag:12:end-lag,:);

    [Beta, FitInfo] = lasso(xx(tmpidx(1:20),:), trhr_prcp(tmpidx(1:20)), 'Lambda', 28 , 'Alpha', 0.01);
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
% -------------------------------------------------------------------
% for plotting comparison between obs and est precipitation anomalies
% -------------------------------------------------------------------
% figure()
% subplot(1,5, 1:2)
% scatter(trhr_prcp(tmpidx(1:20)), tmpfit(tmpidx(1:20)), 'ro')
% hold on
% scatter(trhr_prcp(tmpidx(21:end)), tmpfit(tmpidx(21:end)), 'bd')
% grid
% title('(a)')
% ylabel('Est Precip Anomaly')
% xlabel('Obs Precip Anomaly')
% 
% subplot(1,5, 3:5)
% plot(1:39, trhr_prcp, 'r')
% hold on
% plot(1:39, tmpfit*4, 'b')
% grid
% plot([1 1]*20.5, 2.5*[-1 1], 'k--')
% plot([0 40], [0 0], 'k-')
% ylabel('Precip Anomaly')
% title('(b)')
% set(gca, 'XTick', 1:3:39)
% set(gca, 'XTickLabel', 1981:3:2019)
% 
% ---------------------------------------------------
% repeat analysis with random splitting for rpn times
% ---------------------------------------------------
% rpn = 25;
% cc = zeros(25,rpn);     
% 
% for kkk = 1:rpn
%     tmpidx = randperm(39);
%     for lag = 0:24
%         xx = sstv(24+5-lag:12:end-lag,:);
%         [Beta, FitInfo] = lasso(xx(tmpidx(1:20),:), trhr_prcp(tmpidx(1:20)), 'Lambda', 28 , 'Alpha', 0.01);
%         tmpfit = FitInfo.Intercept + xx*Beta;
%         cc(lag+1,kkk) = corr(tmpfit(tmpidx(21:end)), trhr_prcp(tmpidx(21:end)));
%         lag
%     end
% end
% 
% boxplot( cc')
% grid
% 
%----------------------------------------------
% this for testing logistic regression with elastic net 
% ---------------------------------------------
% acc = zeros(25,1);
% 
% for lag = 0:24
%     xx = sstv(24+5-lag:12:end-lag,:);
% %     
% %     [Beta, FitInfo] = lassoglm(xx(1:20,:), trhr_log(1:20),'binomial', 'Lambda', 25, 'Alpha', 0.01);
% %     tfit = glmval([FitInfo.Intercept; Beta], xx(21:end,:),'logit');
% %     acc(lag+1) = sum((tfit>=0.5) == trhr_log(21:end));
% %     pause = 1;
% %     
%     tmpfit = zeros(39,1);
%     for k = 1:39
%         txx = xx;
%         tyy = trhr_log;
%         txx(k,:) = [];
%         tyy(k,:) = [];
%         [Beta, FitInfo] = lassoglm(txx, tyy,'binomial', 'Lambda', 9, 'Alpha', 0.01);
%         tfit = glmval([FitInfo.Intercept; Beta], xx,'logit');
%         tmpfit(k) = tfit(k);
%         acc(lag+1) = acc(lag+1) + ((tfit(k)>0.5) == trhr_log(k));
%     end
% %     acc(lag+1) = sum((tmpfit>0.5).*(trhr_log==1)) + sum((tmpfit<=0.5).*(trhr_log==0));
%     pause = 1;
%     lag
% end
% 
% figure()
% bar(0:24, acc/39)
% axis([-1 25 0.4 0.75])
% hold on
% plot([-1 25], 0.5*[1 1], 'r--')
% plot([-1 25], 0.5897*[1 1], 'g--')
% grid
% ylabel('Accuracy')
% xlabel('Lead Time')

%----------------------------------------------



