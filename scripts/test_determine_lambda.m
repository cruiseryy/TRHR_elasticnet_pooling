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

figure(1)
hold on
plot(1981:2019, trhr_prcp,'ro-')
plot([1 1]*2000.5, 2.5*[-1 1], 'k--')
grid
ylabel('wet season precip anomaly')
close(1)

% standardization of SST
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
ssta = zeros(360,180,492);
ssta(:,:,1:12*22) = (sst(:,:,1:12*22) - repmat(msst1,1,1,22))./repmat(std_sst1,1,1,22);
ssta(:,:,12*22+1:end) = (sst(:,:,1+12*22:end) - repmat(msst2,1,1,19))./repmat(std_sst2,1,1,19);

sstv = zeros(492,2);
loc = zeros(2,1);
tmark = 1;
for i = 1:360
    for j = 1:180
%         if lsmask(i,j) == 1
        if ocean_mask(i,j) == 1
            sstv(:,tmark) = reshape(ssta(i,j,:),492,1);
            loc(:,tmark) = [i; j];
            tmark = tmark + 1;
        end
    end
end

% pre-determine lambda using full series data and lag of 0 

trhr_log = (trhr_prcp>=0);
pause = 1;
%----------------------------------------------
% 
% lag = 0;
% xx = sstv(24+5-lag:12:end-lag,:);
% testlambda = exp(linspace(log(1), log(200), 20));
% % [Beta, FitInfo] = lasso(xx(1:20,:), trhr_prcp(1:20),'Lambda',testlambda, 'CV', 5, 'Alpha', 0.01);
% % [Beta, FitInfo] = lassoglm(xx, trhr_log, 'binomial', 'Lambda', testlambda, 'CV',5 , 'Alpha', 0.01);
% 
% lassoPlot(Beta,FitInfo,'PlotType','CV');
% % minMSE lambda = 37.5 (~38) for real amplitude predictand
% % minMSE lambda = 12/16 for binary predictand
% pause = 1;
% 
%----------------------------------------------
% 
% lag = 8;
% xx = sstv(24+5-lag:12:end-lag,:);
% [Beta, FitInfo] = lassoglm(xx, trhr_log, 'binomial', 'CV', 3, 'Alpha', 0.01);
% [Beta, FitInfo] = lassoglm(xx(1:19,:), trhr_log(1:19), 'binomial', 'Lambda', 12, 'Alpha', 0.01);
% tfit = glmval([FitInfo.Intercept; Beta], xx,'logit');
% tfitlog = (tfit>=0.5);
% acc = sum(tfitlog(20:end) == trhr_log(20:end));
% pause = 1;
% 
%----------------------------------------------

cc = zeros(25,1);
% tmplambda = exp(linspace(log(17),log(53),10));
% for k = 1:10
for lag = 0:24
    xx = sstv(24+5-lag:12:end-lag,:);

    [Beta, FitInfo] = lasso(xx(1:20,:), trhr_prcp(1:20,:), 'Lambda',38 , 'Alpha', 0.01);
    tmpfit = FitInfo.Intercept + xx*Beta;

%     tmpfit = zeros(39,1);
%     for k = 1:39
%         txx = xx;
%         tyy = trhr_prcp;
%         txx(k,:) = [];
%         tyy(k,:) = [];
%         [Beta, FitInfo] = lasso(txx, tyy, 'Lambda', 10, 'Alpha', 0.01);
%         tfit = FitInfo.Intercept + xx(k,:)*Beta; tmpfit(k) = tfit;
%     end
%     
%     cc(lag+1) = corr(tmpfit, trhr_prcp);
    cc(lag+1) = corr(tmpfit(21:end), trhr_prcp(21:end));
    lag
end
% end
% 
%----------------------------------------------
% 
% acc = zeros(25,1);
% % load('Lambda.mat')
% % mlambda = median(LLambda);
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
%         [Beta, FitInfo] = lassoglm(txx, tyy,'binomial', 'Lambda', 16, 'Alpha', 0.01);
%         tfit = glmval([FitInfo.Intercept; Beta], xx,'logit');
%         tmpfit(k) = tfit(k);
%         acc(lag+1) = acc(lag+1) + ((tfit(k)>0.5) == trhr_log(k));
%     end
% %     acc(lag+1) = sum((tmpfit>0.5).*(trhr_log==1)) + sum((tmpfit<=0.5).*(trhr_log==0));
%     lag
% end
% 
%----------------------------------------------
% 
% trhrw = (trhr_prcp > 0.4147);
% trhrm = (trhr_prcp <= 0.4147).* (trhr_prcp >= -0.4853);
% trhrd = (trhr_prcp < -0.4853);
% yyw = trhrw;
% yym = trhrm;
% yyd = trhrd;
% trhrc = zeros(39,1);
% trhrc(yyw == 1) = 1;
% trhrc(yyd == 1) = 2; 
% trhrc(yym == 1) = 3;
% acc = zeros(25,1);

%----------------------------------------------
% 
% lag = 0;
% xx = sstv(24+5-lag:12:end-lag,:);
% testlambda = exp(linspace(log(1), log(200), 20));
% % [Beta, FitInfo] = lasso(xx(1:20,:), trhr_prcp(1:20),'Lambda',testlambda, 'CV', 5, 'Alpha', 0.01);
% % [Beta, FitInfo] = lassoglm(xx(1:end-1,:), trhr_log(1:end-1), 'binomial', 'Lambda', testlambda, 'CV',5 , 'Alpha', 0.01);
% [Beta, FitInfo] = lassoglm(xx, trhrw,'poisson','Lambda',testlambda, 'CV', 5, 'Alpha', 0.01);
% lassoPlot(Beta,FitInfo,'PlotType','CV');
% % minMSE lambda = 35 for real amplitude predictand
% % minMSE lambda = 12/16 for binary predictand
% pause = 1;
% 
%----------------------------------------------
% 
% for lag = 0:24
%     tic
%     xx = sstv(24+5-lag:12:end-lag,:);
%     tpr = [];
%     for k = 1:39
%         txx = xx; txx(k,:) = [];
%         tyyw = yyw; tyyw(k) = [];
%         tyyd = yyd; tyyd(k) = [];
%         tyym = yym; tyym(k) = [];
%         
%         [b1, f1] = lassoglm(txx, tyym,'poisson','Lambda', 12, 'Alpha', 0.01);
%         [b2, f2] = lassoglm(txx, tyyw,'poisson','Lambda', 12, 'Alpha', 0.01);
%         [b3, f3] = lassoglm(txx, tyyd,'poisson','Lambda', 12, 'Alpha', 0.01);
%         wtm = xx(k,:)*(b2 - b1) + f2.Intercept - f1.Intercept; wtm = exp(wtm);
%         dtm = xx(k,:)*(b3 - b1) + f3.Intercept - f1.Intercept; dtm = exp(dtm);
%         tmp_pr = [wtm dtm 1]./(wtm + dtm + 1); % wet dry mild
%         tpr = [tpr; tmp_pr];
%         tmp_case = find(tmp_pr == max(tmp_pr)); 
%         pause = 1;
%         acc(lag+1) = acc(lag+1) + (tmp_case == trhrc(k));
%     end
%     [lag toc]
% end




