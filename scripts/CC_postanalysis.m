% clear;clc;close all
clear;clc

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

% max pooling
% average pooling 
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


Fs = 1; 
T = 1/Fs;
y = trhr_prcp;
L = length(y);
tt = (0:L-1)*T;

NFFT = 2^nextpow2(L);
Y0 = fft(y,NFFT)/L;
f0 = Fs/2*linspace(0,1,NFFT/2+1);
% loglog(1./f0,2*abs(Y0(1:NFFT/2+1))) 
% grid

% lag = 12+10;
cclag = [];
cc90 = [];
for lag = 22
% tmpY = [];
    cc = zeros(ni, nj);
    cc1 = zeros(ni, nj);
    cc2 = zeros(ni, nj);
    % figure()
    % hold on
    for k = 1:size(sstv,2)
        y = sstv(24+5-lag:12:end-lag,k);
    %     Y = fft(y,NFFT)/L;
    %     tmpY = [tmpY abs(Y(1:NFFT/2+1))];
    %     f = Fs/2*linspace(0,1,NFFT/2+1);

        cc(loc(1,k), loc(2,k)) = corr(y,trhr_prcp);
        cc1(loc(1,k), loc(2,k)) = corr(y(1:20),trhr_prcp(1:20));
        cc2(loc(1,k), loc(2,k)) = corr(y(21:end),trhr_prcp(21:end));
        pause = 1;
    %     p1 = plot(1./f0,2*abs(Y(1:NFFT/2+1)),'g');
    %     p1.Color(4) = 0.1;
    end
    ccv1 = reshape(cc1(cc1~=0),[],1);
    ccv2 = reshape(cc2(cc2~=0),[],1);
    cc90 = [cc90; quantile(reshape(abs(cc(cc~=0)),[],1),0.99)];
    cclag = [cclag; corr(ccv1,ccv2)];
    
    figure()
%     subplot(2,1,1)
    test_pic(cc1)
    title('81-00')
%     subplot(2,1,2)
    figure()
    test_pic(cc2)
    title('01-19')

    pause = 1;
end
%-----------
% cc = [[-0.128070175438596;-0.0929824561403509;0.0456140350877193;-0.0824561403508772;-0.0771929824561403;-0.0736842105263158;-0.105263157894737;0.107017543859649;0.231578947368421;-0.0666666666666667;0.0912280701754386;0.142105263157895;-0.0315789473684211;0.391228070175439;0.610526315789474;0.166666666666667;0.187719298245614;-0.0280701754385965;0.0175438596491228;0.342105263157895;0.0824561403508772;0.243859649122807;0.621052631578947;0.387719298245614;0.182456140350877]];
%-----------
% pause = 1;
% tmpy1 = quantile(tmpY,0.25,2);
% tmpy2 = quantile(tmpY,0.75,2);
% tmpyavg = mean(tmpY,2);
% figure()
% hold on
% plot(1./f0, 2*tmpy1, 'g')
% plot(1./f0, 2*tmpy2, 'g')
% plot(1./f0, 2*tmpyavg, 'b')
% plot(1./f0,2*abs(Y0(1:NFFT/2+1)),'k')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% grid

% figure()
% test_pic(cc1)
% title(strcat('lag=',num2str(lag),'mo\_1981-2000'))
% figure()
% test_pic(cc2)
% title(strcat('lag=',num2str(lag),'mo\_2001-2019'))
% 
% figure()
% test_pic(cc)
% title(strcat('lag=',num2str(lag),'mo\_1981-2019'))

% load('../data/sst/nino34.mat')
% 
% cc_nino = zeros(25,1);
% cc_ninof = zeros(25,1); % first half period 1:20
% cc_ninos = zeros(25,1); % second half period 21:end
% for lag = 0:24
%     tmp_nino_idx = nino34(24+5-lag:12:end-lag);
%     cc_nino(lag+1) = corr(tmp_nino_idx, trhr_prcp);
%     cc_ninof(lag+1) = corr(tmp_nino_idx(1:20), trhr_prcp(1:20));
%     cc_ninos(lag+1) = corr(tmp_nino_idx(21:end), trhr_prcp(21:end));
% end
% 
% cc_nino2 = zeros(25,1);
% cc_ninof2 = zeros(25,1); % first half period 1:20
% cc_ninos2 = zeros(25,1); % second half period 21:end
% for lag = 0:24
%     tmp_nino_idx = nino34_2(24+5-lag:12:end-lag);
%     cc_nino2(lag+1) = corr(tmp_nino_idx, trhr_prcp);
%     cc_ninof2(lag+1) = corr(tmp_nino_idx(1:20), trhr_prcp(1:20));
%     cc_ninos2(lag+1) = corr(tmp_nino_idx(21:end), trhr_prcp(21:end));
% end
% % 
% figure()
% bar(0:24, [cc_ninof,cc_ninos])
% legend('81-00', '01-19')
% axis([-1 25 -0.3 0.4])

pause = 1;

load('cc_realamp.mat')
figure()
bar(0:24, cc_lasso)
hold on
plot(0:24, cclag, 'bd-','MarkerSize',8)
ylabel('Correlation Coef')
xlabel('Lead Time (mo)')
grid
axis([-1 25 -0.2 0.7])

% load('elastic_net_cc_max_pooling.mat')
% bar(0:24,[cc cclag])
% legend('regression','cc analysis') 
    
    
