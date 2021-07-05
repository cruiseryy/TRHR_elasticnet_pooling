clear;clc;close all 

load('sst_prcp_cc.mat');

% SST index 1 for lag = 14 mo, 180E-200E, 42S-18S
lonidxp1 = 180/4+1; lonidx1 = 180 + 1;
lonidxp2 = 200/4; lonidx2 = 200;
latidxp1 = 48/4+1; latidx1 = 48+1;
latidxp2 = 72/4; latidx2 = 72;

% tmap1 = zeros(90,45); 
% tmap1(lonidxp1:lonidxp2,latidxp1:latidxp2) = 1;
% tmap2 = zeros(360,180);
% tmap2(lonidx1:lonidx2,latidx1:latidx2) = 1;
% 
% figure()
% test_pic14(tmap1)
% figure()
% test_pic14(tmap2)

sstidx1 = mean(mean(ssta(lonidx1:lonidx2,latidx1:latidx2,:)));
sstidx1 = reshape(sstidx1,492,1);

sstidxp1 = mean(mean(ssta_maxp(lonidxp1:lonidxp2,latidxp1:latidxp2,:)));
sstidxp1 = reshape(sstidxp1,492,1);

cc1 = zeros(24,1);
ccp1 = zeros(24,1);

for lag = 0:24
    tsst = sstidx1(24+5-lag:12:end-lag);
    tsstp = sstidxp1(24+5-lag:12:end-lag);
    cc1(lag+1) = corr(trhr_prcp, tsst);
    ccp1(lag+1) = corr(trhr_prcp, tsstp);
end

% SST index 2 for lag = 22 mo, 240E-272E, 22S-2N

lonidxp12 = 240/4+1; lonidx12 = 240 + 1;
lonidxp22 = 272/4; lonidx22 = 272;
latidxp12 = 68/4+1; latidx12 = 68+1;
latidxp22 = 92/4; latidx22 = 92;

% tmap3 = zeros(90,45); 
% tmap3(lonidxp12:lonidxp22,latidxp12:latidxp22) = 1;
% tmap4 = zeros(360,180);
% tmap4(lonidx12:lonidx22,latidx12:latidx22) = 1;
% 
% figure()
% test_pic22(tmap3)
% figure()
% test_pic22(tmap4)

sstidx2 = mean(mean(ssta(lonidx12:lonidx22,latidx12:latidx22,:)));
sstidx2 = reshape(sstidx2,492,1);

sstidxp2 = mean(mean(ssta_maxp(lonidxp12:lonidxp22,latidxp12:latidxp22,:)));
sstidxp2 = reshape(sstidxp2,492,1);

cc2 = zeros(24,1);
ccp2 = zeros(24,1);

for lag = 0:24
    tsst = sstidx2(24+5-lag:12:end-lag);
    tsstp = sstidxp2(24+5-lag:12:end-lag);
    cc2(lag+1) = corr(trhr_prcp, tsst);
    ccp2(lag+1) = corr(trhr_prcp, tsstp);
end


figure()
hold on
plot(0:24, cc1, 'rd-')
plot(0:24, ccp1, 'rx-')
plot(0:24, cc2, 'bd-')
plot(0:24, ccp2, 'bx-')

plot([-1 25], [0 0], 'k')
plot([-1 25], [1 1]*0.418, 'k--')
plot([-1 25], -1*[1 1]*0.418, 'k--')
grid

axis([-1 25 -0.3 0.65])

xlabel('lead time')
ylabel('Correlation Coef')
legend('SST index 1','SST index 1 w/ pooling','SST index 2','SST index 2 w/ pooling')




