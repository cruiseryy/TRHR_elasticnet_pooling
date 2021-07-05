clear;clc;close all
rng(0)

% warning('off')

load('../sst.mat')
load('../TRHR_precip.mat')

precip(:,:,457:458) = [];
sprecp = mean(mean(precip));
TRHR = zeros(456/12,1);
for k = 5:10
    TRHR = TRHR + reshape(sprecp(k:12:456),456/12,1);
end

TRHR = (TRHR - mean(TRHR))./std(TRHR);
% spi_phat = gamfit(TRHR);
% tspi = gamcdf(TRHR,spi_phat(1),spi_phat(2));
% TRHR = norminv(tspi);
    

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
% ssta(:,:,468) = [];
% ssta(:,:,1:11) = [];
score1 = zeros(13,1);


for lag = 0:12
    tmp_vector = sst_vector(12+5-lag:12:468-lag,:);
    tfu = zeros(38,1);
    tfv = zeros(38,1);
    
    for kk = 1:38
        tmpxx = tmp_vector;
        tmpxx(kk,:) = [];
        tmpyy = TRHR;
        tmpyy(kk) = [];
        
        [A,B,r,~,~,~] = canoncorr(tmpxx,tmpyy);
%         tU = (tmp_vector-repmat(mean(tmp_vector),38,1))*A;
%         tV = (TRHR-repmat(mean(TRHR),38,1))*B;
        tU = tmp_vector*A;
        tV = TRHR*B;
        tfu(kk) = tU(kk)/B;
        tfv(kk) = tV(kk)/B;
%         pause = 1;
        
    end
    score1(lag+1) = sum(sign(tfu).*sign(TRHR) > 0)/38;
%     pause = 1;
    
    
%     [A,B,r,U,V,stats] = canoncorr(tmp_vector,TRHR);
%     
%     tmap = zeros(360,180);
%     for kk = 1:length(A)
%         tmap(Loc(1,kk),Loc(2,kk)) = A(kk);
%     end
%     figure(1)
%     subplot(1,3,1:2)
%     title(strcat('lag=',num2str(lag)))
%     test_pic(tmap)
%     
%     subplot(1,3,3)
%     plot(V,U,'ro')
%     xlabel('b^tY')
%     ylabel('a^tX')
%     [lag r]
%    
%     ppause = 1;
% %     close(1)
    

    
    
end



% scatter(tfv, tfu)
% xlabel('V(precip)')
% ylabel('U(SSTA)')
% title(strcat('lag=',num2str(lag)))

