clear;clc;close all

load('tmp_coef.mat')

bbeta = [Beta14 Beta22 Betabi14 Betabi22];

title_str = {'lag = 14','lag = 22','lag = 14 bi', 'lag = 22 bi'};

% for i = 1:4
%     
%     figure(i)
%     tmpbeta = bbeta(:,i);
%     tmap = zeros(90,45);
%     for k = 1:size(loc,2)
%         tmap(loc(1,k),loc(2,k)) = tmpbeta(k);
%         
%     end
%     
%     test_pic(tmap)
%     title(title_str{i})
%     
% end

%% lag = 14 mo
figure(1)
tmpbeta = bbeta(:,1);
tmap = zeros(90,45);

for k = 1:size(loc,2)
    tmap(loc(1,k),loc(2,k)) = tmpbeta(k);
end

test_pic14(tmap)
title(title_str{1})

figure(3)
tmpbeta = bbeta(:,3);
tmap = zeros(90,45);

for k = 1:size(loc,2)
    tmap(loc(1,k),loc(2,k)) = tmpbeta(k);
end

test_pic14(tmap)
title(title_str{3})

cc14 = corr(Beta14,Betabi14);

mask14 = (((Beta14~=0) + (Betabi14~=0))~=0);
ccc14 = corr(Beta14(mask14~=0), Betabi14(mask14~=0));
figure(5)
scatter(Beta14, Betabi14, 'r.')

%% lag = 22 mo
figure(2)
tmpbeta = bbeta(:,2);
tmap = zeros(90,45);

for k = 1:size(loc,2)
    tmap(loc(1,k),loc(2,k)) = tmpbeta(k);
end

test_pic22(tmap)
title(title_str{2})

figure(4)
tmpbeta = bbeta(:,4);
tmap = zeros(90,45);

for k = 1:size(loc,2)
    tmap(loc(1,k),loc(2,k)) = tmpbeta(k);
end

test_pic22(tmap)
title(title_str{4})

cc22 = corr(Beta22,Betabi22);

mask22 = (((Beta22~=0) + (Betabi22~=0))~=0);
ccc22 = corr(Beta22(mask22~=0), Betabi22(mask22~=0));
figure(6)
scatter(Beta22, Betabi22, 'r.')