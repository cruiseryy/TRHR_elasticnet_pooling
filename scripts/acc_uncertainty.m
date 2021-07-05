clear;clc;close all

acc = [];

for k = 1:1e4
    
    tmpacc = 0;
    for i = 1:39
        tmpacc = tmpacc + (rand()<=0.5);
    end
    acc = [acc tmpacc/39];
end    

quantile(acc,0.9) %0.5897 p = 0.1   
