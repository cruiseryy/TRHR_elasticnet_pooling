clear;clc;close all 
rng(0)

acc = zeros(1000,1);

for k = 1:1e5
    mark = 0;
    for i = 1:39
        mark = mark + (rand()>=0.5);
    end
    
    acc(k) = mark/39;
end 