clear;clc;close all
rng(0)
n = 40;
tt = linspace(0, 2*pi, n);
tt = tt';
x = sin(tt) + 0.5*randn(n,1);
y = x + 1 + 0.5*randn(n,1);

[b, ~, ~, ~, ~] = regress(y(1:n/2), [ones(n/2,1) x(1:n/2)]);

fit = b(1) + b(2)*x; 

cc = [corr(fit(1:n/2), y(1:n/2)) corr(fit(n/2+1:end), y(n/2+1:end))]

x2 = [x(1:n/2)- mean(x(1:n/2)); x(n/2+1:end)-mean(x(n/2+1:end))];
y2 = [y(1:n/2)- mean(y(1:n/2)); y(n/2+1:end)-mean(y(n/2+1:end))];

[b, ~, ~, ~, ~] = regress(y2(1:n/2), [ones(n/2,1) x2(1:n/2)]);

fit = b(1) + b(2)*x2; 

cc = [corr(fit(1:n/2), y2(1:n/2)) corr(fit(n/2+1:end), y2(n/2+1:end))]


rcc = [];

for k = 1:25

    idx = randperm(n);
    tmpy = y2(idx);
    tmpx = x2(idx); 

    [b, ~, ~, ~, ~] = regress(tmpy(1:n/2), [ones(n/2,1) tmpx(1:n/2)]);

    fit = b(1) + b(2)*tmpx; 

    rcc = [rcc corr(fit(n/2+1:end), tmpy(n/2+1:end))];

end

