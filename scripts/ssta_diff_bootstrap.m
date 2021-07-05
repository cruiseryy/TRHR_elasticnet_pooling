
rng(0)

nn = 1e5;
tdiff = [];
for i = 1:nn
    tmpr = randn(1,41); 
    tdiff(i) = mean(tmpr(1:22)) - mean(tmpr(23:end));
end

[quantile(tdiff,0.05) mean(tdiff) quantile(tdiff,0.95)]