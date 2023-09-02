function barrier_survey

tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.10, 0.075, 20.0];
vtype = 1;


dx = 1/2^level;
xmin = round(0.8/dx) + 1;

i = 1;
lnFe = linspace(-2, 5, 251);
lnv0 = linspace(-2, 5, 251);
v0 = exp(lnv0);

for lv = linspace(-2, 5, 251)
    vpar = [0.6, 0.8, exp(lv)];
[~, ~, ~, ~, ~, ~, prob, ~] = ...
    sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
Pmax = mean(prob(:, end));
Pmin = mean(prob(:, xmin))/Pmax;
% Pmax = Pnx normalized to be 1
lnFe(i) = log((1 - Pmin)/0.2);
i = i+1;
end

clf;
plot(lnv0, lnFe)
options = {'Interpreter', 'latex', 'FontSize', 12};
xlabel("ln($V_0$)", options{:})
ylabel("ln($\bar{F}_e$)", options{:})
title("ln($\bar{F}_e(0.8, 1.0)$) v.s. ln($V_0$)", options{:})

end