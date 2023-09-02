function well_survey

tmax = 0.10;
level =9;
lambda = 0.01;
idtype = 1;
idpar = [0.40, 0.075, 0.0];
vtype = 1;

dx = 1/2^level;
% find indeces of x1 and x2 in the vector
xmin = round(0.6/dx) + 1;
xmax = round(0.8/dx) + 1;

i = 1;
lnFe = linspace(2, 10, 251);
lnv0 = linspace(2, 10, 251);

for lv = linspace(0, 10, 251)
    vpar = [0.6, 0.8, -exp(lv)];

[~, ~, ~, ~, ~, ~, prob, ~] = ...
    sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

Pnx = mean(prob(:, end));
% normalize Pmin and Pmax by Pnx
Pmin = mean(prob(:, xmin))/Pnx;
Pmax = mean(prob(:, xmax))/Pnx;
lnFe(i) = log((Pmax - Pmin)/0.2);
i = i+1;
end

clf;
plot(lnv0, lnFe)
options = {'Interpreter', 'latex', 'FontSize', 12};
xlabel("ln($|V_0|$)", options{:})
ylabel("ln($\bar{F}_e$)", options{:})
title("ln($\bar{F}_e(0.6, 0.8)$) v.s. ln($V_0$)", options{:})
end