
format long

% convtest 1 
idpar = [3];
tmax = 0.25;
lambda = 0.1;
idtype = 0;
vtype = 0;
vpar = [];    
l = 6 : 9;

x = cell(4, 1);
t = cell(4, 1);
psi = cell(4, 1);
relerror = cell(4, 1);
exerror = cell(4, 1);

for i = 1 : 4
[x{i}, t{i}, psi{i}, ~, ~, ~] = ...
    sch_1d_cn(tmax, l(i), lambda, idtype, idpar, vtype, vpar);
end

for i = 1 : 3
    relerror{i} = rms(psi{i+1}(1:2:end, 1:2:end) - psi{i}, 2);
    exerror{i} = rms(exact(x{i}, t{i}) - psi{i}, 2);
end

figure(1)
plot(t{1}, relerror{1})
hold on
plot(t{2}, 4*relerror{2})
plot(t{3}, 4^2*relerror{3})
legend_opt = {"interpreter",'latex', 'location', 'southeast',"FontSize", 10};
axis_opt = {"interpreter",'latex', "FontSize", 12};
title_opt = {"interpreter",'latex', "FontSize", 14};
legend('$||d\psi^6||_2$','$4||d\psi^7||_2$','$4^2||d\psi^8||_2$',...
        legend_opt{:})
ylabel("$||d\psi^l||_2(t^n)$", axis_opt{:})
xlabel("t", axis_opt{:})
title("Convergence Test 1a: Exact Family", title_opt{:})

figure(2)
clf;
plot(t{1}, exerror{1})
hold on
plot(t{2}, 4*exerror{2})
plot(t{3}, 4^2*exerror{3})
legend('$||E\psi^6||_2$','$4||E\psi^7||_2$','$4^2||E\psi^8||_2$',...
        legend_opt{:})
ylabel("$||E\psi^l||_2(t^n)$", axis_opt{:})
xlabel("t", axis_opt{:})
title("Convergence Test 1b: Exact Family", title_opt{:})
 
% contest 2
idpar = [0.50 0.075 0.0];
tmax = 0.01;
lambda = 0.01;
idtype = 1;

for i = 1 : 4
[x{i}, t{i}, psi{i}, ~, ~, ~] = ...
    sch_1d_cn(tmax, l(i), lambda, idtype, idpar, vtype, vpar);
end

for i = 1 : 3
    relerror{i} = rms(psi{i+1}(1:2:end, 1:2:end) - psi{i}, 2);
end

figure(3)
clf;
plot(t{1}, relerror{1})
hold on
plot(t{2}, 4*relerror{2})
plot(t{3}, 4^2*relerror{3})
legend('$||d\psi^6||_2$','$4||d\psi^7||_2$','$4^2||d\psi^8||_2$',...
        legend_opt{:})
ylabel("$||d\psi^l||_2(t^n)$", axis_opt{:})
xlabel("t", axis_opt{:})
title("Convergence Test 2: Boosted Gaussian", title_opt{:})


function epsi = exact(x, t)
        [X,T] = meshgrid(x,t);
        epsi = exp(-1i * 3^2 * pi^2 * T) .* sin(3*pi * X);
end