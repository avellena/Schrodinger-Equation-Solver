function ctest_2d
format short
idtype = 0;
vtype = 0;
vpar = [];
idpar = [2, 3];
tmax = 0.05;
lambda = 0.05;

% function to calculate the exact solutions
function epsi = exact(x, y, t)
        epsi = zeros(length(t), length(x), length(y));
        [Y, X] = meshgrid(x,y);
        page =  sin(2*pi * X) .* sin(3*pi * Y);
        i = 1;
        for ti = t
            epsi(i, :, :) = exp(-1i*(3^2 + 2^2)*pi^2 * ti) .* page;
            i = i+1;
        end
    end

[x6, y6, t6, psi6, ~, ~, ~, ~] = ...
sch_2d_adi(tmax, 6, lambda, idtype, idpar, vtype, vpar);
[x7, y7, t7, psi7, ~, ~, ~, ~] = ...
sch_2d_adi(tmax, 7, lambda, idtype, idpar, vtype, vpar);
[x8, y8, t8, psi8, ~, ~, ~, ~] = ...
sch_2d_adi(tmax, 8, lambda, idtype, idpar, vtype, vpar);
[x9, y9, t9, psi9, ~, ~, ~, ~] = ...
sch_2d_adi(tmax, 9, lambda, idtype, idpar, vtype, vpar);


clf;
% level-to-level convergence test
figure(1)
reld6 = rms(psi7(1:2:end, 1:2:end, 1:2:end) - psi6, [2,3]);
reld7 = rms(psi8(1:2:end, 1:2:end, 1:2:end) - psi7, [2,3]);
reld8 = rms(psi9(1:2:end, 1:2:end, 1:2:end) - psi8, [2,3]);
plot(t6, reld6);
hold on
plot(t7, 4*reld7);
plot(t8, 16*reld8);
legend_opt = {"interpreter",'latex', 'location', 'southeast',"FontSize", 10};
axis_opt = {"interpreter",'latex', "FontSize", 12};
title_opt = {"interpreter",'latex', "FontSize", 14};
legend('$||d\psi^6||_2$','$4||d\psi^7||_2$','$4^2||d\psi^8||_2$',...
        legend_opt{:})
ylabel("$||d\psi^l||_2(t^n)$", axis_opt{:})
xlabel("t", axis_opt{:})
title("2D Convergence Test 1", title_opt{:})


% exact error convergence test
figure(2)
clf;
d6 = rms(psi6 - exact(x6, y6, t6), [2 3]);
plot(t6, d6);
hold on
d7 = rms(psi7 - exact(x7, y7, t7), [2 3]);
plot(t7, 4*d7);
d8 = rms(psi8 - exact(x8, y8, t8), [2 3]);
plot(t8, 16*d8);
legend('$||E\psi^6||_2$','$4||E\psi^7||_2$','$4^2||E\psi^8||_2$',...
        legend_opt{:})
ylabel("$||E\psi^l||_2(t^n)$", axis_opt{:})
xlabel("t", axis_opt{:})
title("2D Convergence Test 2", title_opt{:})
end