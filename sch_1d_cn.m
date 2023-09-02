function [x, t, psi, psire, psiim, psimod, prob, v] = ...
    sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Inputs
%
% tmax: Maximum integration time
% level: Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar: Vector of initial data parameters
% vtype: Selects potential type
% vpar: Vector of potential parameters
%
% Outputs
%
% x: Vector of x coordinates [nx]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx]
% psire Array of computed psi_re values [nt x nx]
% psiim Array of computed psi_im values [nt x nx]
% psimod Array of computed sqrt(psi psi*) values [nt x nx]
% prob Array of computed running integral values [nt x nx]
% v Array of potential values [nx]

nx = 2^level + 1;
x = linspace(0.0, 1.0, nx);
dx = x(2) - x(1);
dt = lambda * dx;
nt = round(tmax / dt) + 1;
t = [0 : nt-1] * dt;

% Initialize solution, and set initial data
psi = zeros(nt, nx);
v = zeros(nx, 1);
prob = zeros(nt, nx);
if idtype == 0
    m = idpar(1);
    psi(1, :) = sin(m*pi*x) + 0i;
elseif idtype == 1
    x0 = idpar(1);
    d = idpar(2);
    p = idpar(3);
    psi(1,:) = exp(1i*p*x).*exp(-((x-x0)/d).^2);
else
    fprintf('diff_1d_imp: Invalid idtype=%d\n', idtype);
    return
end

% Initialize potential
if vtype == 0
elseif vtype == 1
    from = ceil(vpar(1)/dx) + 1;
    to = ceil(vpar(2)/dx);
    v(from:to, 1) = vpar(3);
else
    fprintf('sch_1d_cn: Invalid vtype=%d\n', vtype);
    return
end


% Initialize storage for sparse matrix and RHS ...
dl = zeros(nx,1);
d  = zeros(nx,1);
du = zeros(nx,1);
f  = zeros(nx,1);

% Set up tridiagonal system ...
dl = 1.0 / (2*dx^2) * ones(nx, 1);
d  = (1i / dt - 1.0 / dx^2) * ones(nx,1) - v/2;
du = dl;
% Fix up boundary cases ...
d(1) = 1.0;
du(2) = 0.0;
dl(nx-1) = 0.0;
d(nx) = 1.0;
% Define sparse matrix ...
A = spdiags([dl d du], -1:1, nx, nx);

% Set up sparse matrix for RHS
dr  = (1i / dt + 1.0 / dx^2) * ones(nx,1) + v/2;
dr(1) = 1.0;
dr(nx) = 1.0;
B = spdiags([-dl dr -du], -1:1, nx, nx);

for n = 1 : nt-1
    % Define RHS of linear system ...
    f = B * psi(n, :).';
    % Solve system, thus updating approximation to next time
    % step ...
    psi(n+1, :) = A \ f;
end

psire = real(psi);
psiim = imag(psi);
psimod = abs(psi);
% calculate running integral along x dimension, produces nt x nx array
prob = cumtrapz(psi .* conj(psi), 2);

end


