function [x, y, t, psi, psire, psiim, psimod, v] = ...
sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)
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
% y: Vector of y coordinates [ny]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx x ny]
% psire Array of computed psi_re values [nt x nx x ny]
% psiim Array of computed psi_im values [nt x nx x ny]
% psimod Array of computed sqrt(psi psi*) values [nt x nx x ny]
% v Array of potential values [nx x ny]

nx = 2^level + 1;
ny = nx;
x = linspace(0.0, 1.0, nx)';
y = linspace(0.0, 1.0, ny);
dx = 2^(-level);
dy = dx;
dt = lambda * dx;
nt = round(tmax / dt) + 1;
t = [0 : nt-1] * dt;

% Initialize solution, and set initial data
psi = zeros(nt, nx, ny);
v = zeros(nx, ny);
if idtype == 0
    mx = idpar(1);
    my = idpar(2);
    psi(1, :, :) = sin(mx*pi*x) .* sin(my*pi*y);
elseif idtype == 1
    x0 = idpar(1);
    y0 = idpar(2);
    delx = idpar(3);
    dely = idpar(4);
    px = idpar(5);
    py = idpar(6);
    psi(1,:, :) = exp(1i*px*x).* exp(1i*py*y) .* exp(-(((x-x0).^2/delx^2) + ...
        (y-y0).^2/dely^2));
else
    fprintf('diff_1d_imp: Invalid idtype=%d\n', idtype);
    return
end
% Initialize potential
if vtype == 0
elseif vtype == 1
    xfrom = ceil(vpar(1)/dx) + 1;
    xto = ceil(vpar(2)/dx);
    yfrom = ceil(vpar(3)/dy) + 1;
    yto = ceil(vpar(4)/dy);
    v(xfrom:xto, yfrom:yto) = vpar(5);
elseif vtype == 2
    jfrom = (ny - 1)/4 + 1;
    x1 = ceil(vpar(1)/dx) + 1;
    x2 = ceil(vpar(2)/dx);
    x3 = ceil(vpar(3)/dx) + 1;
    x4 = ceil(vpar(4)/dx);
    v(:, jfrom:jfrom+1) = vpar(5);
    v(x1: x2, :) = 0;
    v(x3: x4, :) = 0;
else
    fprintf('sch_1d_cn: Invalid vtype=%d\n', vtype);
    return
end

% Set up tridiags for stage 1 
dl = -1i*dt / (2*dx^2) * ones(nx, 1);
d  = (1i*dt/dx^2 + 1) * ones(nx,1);
du = dl;
% Fix up boundary cases ...
d(1) = 1.0;
du(2) = 0.0;
dl(nx-1) = 0.0;
d(nx) = 1.0;
% Define sparse matrix ...
A = spdiags([dl d du], -1:1, nx, nx);

% set up RHS matrix for step 1
dr = (1 - 1i*dt/dy^2) *ones(nx, 1); 
dr(1) = 1;
dr(end) = 1;
ARHS = spdiags([-dl dr -du], -1:1, nx, nx);

for n = 1 :  nt-1
    psi_ij = squeeze(psi(n, :, :));  % nx x ny array
    F = zeros(nx, ny);      % initialize RHS for stage 1
    

    for xi = 2 : nx - 1
        F(xi, 2:ny-1) = 1i*dt/(2*dy^2) * psi_ij(xi, 1:ny-2) + ...
                  (1 - 1i*dt/2*v(xi, 2: ny-1) - 1i*dt/dy^2) .* psi_ij(xi, 2: ny-1) + ...
                  1i*dt/(2*dy^2) * psi_ij(xi, 3:ny);
    end    
    
    % final RHS matrix for Stage 1 tridiag solve
    F = ARHS * F;    
    psi_temp = zeros(nx, ny); 
    
    % Stage 1 tridiag solves, don't need to loop over y because A does not
    % depend on V; can solve for nx x ny matrix directly
    psi_temp(:, 2:ny-1) = A \ F(:, 2:ny-1);
   
    
    % Stage 2: B needs different V at different x; loop over x and solve
    % for all y's at each x
    for xi = 2 : nx-1
        dB = d + 1i*dt*v(xi, :)'/2;
        dB(1) = 1;
        dB(end) = 1;
        B = spdiags([dl dB du], -1:1, nx, nx);
        psi(n+1, xi, :) = B \ (psi_temp(xi, :).');
    end
    
end

psire = real(psi);
psiim = imag(psi);
psimod = abs(psi);

end








