clear
close all
format long

%% user input

h = 0.1; %size of a single grid
omega = 1.813 - 0.0025i; %resonant frequency
Ld = 5; %length design region
Lp = 0.5; %thickness of pml
radius = 1; %radius of cylinder 
res = radius / h;
epsr = 12; %permittivity of dielectric
beta = 20; %strength of pml

fprintf('Nx = %4.2f \n', res)
%% generate epsilon

L = Ld + 2*Lp;
xv = -L/2:h:L/2;
yv = xv;
N = length(xv);
[x,y] = meshgrid(xv,yv);
epsilon = ones(size(x));
flag = zeros(size(x));
for i = 1:N
    for j = 1:N
        if x(i,j)^2 + y(i,j)^2 < radius^2
            epsilon(i,j) = epsr;
            flag(i,j) = 1;
        end
    end
end

%% find mode

design_region = [Ld;Ld];
pml_thickness = [Lp;Lp];
dim = design_region + 2.*pml_thickness;
Nx = round(dim(1) / h) + 1; %num of x dim grid points
Ny = round(dim(2) / h) + 1; %num of y dim grid points

BC = {{'pml', [pml_thickness(1), beta]}, {'pml', ...
    [pml_thickness(2), beta]}}; %%boundry condition {x,y}

num_mode = 1; %% we will find # of modes with frequncy closest to the pre-defined omega
[rigeig,rigeigval] = ModeSolverTE(h,dim,BC,epsilon,num_mode,omega);
omega_mode = sqrt(rigeigval(num_mode,num_mode));

figure;
set(gcf,'position',[100,100,1000,300])
subplot(1,2,1);imagesc(xv,yv,real(reshape(rigeig(:,num_mode),Nx,Ny).')); 
colorbar; 
axis equal; 
title('Real Part of Electric Field')
xlabel('x / \lambda'); 
ylabel('y / \lambda');
subplot(1,2,2);imagesc(xv,yv,imag(reshape(rigeig(:,num_mode),Nx,Ny).')); 
colorbar; 
axis equal; 
title('Imaginary Part of Electric Field')
xlabel('x / \lambda'); 
ylabel('y / \lambda');

wavelength = (2*pi) / real(omega_mode);
display(['omega is ' num2str(omega_mode)]);
display(['wavelength is ' num2str(wavelength)]);

%% Analtyical & VIE Q factor

chi = epsr - 1;
xi = -1 ./ chi;
xy = [x(flag(:)==1) y(flag(:)==1)]; %input for green's functions
p_current = (epsr - 1) .* rigeig; %polarization current from radiated field
p = p_current(flag(:) == 1); %extract just those points within cylinder

% Greens matrices
dS = h^2;
G0EE = dS * cal_G0(xy,omega_mode);
G0Hx = dS * cal_G0HE_xz(xy, omega_mode);
G0Hy = dS * cal_G0HE_yz(xy, omega_mode);

% incident field
e_rad = G0EE * p;
e_tot = p / chi;
hx_rad = G0Hx * p;
hy_rad = G0Hy * p;

% Q VIE  & Q_ana(for reference)
Q_ana = -0.5 * (real(omega_mode) / imag(omega_mode));
U_stored = 1/4 * (epsr * norm(e_rad)^2 + 1 * (norm(hx_rad)^2 + norm(hy_rad)^2)) * dS;
P_rad = 1/2 * imag(omega_mode'*p'*G0EE*p) * dS;
Q_VIE = real(omega_mode) * (U_stored / P_rad);

%display results
error = norm(e_rad - e_tot) / norm(e_rad) * 100; %calculate percent error
fprintf('\nError in VIE = %4.2f%%\n', error); 
disp(['Prad = ', num2str(P_rad)]);
disp(['Ustore = ', num2str(U_stored)])
fprintf('Q_VIE = %4.2f\n', Q_VIE)
fprintf('Q_ana = %4.2f\n', Q_ana)
fprintf('\n')

%% begin CVX

ND = 5; %user input number of constraints
wr = real(omega_mode);
wi_vals = linspace(imag(omega_mode),0,5);
w_vals = wr + (1j .* wi_vals);
Q_bd = zeros(1,length(wi_vals));
Dmax = ones(1,length(wi_vals));
cs = cell(1,ND);

for i = 1:length(wi_vals)
    fprintf('Testing wi = %f\n',wi_vals(i))
    cs{i} = cal_DmatrixBound(epsr,ND,xy,w_vals(i),h);
    Q_bd(i) = -0.5 * (wr / wi_vals(i));
    fprintf('Qfactor Bound is %f\n',Q_bd(i))
    for j = 1:length(cs(i))
        if strcmp('Infeasible',cs{i}{j}) == 1 || strcmp('Failed',cs{i}{j}) == 1
            Dmax(i) = j-1;
            fprintf('Dmax is %d\n',Dmax(i))
        end
    end
end

figure
plot(Dmax,-1.*wi_vals)
xlabel('Max # of local power constraints yielding feasibility')
ylabel('-omega_i') 

figure
plot(-1.*wi_vals,Q_bd)
xlabel('omega_i')
ylabel('Q_bd')
