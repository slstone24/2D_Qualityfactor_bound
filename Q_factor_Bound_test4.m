clear
close all
format long

%% user input
h = 0.1; %size of a single grid
% omega = 1.813 - 0.0025i; %resonant frequency
omega = 2.4;
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


%% find mode(s) & resonant frequencies

design_region = [Ld;Ld];
pml_thickness = [Lp;Lp];
dim = design_region + 2.*pml_thickness;
Nx = round(dim(1) / h) + 1; %num of x dim grid points
Ny = round(dim(2) / h) + 1; %num of y dim grid points

BC = {{'pml', [pml_thickness(1), beta]}, {'pml', ...
    [pml_thickness(2), beta]}}; %%boundry condition {x,y}

num_mode = 10; %number of resonant frequencies to solve for
[rigeig,rigeigval] = ModeSolverTE(h,dim,BC,epsilon,num_mode,omega);

omega_vals = zeros(1,num_mode);
for i = 1:num_mode
    omega_vals(i) = sqrt(rigeigval(i,i));
end

omega_vals = sort(omega_vals);
%% Analtyical & VIE Q factor

xy = [x(flag(:)==1) y(flag(:)==1)]; %input for green's functions
Q_VIE = ones(1,num_mode);
Q_ana = ones(1,num_mode);

fprintf('Now displaying analytical and VIE Q factor results \n\n\n')

for i = 1:num_mode
    fprintf('\nOmega = %4.4f%4.4fi \n', real(omega_vals(i)),imag(omega_vals(i)))
    Q_VIE(i) = QVIE(epsr,rigeig(:,i),xy,h,omega_vals(i),flag);
    Q_ana(i) = -0.5 * (real(omega_vals(i)) / imag(omega_vals(i)));
    fprintf('Q_ana = %4.4f \n', Q_ana(i))
end

%% begin CVX

ND = 4; %user input number of constraints
n_samples = 5;
cs = cell(num_mode,n_samples);
w = zeros(num_mode,n_samples);
wr = sort(real(omega_vals));
wi = zeros(num_mode,n_samples);
wi(:,end) = imag(omega_vals)';

for i = 1:num_mode
    wi(i,:) = 1i.*linspace(0,imag(omega_vals(i)),n_samples);
end

for i = 1:num_mode
    fprintf('\n\n')
    for j = 1:n_samples
        w(i,j) = wr(i) + wi(i,j);
        fprintf('w = %f %fi\n', real(w(i,j)),imag(w(i,j)))
        fprintf('Testing wi = %f\n',imag(wi(i,j)))
        cs{i,j} = cal_DmatrixBound(epsr,ND,xy,w(i,j),h);
        fprintf('\n')
    end
end

%% extract Q bound

Q_bd = zeros(1,num_mode);
fprintf('\n')

% loop through feasibility values and set to empty if infeasible or failed
for i = 1:num_mode
    for j = 2:n_samples
            if any(strcmp(cs{i,j}, 'Infeasible')) || any(strcmp(cs{i,j}, 'Failed')) || isempty(cs{i,j})
               cs{i,j} = {};
            end
    end
end

 
%calculate Q_bd from minimum feasible wi
for i = 1:num_mode
    for j = 2:n_samples
        if ~isempty(cs{i,j})
            Q_bd(i) = -0.5 * ((wr(i)) / imag(wi(i,j)));
            fprintf('Q_bd = %f\n',Q_bd(i))
            break
        else
            fprintf('its empty\n')
            continue
        end     
    end
end

% %plot results
figure
plot(sort(real(omega_vals)),Q_ana)
hold on
plot(sort(real(omega_vals)),Q_bd,'--ro')
plot(sort(real(omega_vals)),Q_VIE,'--go');
xlabel('Real(Omega)')
ylabel('Qfactor')
