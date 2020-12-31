clear
close all
format long

%% user input

h = 0.1; %size of a single grid
omega = 1.813 - 0.0025i; %resonant frequency
% omega = 0;
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

num_mode = 20; %number of resonant frequencies to solve for
[rigeig,rigeigval] = ModeSolverTE(h,dim,BC,epsilon,num_mode,omega);

omega_vals = zeros(1,num_mode);
for i = 1:num_mode
    omega_vals(i) = sqrt(rigeigval(i,i));
end

%% Analtyical & VIE Q factor

xy = [x(flag(:)==1) y(flag(:)==1)]; %input for green's functions
Q_VIE = ones(1,num_mode);
Q_ana = ones(1,num_mode);
for i = 1:num_mode
    Q_VIE(i) = QVIE(epsr,rigeig(:,i),xy,h,omega_vals(i),flag);
    Q_ana(i) = -0.5 * (real(omega_vals(i)) / imag(omega_vals(i)));
end

%% begin CVX

ND = 5; %user input number of constraints
step = 20;
cs = cell(num_mode,step);
w = zeros(num_mode,step);
wr = sort(real(omega_vals));
wi = 1i.*linspace(0,imag(omega),step);

for i = 1:num_mode
    fprintf('\n\n')
    for j = 1:length(wi)
        w(i,j) = wr(i) + wi(j);
        fprintf('w = %f %fi\n', real(w(i,j)),imag(w(i,j)))
        fprintf('Testing wi = %f\n',imag(wi(j)))
        cs{i,j} = cal_DmatrixBound(epsr,ND,xy,w(i,j),h);
        fprintf('\n')
    end
end

%% extract Q bound

Q_bd = zeros(num_mode,step);
fprintf('\n')

% loop through feasibiltiy values and store calculated Q
for i = 1:num_mode
    fprintf('\n')
    for j = 1:length(wi)
            if ~any(strcmp(cs{i,j},'Infeasible')) || ~any(strcmp(cs{i,j},'Failed'))
                Q_bd(i,j) = -0.5 * (real(w(i,j)) / imag(wi(j)));
                fprintf('Q_bd = %f\n',Q_bd(i,j))
            else
                Q_bd(i,j) = Q_ana(i);
                fprintf('Q_bd = Q_ana = %f\n',Q_bd(i,j))
            end
    end
end

% determine minimum wi that has returned feasible
Dmax = zeros(num_mode,step); %flag when infeasible or failed
for i = 1:num_mode
    for j = 2:step
        for k = 1:ND
            if strcmp('Infeasible', cs{i,j}{k}) || strcmp('Failed', cs{i,j}{k})
                Dmax(i,j) = 0;
            elseif strcmp('Solved',cs{i,j}{k})
                Dmax(i,j) = 1;
            end
        end
    end
end

%extract upper Q bound for each resonant frequency
Q = Q_ana;
for i = 1:num_mode
    for j = 2:step
            if Dmax(i,j) == 1
                Q(i) = Q_bd(i,j);
                break
            else
                continue
            end
    end
end

for i = 1:num_mode
    if Q(i) < Q_ana(i)
        Q(i) = Q_ana(i);
    end
end

%plot results
figure
plot(sort(real(omega_vals)),Q_ana)
hold on
plot(sort(real(omega_vals)),Q)
xlabel('Real(Omega)')
ylabel('Qfactor')
