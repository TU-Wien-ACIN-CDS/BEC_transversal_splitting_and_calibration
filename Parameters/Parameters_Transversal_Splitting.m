%% Loads necessary constants and parameters

% universal constants
hbar        = 1.054589e-34;    % Constant of Planck (J*s).
kB          = 1.38066e-23;     % Boltzmann constant (J/K). 
gEarth      = 9.81;            % Accelaration of the earth (m/s^2).
Ggrav       = 6.673e-11;       % Gravitation konstant (N*m^2/Kg^2).
cLight      = 2.99792458e8;    % Speed of light (m/s).
qelectron   = -1.602189e-19;   % Charge of electron (C).
melectron   = 9.10953e-31;     % Mass of electron (Kg).
Nmol        = 6.02205e23;      % Number of Avogadro (mol).
unucleus    = 1.6605655e-27;   % 1/12 of mass of C (Kg).
epsilon0c   = 8.854e-12;       % Dielectric constant.
mue0c       = pi*4e-7;         % Magnetic constant.
a0          = 5.291771e-11;    % Bohr radius (m).
muBohr      = 9.27408e-28;     % Bohr magneton (J/G).
gJ          = 2*1.001159657;   % Electron spin g factor.

% Features of Rb87.

MRb87   = 1.44316060e-25;      % Mass of the atom of rubidium 87 (Kg).
a21     = 95.47*a0;            % Three-dimensional s-scattering length |2,1>.
a11     = 100.44*a0;           % Three-dimensional s-scattering length |1,-1>.
a12     = 98.09*a0;            % Three-dimensional s-scattering length among |1,-1> and |2,1>.
mF      = 1;                   % Magnetic quantum number mF of the state 5S_{1/2}.
gF      = 0.5;                 % Abs Lande factor of the state 5S_{1/2} for F=1,2.
whf     = 6834.68261090434;    % Ground hyperfine splitting (MHz).
alpha0  = hbar*0.0794*(2*pi);  % Ground-state polarizability (Hz*(cm/V)^2).

% Features of Rb85.

MRb85   = 1.40999319e-25;      % Mass of the atom of rubidium 85 (Kg).
a22     = -369*a0;             % Three-dimensional t-scattering length |2,-2>.

% Features of Na.

mNa     = 3.8410e-26;          % Mass of the Na particles.
asNa    = 27.5e-10;            % Scattering length of Na.


par.m       = MRb87/hbar*(1e-6)^2/1e-3;

% spatial discretization, x in micrometers
dx  = 0.01;
x   = [-3:dx:3];
N_x = numel(x);
N_x = 2^(nextpow2(N_x)); % as a power of 2 for higher efficiency in fft
x   = linspace(min(x),max(x),N_x)';
dx  = x(2) - x(1);
L   = x(end)-x(1);

% write to parameter file
par.L       = L;
par.N_x     = N_x;
par.dx      = dx;
par.x       = x;

dt          = 5e-4;% standard time discretization
par.dt      = dt;

%% calculating the laplace operator (second derivative) and the version for the fourier space (used in split-step)
n_order         = 4;
M_matrix        = zeros(n_order+1);
for i = -n_order/2:1:n_order/2
    for BFGS_iter = 0:1:n_order
        M_matrix(i+n_order/2+1,BFGS_iter+1) = (i)^(BFGS_iter)/factorial(BFGS_iter);
    end
end
finit_difference_matrix_norm    = inv(M_matrix);
finit_difference_matrix_x       = inv(M_matrix).*(dx.^(-(0:1:n_order)'));
finit_difference_matrix_t       = inv(M_matrix).*(dt.^(-(0:1:n_order)'));
lapFilter                       = finit_difference_matrix_x(3,:);
lap4                            = speye(N_x)*lapFilter(n_order/2+1);
for i = 1:n_order/2
    indx1   = 1:1:(N_x-i);
    indx2   = (1+i):1:N_x;
    lap4(sub2ind([N_x,N_x],indx1,indx2)) = lapFilter(n_order/2+1+i);
    lap4(sub2ind([N_x,N_x],indx2,indx1)) = lapFilter(n_order/2+1+i);
end
par.lap4                        = lap4; % laplace matrix in space
temp_vals                       = 0;
wav                             = 2*pi*((1:N_x)-1)'/N_x;
for j = 1:n_order/2+1
    temp_vals   = temp_vals + finit_difference_matrix_x(3,j)*cos((n_order/2-j+1)*wav) * (1 + ((n_order/2-j+1)>0));
end
ilap                            = temp_vals;
par.ilap                        = ilap; % laplace operator in fourier space

lap_scaled                      = -1/2/par.m*lap4;
ilap_scaled                     = -1/2/par.m*ilap;

%% calculate interaction constant g_perp/g1D
N       = 4000;                         % expected number of atoms

as      = a11;                          % 3D atom-atom scattering length (m).

wx      = 2*pi*16.3;                    % Axial frequency of the trap (Hz). (for a dressing amplitude of 0.3)
wy      = 2*pi*1.83e3;                  % Frequency along the (transverse) y direction (Hz).
wz      = 2*pi*2.58e3;                  % Frequency along the (transverse) z direction (Hz).
wr      = sqrt(wy*wz);                  % Transverse (or radial) mean frequency (Hz).
ax      = sqrt(hbar/MRb87/wx);          % Axial oscillator length (m).
az      = sqrt(hbar/MRb87/wz);          % Oscillator length along the (transverse) z direction (m).
ar      = sqrt(hbar/MRb87/wr);          % Transverse mean oscillator length (m).

falpha  = @(a)a^3*(a+5)^2-(15*N*ar*as/ax^2)^2;  
alpha   = fzero(falpha,0.5);

L       = ax^2*sqrt(alpha)/ar;          % Condensate length (m).

% Coupling constant g1D (J*m) for the GPE along the x direction. 
g1DJM   = 2*sqrt(2*pi)*hbar^2/MRb87/az*alpha^2*L/315*(alpha^2+9*alpha+21)/as/N^2;

g1D     = g1DJM*1e-3*1e6; % in kJ * um
% disp(['g1D = ' num2str(g1D) 'kJ*/mu m']); 

%% normalization values
T_norm              = 1e-3;
Z_norm              = 1e-6;
G_norm              = hbar * Z_norm / (T_norm * N);
E_norm_J            = N*hbar/T_norm;
E_norm_Hz           = N/T_norm/2/pi;
Psi_norm            = sqrt(N/Z_norm);


g1D_first_principle = g1D * 1e-3 / G_norm; % times 1e-3 to get to J/m units

