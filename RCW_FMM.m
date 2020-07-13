
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RIGOROUS COUPLED WAVE ANALYSIS 
% DEVELOPED BY: IBRAHIM ISSAH. 
% DATE: 2ND JUNE, 2020
% SCHOOL:TAU UNIVERSITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all; 
clear;
clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L  = 1;       % Lattice constant 
NH = 3;       % Number of spatial harmonics %% always choose an odd number
N  = NH^2;    % Size of matrices



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DS_BOARD && GRID DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
um   = 1e-6;     % microns 
lam0 = 0.02;     % wavelength 
Lx   = 0.0175*L; % periodicity x 
Ly   = 0.0150*L;  

Nx   = 512;      % Number of grids 
Ny   = 512;
dx   = Lx/Nx; 
dy   = Ly/Ny; 

xa   = (0:Nx-1)*dx; xa = xa-mean(xa); 
ya   = (0:Ny-1)*dx; ya = ya-mean(ya);

urR  = 1.0; 
erR  = 2.0;      % dielectric of RTN
urT  = 1.0;      % permeability of URT
erT  = 9.0;      % dielectric of TRN

LA   = 2.0;      % thickness 
Ldim = [0.003 0.005]; % Thickness length; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNIT CELLS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
triangle       = 1;        % triangle unit cell 
square         = 0;        % unit cell 
ring_hollow    = 0;        %ring-hollow unit cell
ring_resonator = 0;        %square-unit cell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL % CONVMAT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
urd  = 1.0;
erd  = 6.0;
erd2 = 2.0;
 
%Initialize Layers to device
URA = urd * ones(Nx,Ny);
ER1 = erd * ones(Nx,Ny);
ERA = erd * ones(Nx,Ny);

if (square ==1 && triangle == 1 && ring_hollow == 1  && ring_resonator == 1) | (square ==0 && triangle == 0 && ring_hollow == 0  && ring_resonator == 0)
    disp('Select one lattice structure'); 
end

if square == 1 
    [M, NR] = meshgrid(ya, xa); 
    m = .005; 
    n = .006;
    r = 1.5; 
    ERA = (m./M).^2 + (n./NR).^2 < r; 
    ERA = erd + (erd2-erd)*double(ERA); 
end 

if triangle == 1
%Build Device (Layer 1)
w     = 0.8 * Ly;
h     = 0.5 * sqrt(3) * w;
ny    = round(h/dy);
ny1   = round((Ny-ny)/2);
ny2   = ny1+ny- 1;

for ny = ny1:ny2
  f = 1-(ny - ny1) / (ny2 - ny1);
  nx = round(f*w/dx);
  nx1 = 1 + floor((Nx- nx)/2);
  nx2 = nx1 + nx;
  ERA(nx1:nx2,ny) = erd2;
end
end

if ring_hollow == 1 
c0 = 3e-8; 
NR = [Nx, Ny];    
xrange = [-2 2];  % x boundaries in L0
yrange = 1*[-2 2];  % y boundaries in L0
wvlen_scan = linspace(1,2.6, 20);
wvlen = 1.7;
k0 = 2*pi/wvlen;
omega_p = 0.72*pi*1e15;%3e15; %omega_p was 3 which gave us good results...
gamma = 400e12; %20e12; % (5.5e12 is the default)
omega = 2*pi*c0/wvlen*1e6;
%epsilon_diel = 16;
epsilon_metal =  1 - omega_p^2./(omega^2-1i*gamma*omega); 
epsilon_diel = epsilon_metal;
thickness = 0.15;
fill_factor = 0.1; %half metal, half dielectric

delta_arc = 6*pi/180;
inner_radius = 0.5; outer_radius = 0.7;
eps = ones(NR);
eps = curved_stripe(eps, NR,xrange, yrange, ...
    inner_radius, outer_radius, delta_arc, epsilon_metal, epsilon_diel);

delta_arc_2 = 3*pi/180;
inner_rad_2 = 1.3; outer_rad_2 = 1.5;
eps = curved_stripe(eps, NR,xrange, yrange, ...
    inner_rad_2, outer_rad_2, delta_arc_2, epsilon_metal, epsilon_diel);
ERA = erd + (erd2-erd)*double(eps); 
end 

if ring_resonator == 1 
    NR = [Nx, Ny]; 
    inner_rad = 60; outer_rad = 90;
    xc = (round(NR(1)/2)); yc = (round(NR(2)/2));
    xa = -xc+1:xc; ya = -yc+1:yc;
    [X,Y] = meshgrid(xa,ya);
    eps_ring = ((X.^2+Y.^2)<outer_rad^2);
    eps_inner = ((X.^2+Y.^2)<inner_rad^2);
    eps_ring = eps_ring-eps_inner;
    eps_ring(eps_ring == 1) = erd;
    eps_ring(eps_ring == 0) = erd2;
    ERA = eps_ring; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STACK LAYERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ER(:, :, 1) = ERA; 
ER(:, :, 2) = ER1; 
%ER(:, :, 3) = ERA; 
%ER(:, :, 4) = ER1;

UR(:, :, 1) = URA; 
UR(:, :, 2) = URA; 
%UR(:, :, 3) = URA; 
%UR(:, :, 4) = URA; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD CONVMAT_ MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Device.ERC(N, N, length(Ldim)) = 0; 
Device.URC(N, N, length(Ldim)) = 0; 

for ll = 1: length(Ldim)
Device.ERC(:, :, ll) = convmat_PWMEM(ER(:, :, ll), NH, NH); 
Device.URC(:, :, ll) = convmat_PWMEM(UR(:, :, ll), NH, NH);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOURCE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degrees   = pi / 180;
src.theta = 0 * degrees;     % Elevation Incident Angle
src.phi   = 0 * degrees;     % Azimuthal Incident Angle
src.pte   = 1;               % Source Polarization
src.ptm   = 0;               % Source Polarization
% POLARIZATION UNIT VECTOR
n = [0; 0; 1];
k = [sin(src.theta)*cos(src.phi); sin(src.theta)*sin(src.phi); cos(src.theta)];
if(src.theta == 0)
   src.ate = [0; 1; 0]; 
else  %multiply by -1 to match bench...   check this
    src.ate = -1 * cross(n,k) ./ norm(cross(n,k));
end
src.atm = cross(src.ate,k) / norm(cross(src.ate,k));
%Compute source vector
esrcShape = zeros(sqrt(N),sqrt(N));
center = ceil(sqrt(N)/2);
%Plane Wave Example
esrcShape(center,center) = 1;
esrcShape = reshape(esrcShape,[N,1]);
%normalize
esrcShape = esrcShape / norm(esrcShape);
Esrcx = src.ate(1)*src.pte + src.atm(1)*src.ptm;
Esrcy = src.ate(2)*src.pte + src.atm(2)*src.ptm;
src.esrc = [Esrcx .* esrcShape; Esrcy .* esrcShape];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD K-MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ninc = sqrt(erR * urR);
  kinc.x = ninc * k(1);
  kinc.y = ninc * k(2);
  kinc.z = ninc * k(3);
  k0 = 2*pi / lam0;
  
  %assume equal and odd number of spatial harmonics in x,y directions
  M = (-floor(NH/2):floor(NH)/2);
  N = length(M)^2; %Size of matrices 
  
  kx = kinc.x * ones(1,sqrt(N));
  ky = kinc.y * ones(1,sqrt(N));
  for m = 1:sqrt(N)
    kx(m) = kinc.x - (2*pi*M(m))/(k0*Lx);
    ky(m) = kinc.y - (2*pi*M(m))/(k0*Ly);
  end

  Kx = zeros(N,N); 
  Ky = zeros(N,N);
  
  for m = 1:sqrt(N)
    for n = 1:sqrt(N)
      ind = (m-1)*sqrt(N) + n;
      Kx(ind,ind) = kx(n);
      Ky(ind,ind) = ky(m);
    end
  end
  
  %CHECK THIS
  KzT = conj(sqrt(erT*urT*eye(N) - Kx.^2 - Ky.^2));
  KzR = -1 * conj(sqrt(erR*urR*eye(N) - Kx.^2 - Ky.^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE INITIAL EIGENMODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kz = conj(sqrt(eye(N) - Kx.^2 - Ky.^2));

x1 = 1:N;
x2 = N+1:2*N;
y1 = x1;
y2 = x2;

Q0 = zeros(2*N,2*N);

Q0(x1,y1) = Kx * Ky;
Q0(x1,y2) = eye(N) - Kx.^2; 
Q0(x2,y1) = Ky.^2 - eye(N);
Q0(x2,y2) = -Kx * Ky;

W0        = zeros(2*N,2*N);
W0(x1,y1) = eye(N);
W0(x2,y2) = eye(N);

LAM = W0;
LAM(x1,y1) = (1i*Kz);
LAM(x2,y2) = (1i*Kz);  
V0 = Q0/(LAM);

SG.S11 = zeros(2*N,2*N);
SG.S21 = eye(2*N);
SG.S12 = eye(2*N);
SG.S22 = zeros(2*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE REF AND TRANS EIGENMODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q1 = zeros(2*N,2*N);
Q1(x1,y1) = Kx * Ky;
Q1(x1,y2) = erR*urR*eye(N) - Kx.^2; %Check this...
Q1(x2,y1) = Ky.^2 - erR*urR*eye(N);
Q1(x2,y2) = -Ky * Kx;
Q1 = Q1/urR;
Wref = zeros(2*N,2*N);
Wref(x1,y1) = eye(N);
Wref(x2,y2) = eye(N);

%The -1 are fudge factors
%kind of makes sense for going opposite direction...
lam = W0;
lam(x1,y1) = -1*(sqrt(-1)*KzR);
lam(x2,y2) = -1*(sqrt(-1)*KzR);

Vref = Q1/(lam);

A = W0\Wref + V0\Vref;
B = W0\Wref - V0\Vref;

Sref.S11 =-(A\B);
Sref.S12 = 2*(A)^(-1);
Sref.S21 = 0.5*(A - B/A*B);
Sref.S22 = (B/A); 

% TRANSMITTED
Q2 = zeros(2*N,2*N);
Q2(x1,y1) = Kx * Ky;
Q2(x1,y2) = erT*urT*eye(N) - Kx.^2; %Check this...
Q2(x2,y1) = Ky.^2 - erT*urT*eye(N);
Q2(x2,y2) = -Ky * Kx;

Q2 = Q2/urT;

Wtrn = zeros(2*N,2*N);
Wtrn(x1,y1) = eye(N);
Wtrn(x2,y2) = eye(N);
lam = W0;
lam(x1,y1) = sqrt(-1)*KzT;
lam(x2,y2) = sqrt(-1)*KzT;

Vtrn = Q2/(lam);

A = W0\Wtrn + V0\Vtrn;
B = W0\Wtrn - V0\Vtrn;
Strn.S11 = B/A;
Strn.S12 = 0.5 * (A - B/A*B);
Strn.S21 = 2*(A)^(-1);
Strn.S22 = -(A\B); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCATTERING MATRICES AND GLOBAL SCATTERING DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize Global Scattering Matrix
SG.S11 = zeros(2*N,2*N);
SG.S12 = eye(2*N);
SG.S21 = eye(2*N);
SG.S22 = zeros(2*N,2*N);


for layer = 1:length(Ldim) 
%Get parameters for current layer
  ERC = Device.ERC(:,:,layer);
  URC = Device.URC(:,:,layer); 

  Q = zeros(2*N,2*N);
  P = zeros(2*N,2*N);
  
  %Compute PQ Matrices
  P(x1,y1) = Kx/(ERC)*Ky;
  P(x1,y2) = URC - Kx/(ERC) * Kx;
  P(x2,y1) = Ky/(ERC)*Ky - URC;
  P(x2,y2) = - Ky/(ERC) * Kx;
  
  Q(x1,y1) = Kx/(URC) * Ky;
  Q(x1,y2) = ERC - Kx/(URC) * Kx;
  Q(x2,y1) = Ky/(URC) * Ky - ERC;
  Q(x2,y2) = - Ky/(URC) * Kx;
    
  %Compute Eigenmodes of layer
  %Eigenvalues /vectors can be out of order, but as long as they are
  %together, the end Scattering Matrices will be the same.
  OM2 = P*Q;
  [W, LAM2] = eig(OM2);

  LAM = sqrt(LAM2);
  V = Q*W/(LAM);

  
  %Compute S-matrix temp variables
  A = (W)\ W0 + (V)\V0;
  B = (W)\ W0 - (V)\V0;
  X = expm(-LAM * k0 * Ldim(layer));
  
  Si.S11 = (A - X*B/A*X*B)\(X*B/A*X*A - B);
  Si.S12 = (A - X*B/A*X*B)\X*(A - B/A*B);
  Si.S21 = Si.S12;
  Si.S22 = Si.S11;
  
  % SCATTERING MATRICES MULTIPLICATION
  SG = star(SG,Si);

end

SG = star(Sref,SG);
SG = star(SG,Strn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate Source Mode Coefficients
 csrc = Wref\src.esrc;

%**********************
%%   REFLECTED FIELDS
%**********************
eref = Wref * SG.S11 * csrc;
rx = eref(1:N);
ry = eref(N+1:end);
rz = -(KzR)\(Kx*rx + Ky*ry);
R = abs(rx).^2 + abs(ry).^2 + abs(rz).^2;
%Gamma = rx + ry + rz;
kzInc = cos(src.theta)*sqrt(erR * urR);
R = real(-KzR / kzInc) * R;
Gamma = reshape(R,[NH,NH]);
REF = sum(Gamma(:));


%***********************
%%   TRANSMITTED FIELDS
%***********************
etrn = Wtrn * SG.S21 * csrc;
tx = etrn(1:N);
ty = etrn(N+1:end);
tz = -(KzT) \ (Kx*tx + Ky*ty);
%Tau = tx + ty + tz; 
T = abs(tx).^2 + abs(ty).^2 + abs(tz).^2;
T = real((urR / urT) * KzT/kzInc) * T;
Tau = reshape(T,[NH,NH]);
TRN = sum(Tau(:));

%CONSERVATION LAW; 

CON = REF+TRN; 
 
h = imagesc(xa, ya, real(ER(:, :, 1))');
hh = get(h, 'Parent'); 
set(hh, 'YDir', 'normal'); 
title(['REF = ', num2str(REF),' \rightarrow', ' TRN = \rightarrow ', num2str(TRN), '\rightarrow ', ' CON = \rightarrow', num2str(CON)], 'Fontsize', 14); 
axis equal tight
colormap('jet')
