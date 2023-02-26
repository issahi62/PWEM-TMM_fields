function C = convmat_PWMEM(A, P, Q, R)
% CONVMAT  Rectangular Convolution Matrix 
% A = matrix for fiding its convolution 
% P, Q, R are the spatial harmonics; That is the orders of each mode (-1 0,
% 1).. Remember to choose odd numbers. 
%
% C = convmat_PWMEM(A, P); for 1D problems  
% C = convmat_PWMEM(A, P, Q); for 2D problems 
% C = convmat_PWMEM(A, P, Q, R); for 3D problems 

% This function constructs convolution matrices 
% real space grid. 

%   DETERMINE THE SIZE OF A 
[Nx, Ny, Nz] = size(A); 

% HANDLE NUMBER OF HARMONICS FOR ALL DISPERSIONS
if nargin == 2
    Q =1; 
    R = 1; 
    
elseif nargin == 3
    R = 1; 
end 
% COMPUTE INDICES OF SPATIAL HARMONICS
NH = P*Q*R;                   %total number

p= [-floor(P/2): floor(P/2)]; %indices along x
q= [-floor(Q/2): floor(Q/2)]; % indices along y
r= [-floor(R/2): floor(R/2)]; %indices along z


%COMPUTE FOURIER COEFFICIETNTS OF A
A = fftshift(fftn(A))/(Nx*Ny*Nz); 

% COMPUTE ARRAY INDICES OF CENTER HARMONIC
pO = 1+floor(Nx/2); 
qO = 1+floor(Ny/2); 
rO = 1+floor(Nz/2); 

% FILLING ELEMENTS OF CONVOLUTION MATRIX

for rrow = 1:R
for qrow = 1:Q
for prow = 1:P
    row = (rrow-1)*Q*P + (qrow-1)*P + prow;
    for rcol = 1:R
    for qcol = 1:Q
    for pcol = 1:P
        col = (rcol-1)*Q*P + (qcol-1)*P + pcol;
        pfft = p(prow)-p(pcol); 
        qfft = q(qrow)-q(qcol); 
        rfft = r(rrow)-r(rcol); 
        C(row,col) = A(pO+pfft, qO+qfft, rO+rfft);
    end
    end
    end
end
end
end 
end 
