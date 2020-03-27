% Hamiltonian for a zigzag nanoribbon
% Input:
% % x, composition
% % n, strain parameter
% % f, flatness parameter in [0, 1] 
% % uc, # of unit cells
% % stagg, staggered potential
% % kpoints, # of kpoints along x axis
% Output:
% % k (vector), kpoints
% % V (matrix), eigenvalues for kpoints

function [k, V] = HamiltonianZigzag(x, n, f, uc, stagg, kpoints)
% 
tor = 0;
N_at = 4*uc + 4;
V = zeros(N_at*8, kpoints);

% % Geometrical parameters
dz0 = (1 - x)*1.59 + x*1.64; % unstrained perpendicular bond component
a0 = (1 - x)*4.5332 + x*4.3007; % unstrained lattice constant

dz = f*((1 - x)*1.59 + x*1.64);
a = (1 - x)*(-0.3*sqrt(3)*dz + 5.35) + x*(-0.24*sqrt(3)*dz + 5.00);
g = (1 - x)*1.3861 + x*1.4610;
dx = a/sqrt(3);
bondd = sqrt(dz^2 + dx^2);
scale = (bondd/sqrt((a0/sqrt(3))^2 + dz0^2))^(-n);
c = (1 - x)*11.7967 + x*(11.2221);
mu = (1 - x)*0.2341 + x*(0.2336);

a1 = [(-1/2)*a, (-sqrt(3)/6)*a, ((1 - x)*2.4731 + x*2.28091)* ...
    sqrt(bondd^2 - dx^2)];
a2 = [(1/2)*a, (-sqrt(3)/6)*a, ((1 - x)*2.4731 + x*2.28091)* ...
    sqrt(bondd^2 - dx^2)];
a3 = [0, (sqrt(3)/3)*a, ((1 - x)*2.4731 + x*2.28091)* ...
    sqrt(bondd^2 - dx^2)];
d = [0, 0, 2*mu]*3*((1 - x)*2.4731 + x*2.28091)*sqrt(bondd^2 - dx^2);
v_for_direc_cos = a2 - d;
alpha = v_for_direc_cos(1)./norm(v_for_direc_cos);
beta = v_for_direc_cos(2)./norm(v_for_direc_cos);
gamma = v_for_direc_cos(3)./norm(v_for_direc_cos);

% % Hopping terms
lambda = (1 - x)*1.5 + x*(0.6);
Es = (1 - x)*(-10.906) + x*(-10.068);
Ep = (1 - x)*(-0.486) + x*(-0.926);
Vss_sigma = ((1 - x)*(-0.608) + x*(-0.694))*scale;
Vsp_sigma = ((1 - x)*1.320 + x*(1.554))*scale;
Vpp_sigma = ((1 - x)*1.854 + x*(2.342))*scale;
Vpp_pi = ((1 - x)*(-0.600) + x*(-0.582))*scale;
Vss_sigma_p = ((1 - x)*0 + x*0)*scale;
Vsp_sigma_p = ((1 - x)*0 + x*0)*scale;
Vpp_sigma_p = ((1 - x)*0 + x*0)*scale;
Vpp_pi_p = ((1 - x)*0 + x*(0))*scale;
Vss_sigma_pp =((1 - x)*0 + x*(0))*scale;
Vsp_sigma_pp =((1 - x)*0 + x*(0))*scale;
Vpp_sigma_pp =((1 - x)*0.156 + x*(0.352))*scale;
Vpp_pi_pp = ((1 - x)*0 + x*0)*scale;

% % Hamiltonian: on-site
TYP0=zeros(8,8);
TYP0(1,1)=Es;
TYP0(2,2)=Es;
TYP0(3,3)=Ep;
TYP0(3,4)=-1i.*1/3*lambda;
TYP0(3,8)=1/3*lambda;
TYP0(4,4)=Ep;
TYP0(4,3)=1i.*1/3*lambda;
TYP0(4,8)=-1i.*1/3*lambda;
TYP0(5,5)=Ep;
TYP0(5,6)=-1/3*lambda;
TYP0(5,7)=1i.*1/3*lambda;
TYP0(6,6)=Ep;
TYP0(6,5)=-1/3*lambda;
TYP0(6,7)=1i.*1/3*lambda;
TYP0(7,7)=Ep;
TYP0(7,6)=-1i.*1/3*lambda;
TYP0(7,5)=-1i.*1/3*lambda;
TYP0(8,8)=Ep;
TYP0(8,3)=1/3*lambda;
TYP0(8,4)=1i.*1/3*lambda;

TYP1=zeros(8,8);
TYP1(1,1)=1.*Vss_sigma;
TYP1(1,3)=-alpha.*Vsp_sigma;
TYP1(1,4)=beta.*Vsp_sigma;
TYP1(1,5)=gamma.*Vsp_sigma;
TYP1(2,2)=1.*Vss_sigma;
TYP1(2,6)=-alpha.*Vsp_sigma;
TYP1(2,7)=beta.*Vsp_sigma;
TYP1(2,8)=gamma.*Vsp_sigma;
TYP1(3,1)=alpha.*Vsp_sigma;
TYP1(3,3)=(alpha.^2).*Vpp_sigma + (1-alpha.^2).*Vpp_pi;
TYP1(3,4)=(-alpha.*beta).*(Vpp_sigma-Vpp_pi);
TYP1(3,5)=(-alpha.*gamma).*(Vpp_sigma-Vpp_pi);
TYP1(4,1)=-beta.*Vsp_sigma;
TYP1(4,3)=(-alpha.*beta).*(Vpp_sigma-Vpp_pi);
TYP1(4,4)=(beta.^2).*Vpp_sigma + (1-beta.^2).*Vpp_pi;
TYP1(4,5)=(beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP1(5,1)=-gamma.*Vsp_sigma;
TYP1(5,3)=(-alpha.*gamma).*(Vpp_sigma-Vpp_pi);
TYP1(5,4)=(beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP1(5,5)=(gamma.^2).*Vpp_sigma + (1-gamma.^2).*Vpp_pi;
TYP1(6,2)=alpha.*Vsp_sigma;
TYP1(6,6)=(alpha.^2).*Vpp_sigma + (1-alpha.^2).*Vpp_pi;
TYP1(6,7)=(-alpha.*beta).*(Vpp_sigma-Vpp_pi);
TYP1(6,8)=(-alpha.*gamma).*(Vpp_sigma-Vpp_pi);
TYP1(7,2)=-beta.*Vsp_sigma;
TYP1(7,6)=(-alpha.*beta).*(Vpp_sigma-Vpp_pi);
TYP1(7,7)=(beta.^2).*Vpp_sigma + (1-beta.^2).*Vpp_pi;
TYP1(7,8)=(beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP1(8,2)=-gamma.*Vsp_sigma;
TYP1(8,6)=(-alpha.*gamma).*(Vpp_sigma-Vpp_pi);
TYP1(8,7)=(beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP1(8,8)=(gamma.^2).*Vpp_sigma + (1-gamma.^2).*Vpp_pi;

TYP2=zeros(8,8);
TYP2(1,1)=1.*Vss_sigma;
TYP2(1,3)=alpha.*Vsp_sigma;
TYP2(1,4)=beta.*Vsp_sigma;
TYP2(1,5)=gamma.*Vsp_sigma;
TYP2(2,2)=1.*Vss_sigma;
TYP2(2,6)=alpha.*Vsp_sigma;
TYP2(2,7)=beta.*Vsp_sigma;
TYP2(2,8)=gamma.*Vsp_sigma;
TYP2(3,1)=-alpha.*Vsp_sigma;
TYP2(3,3)=(alpha.^2).*Vpp_sigma + (1-alpha.^2).*Vpp_pi;
TYP2(3,4)=(alpha.*beta).*(Vpp_sigma-Vpp_pi);
TYP2(3,5)=(alpha.*gamma).*(Vpp_sigma-Vpp_pi);
TYP2(4,1)=-beta.*Vsp_sigma;
TYP2(4,3)=(alpha.*beta).*(Vpp_sigma-Vpp_pi);
TYP2(4,4)=(beta.^2).*Vpp_sigma + (1-beta.^2).*Vpp_pi;
TYP2(4,5)=(beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP2(5,1)=-gamma.*Vsp_sigma;
TYP2(5,3)=(alpha.*gamma).*(Vpp_sigma-Vpp_pi);
TYP2(5,4)=(beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP2(5,5)=(gamma.^2).*Vpp_sigma + (1-gamma.^2).*Vpp_pi;
TYP2(6,2)=-alpha*Vsp_sigma;
TYP2(6,6)=(alpha.^2).*Vpp_sigma + (1-alpha.^2).*Vpp_pi;
TYP2(6,7)=(alpha.*beta).*(Vpp_sigma-Vpp_pi);
TYP2(6,8)=(alpha.*gamma).*(Vpp_sigma-Vpp_pi);
TYP2(7,2)=-beta*Vsp_sigma;
TYP2(7,6)=(alpha.*beta).*(Vpp_sigma-Vpp_pi);
TYP2(7,7)=(beta.^2).*Vpp_sigma + (1-beta.^2).*Vpp_pi;
TYP2(7,8)=(beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP2(8,2)=-gamma*Vsp_sigma;
TYP2(8,6)=(alpha.*gamma).*(Vpp_sigma-Vpp_pi);
TYP2(8,7)=(beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP2(8,8)=(gamma.^2).*Vpp_sigma + (1-gamma.^2).*Vpp_pi;

TYP3=zeros(8,8);
TYP3(1,1)=1*Vss_sigma;
TYP3(1,3)=0*Vsp_sigma;
TYP3(1,4)=-2*beta*Vsp_sigma;
TYP3(1,5)=gamma*Vsp_sigma;
TYP3(2,2)=1*Vss_sigma;
TYP3(2,6)=0*Vsp_sigma;
TYP3(2,7)=-2*beta*Vsp_sigma;
TYP3(2,8)=gamma*Vsp_sigma;
TYP3(3,1)=0*Vsp_sigma;
TYP3(3,3)=Vpp_pi;
TYP3(3,4)=0;
TYP3(3,5)=0;
TYP3(4,1)=2*beta*Vsp_sigma;
TYP3(4,3)=0;
TYP3(4,4)=(4*beta.^2).*Vpp_sigma + (1-4*beta.^2).*Vpp_pi;
TYP3(4,5)=(-2*beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP3(5,1)=-gamma*Vsp_sigma;
TYP3(5,3)=0;
TYP3(5,4)=(-2*beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP3(5,5)=(gamma.^2).*Vpp_sigma + (1-gamma.^2).*Vpp_pi;
TYP3(6,2)=0*Vsp_sigma;
TYP3(6,6)=Vpp_pi;
TYP3(6,7)=0;
TYP3(6,8)=0;
TYP3(7,2)=2*beta*Vsp_sigma;
TYP3(7,6)=0;
TYP3(7,7)=(4*beta.^2).*Vpp_sigma + (1-4*beta.^2).*Vpp_pi;
TYP3(7,8)=(-2*beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP3(8,2)=-gamma*Vsp_sigma;
TYP3(8,6)=0;
TYP3(8,7)=(-2*beta.*gamma).*(Vpp_sigma-Vpp_pi);
TYP3(8,8)=(gamma.^2).*Vpp_sigma + (1-gamma.^2).*Vpp_pi;

TYP1N=zeros(8,8);
TYP1N(1,1)=Vss_sigma_pp;
TYP1N(1,3)=-Vsp_sigma_pp;
TYP1N(1,4)=0*Vsp_sigma_pp;
TYP1N(2,2)=Vss_sigma_pp;
TYP1N(2,6)=-Vsp_sigma_pp;
TYP1N(2,7)=0*Vsp_sigma_pp;
TYP1N(3,1)=Vsp_sigma_pp;
TYP1N(3,3)=1*Vpp_sigma_pp+0*Vpp_pi_pp;
TYP1N(3,4)=0*(Vpp_sigma_pp-Vpp_pi_pp);
TYP1N(4,1)=0*Vsp_sigma_pp;
TYP1N(4,3)=0*(Vpp_sigma_pp-Vpp_pi_pp);
TYP1N(4,4)=0*Vpp_sigma_pp+1*Vpp_pi_pp;
TYP1N(5,5)=Vpp_pi_pp;
TYP1N(6,2)=1*Vsp_sigma_pp;
TYP1N(6,6)= 1*Vpp_sigma_pp+0*Vpp_pi_pp;
TYP1N(6,7)=0*(Vpp_sigma_pp-Vpp_pi_pp);
TYP1N(7,2)=0*Vsp_sigma_pp;
TYP1N(7,6)=0*(Vpp_sigma_pp-Vpp_pi_pp);
TYP1N(7,7)=0*Vpp_sigma_pp+1*Vpp_pi_pp;
TYP1N(8,8)=Vpp_pi_pp;

TYP2N=zeros(8,8);
TYP2N(1,1)=Vss_sigma_pp;
TYP2N(1,3)=Vsp_sigma_pp;
TYP2N(1,4)=0*Vsp_sigma_pp;
TYP2N(2,2)=Vss_sigma_pp;
TYP2N(2,6)=Vsp_sigma_pp;
TYP2N(2,7)=0*Vsp_sigma_pp;
TYP2N(3,1)=-Vsp_sigma_pp;
TYP2N(3,3)=1*Vpp_sigma_pp+0*Vpp_pi_pp;
TYP2N(3,4)=0*(Vpp_sigma_pp-Vpp_pi_pp);
TYP2N(4,1)=0*Vsp_sigma_pp;
TYP2N(4,3)=0*(Vpp_sigma_pp-Vpp_pi_pp);
TYP2N(4,4)=0*Vpp_sigma_pp+1*Vpp_pi_pp;
TYP2N(5,5)=Vpp_pi_pp;
TYP2N(6,2)=-1*Vsp_sigma_pp;
TYP2N(6,6)= 1*Vpp_sigma_pp+0*Vpp_pi_pp;
TYP2N(6,7)=0*(Vpp_sigma_pp-Vpp_pi_pp);
TYP2N(7,2)=0*Vsp_sigma_pp;
TYP2N(7,6)=0*(Vpp_sigma_pp-Vpp_pi_pp);
TYP2N(7,7)=0*Vpp_sigma_pp+1*Vpp_pi_pp;
TYP2N(8,8)=Vpp_pi_pp;

TYP3N=zeros(8,8);
TYP3N(1,1)=Vss_sigma_pp;
TYP3N(1,3)=1/2*Vsp_sigma_pp;
TYP3N(1,4)=-sqrt(3)/2*Vsp_sigma_pp;
TYP3N(2,2)=Vss_sigma_pp;
TYP3N(2,6)=1/2*Vsp_sigma_pp;
TYP3N(2,7)=-sqrt(3)/2*Vsp_sigma_pp;
TYP3N(3,1)=-1/2*Vsp_sigma_pp;
TYP3N(3,3)=1/4*Vpp_sigma_pp+3/4*Vpp_pi_pp;
TYP3N(3,4)=-sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP3N(4,1)=sqrt(3)/2*Vsp_sigma_pp;
TYP3N(4,3)=-sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP3N(4,4)=3/4*Vpp_sigma_pp+1/4*Vpp_pi_pp;
TYP3N(5,5)=Vpp_pi_pp;
TYP3N(6,2)=-1/2*Vsp_sigma_pp;
TYP3N(6,6)= 1/4*Vpp_sigma_pp+3/4*Vpp_pi_pp;
TYP3N(6,7)=-sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP3N(7,2)=sqrt(3)/2*Vsp_sigma_pp;
TYP3N(7,6)=-sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP3N(7,7)=3/4*Vpp_sigma_pp+1/4*Vpp_pi_pp;
TYP3N(8,8)=Vpp_pi_pp;

TYP4N=zeros(8,8);
TYP4N(1,1)=Vss_sigma_pp;
TYP4N(1,3)=-1/2*Vsp_sigma_pp;
TYP4N(1,4)=sqrt(3)/2*Vsp_sigma_pp;
TYP4N(2,2)=Vss_sigma_pp;
TYP4N(2,6)=-1/2*Vsp_sigma_pp;
TYP4N(2,7)=sqrt(3)/2*Vsp_sigma_pp;
TYP4N(3,1)=1/2*Vsp_sigma_pp;
TYP4N(3,3)=1/4*Vpp_sigma_pp+3/4*Vpp_pi_pp;
TYP4N(3,4)=-sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP4N(4,1)=-sqrt(3)/2*Vsp_sigma_pp;
TYP4N(4,3)=-sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP4N(4,4)=3/4*Vpp_sigma_pp+1/4*Vpp_pi_pp;
TYP4N(5,5)=Vpp_pi_pp;
TYP4N(6,2)=1/2*Vsp_sigma_pp;
TYP4N(6,6)= 1/4*Vpp_sigma_pp+3/4*Vpp_pi_pp;
TYP4N(6,7)=-sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP4N(7,2)=-sqrt(3)/2*Vsp_sigma_pp;
TYP4N(7,6)=-sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP4N(7,7)=3/4*Vpp_sigma_pp+1/4*Vpp_pi_pp;
TYP4N(8,8)=Vpp_pi_pp;

TYP5N=zeros(8,8);
TYP5N(1,1)=Vss_sigma_pp;
TYP5N(1,3)=-1/2*Vsp_sigma_pp;
TYP5N(1,4)=-sqrt(3)/2*Vsp_sigma_pp;
TYP5N(2,2)=Vss_sigma_pp;
TYP5N(2,6)=-1/2*Vsp_sigma_pp;
TYP5N(2,7)=-sqrt(3)/2*Vsp_sigma_pp;
TYP5N(3,1)=1/2*Vsp_sigma_pp;
TYP5N(3,3)=1/4*Vpp_sigma_pp+3/4*Vpp_pi_pp;
TYP5N(3,4)=sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP5N(4,1)=sqrt(3)/2*Vsp_sigma_pp;
TYP5N(4,3)=sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP5N(4,4)=3/4*Vpp_sigma_pp+1/4*Vpp_pi_pp;
TYP5N(5,5)=Vpp_pi_pp;
TYP5N(6,2)=1/2*Vsp_sigma_pp;
TYP5N(6,6)= 1/4*Vpp_sigma_pp+3/4*Vpp_pi_pp;
TYP5N(6,7)=sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP5N(7,2)=sqrt(3)/2*Vsp_sigma_pp;
TYP5N(7,6)=sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP5N(7,7)=3/4*Vpp_sigma_pp+1/4*Vpp_pi_pp;
TYP5N(8,8)=Vpp_pi_pp;

TYP6N=zeros(8,8);
TYP6N(1,1)=Vss_sigma_pp;
TYP6N(1,3)=1/2*Vsp_sigma_pp;
TYP6N(1,4)=sqrt(3)/2*Vsp_sigma_pp;
TYP6N(2,2)=Vss_sigma_pp;
TYP6N(2,6)=1/2*Vsp_sigma_pp;
TYP6N(2,7)=sqrt(3)/2*Vsp_sigma_pp;
TYP6N(3,1)=-1/2*Vsp_sigma_pp;
TYP6N(3,3)=1/4*Vpp_sigma_pp+3/4*Vpp_pi_pp;
TYP6N(3,4)=sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP6N(4,1)=-sqrt(3)/2*Vsp_sigma_pp;
TYP6N(4,3)=sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP6N(4,4)=3/4*Vpp_sigma_pp+1/4*Vpp_pi_pp;
TYP6N(5,5)=Vpp_pi_pp;
TYP6N(6,2)=-1/2*Vsp_sigma_pp;
TYP6N(6,6)= 1/4*Vpp_sigma_pp+3/4*Vpp_pi_pp;
TYP6N(6,7)=sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP6N(7,2)=-sqrt(3)/2*Vsp_sigma_pp;
TYP6N(7,6)=sqrt(3)/4*(Vpp_sigma_pp-Vpp_pi_pp);
TYP6N(7,7)=3/4*Vpp_sigma_pp+1/4*Vpp_pi_pp;
TYP6N(8,8)=Vpp_pi_pp;

k = linspace(-pi/a, pi/a, kpoints);
for i = 1:kpoints
    H_NN_TYP0=eye(N_at,N_at);
    [s1,s2]=size(TYP0);
    H_NN_TYP0_full=zeros(s1*N_at,s1*N_at);
    H_NN_TYP0_full=kron(H_NN_TYP0,TYP0);

    H_NN_TYP1=zeros(N_at,N_at);
    H_NN_TYP1(2,1)=1;
    H_NN_TYP1(4,3)=1.*exp(-1i*k(i)*a);
    for j=1:uc
      H_NN_TYP1(6+4*(j-1),5+4*(j-1))=1;
      H_NN_TYP1(8+4*(j-1),7+4*(j-1))=1.*exp(-1i*k(i)*a);
    end
    H_NN_TYP1_full=zeros(s1*N_at,s1*N_at);
    H_NN_TYP1_full=kron(H_NN_TYP1,TYP1)+kron(H_NN_TYP1',TYP1');
    
    H_NN_TYP2=zeros(N_at,N_at);
    H_NN_TYP2(2,1)=1.*exp(1i*k(i)*a);
    H_NN_TYP2(4,3)=1;
    for j=1:uc
        H_NN_TYP2(6+4*(j-1),5+4*(j-1))=1.*exp(1i*k(i)*a);
        H_NN_TYP2(8+4*(j-1),7+4*(j-1))=1;
    end
    H_NN_TYP2_full=zeros(s1*N_at,s1*N_at);
    H_NN_TYP2_full=kron(H_NN_TYP2,TYP2)+kron(H_NN_TYP2',TYP2');

    H_NN_TYP3=zeros(N_at,N_at);
    H_NN_TYP3(2,3)=1;
    for j=1:uc
        H_NN_TYP3(4+4*(j-1),5+4*(j-1))=1;
        H_NN_TYP3(6+4*(j-1),7+4*(j-1))=1;
    end
    H_NN_TYP3(8+(j-1)*4,1)=1.*tor;
    H_NN_TYP3_full=zeros(s1*N_at,s1*N_at);
    H_NN_TYP3_full=kron(H_NN_TYP3,TYP3)+kron(H_NN_TYP3',TYP3');
    H_NN=H_NN_TYP0_full+H_NN_TYP1_full+H_NN_TYP2_full+H_NN_TYP3_full;

% %Hamiltonian: NNN
    H_NNN_TYP1N=eye(N_at,N_at).*exp(-1i*k(i)*a);
    H_NNN_TYP1N_full=zeros(s1*N_at,s1*N_at);
    H_NNN_TYP1N_full=kron(H_NNN_TYP1N,TYP1N);
    H_NNN_TYP2N=eye(N_at,N_at).*exp(1i*k(i)*a);
    H_NNN_TYP2N_full=zeros(s1*N_at,s1*N_at);
    H_NNN_TYP2N_full=kron(H_NNN_TYP2N,TYP2N);
    H_NNN_TYP3N=zeros(N_at,N_at);
    H_NNN_TYP3N(3,1)=1.*exp(1i*k(i)*a);
    H_NNN_TYP3N(4,2)=1;
    for j=1:uc
        H_NNN_TYP3N(5+4*(j-1),3+4*(j-1))=1;
        H_NNN_TYP3N(6+4*(j-1),4+4*(j-1))=1.*exp(1i*k(i)*a);
        H_NNN_TYP3N(7+4*(j-1),5+4*(j-1))=1.*exp(1i*k(i)*a);
        H_NNN_TYP3N(8+4*(j-1),6+4*(j-1))=1;
    end
    H_NNN_TYP3N(1,7+(j-1)*4)= 1.*tor;
    H_NNN_TYP3N(2,8+(j-1)*4)=1.*exp(1i*k(i)*a).*tor ;
    H_NNN_TYP3N_full=zeros(s1*N_at,s1*N_at);
    H_NNN_TYP3N_full=kron(H_NNN_TYP3N,TYP3N);
    
    H_NNN_TYP4N=zeros(N_at,N_at);
    H_NNN_TYP4N(1,3)=1.*exp(-1i*k(i)*a);
    H_NNN_TYP4N(2,4)=1;
    for j=1:uc
        H_NNN_TYP4N(3+4*(j-1),5+4*(j-1))=1;
        H_NNN_TYP4N(4+4*(j-1),6+4*(j-1))=1.*exp(-1i*k(i)*a);
        H_NNN_TYP4N(5+4*(j-1),7+4*(j-1))=1.*exp(-1i*k(i)*a);
        H_NNN_TYP4N(6+4*(j-1),8+4*(j-1))=1;
    end
    H_NNN_TYP4N(8+(j-1)*4,2)=1.*exp(-1i*k(i)*a).*tor;
    H_NNN_TYP4N(7+(j-1)*4,1)=1.*tor;
    H_NNN_TYP4N_full=zeros(s1*N_at,s1*N_at);
    H_NNN_TYP4N_full=kron(H_NNN_TYP4N,TYP4N);

    H_NNN_TYP5N=zeros(N_at,N_at);
    H_NNN_TYP5N(3,1)=1;
    H_NNN_TYP5N(4,2)=1.*exp(-1i*k(i)*a);
    for j=1:uc
        H_NNN_TYP5N(5+4*(j-1),3+4*(j-1))=1.*exp(-1i*k(i)*a);
        H_NNN_TYP5N(6+4*(j-1),4+4*(j-1))=1;
        H_NNN_TYP5N(7+4*(j-1),5+4*(j-1))=1;
        H_NNN_TYP5N(8+4*(j-1),6+4*(j-1))=1.*exp(-1i*k(i)*a);
    end
    H_NNN_TYP5N(1,7+(j-1)*4)=1.*exp(-1i*k(i)*a).*tor;
    H_NNN_TYP5N(2,8+(j-1)*4)=1.*tor;
    H_NNN_TYP5N_full=zeros(s1*N_at,s1*N_at);
    H_NNN_TYP5N_full=kron(H_NNN_TYP5N,TYP5N);


    H_NNN_TYP6N=zeros(N_at,N_at);
    H_NNN_TYP6N(1,3)=1;
    H_NNN_TYP6N(2,4)=1.*exp(1i*k(i)*a);
    for j=1:uc
        H_NNN_TYP6N(3+4*(j-1),5+4*(j-1))=1.*exp(1i*k(i)*a);
        H_NNN_TYP6N(4+4*(j-1),6+4*(j-1))=1;
        H_NNN_TYP6N(5+4*(j-1),7+4*(j-1))=1;
        H_NNN_TYP6N(6+4*(j-1),8+4*(j-1))=1.*exp(1i*k(i)*a);
    end
    H_NNN_TYP6N(8+(j-1)*4,2)=1.*tor;
    H_NNN_TYP6N(7+(j-1)*4,1)=1.*exp(1i*k(i)*a).*tor;
    H_NNN_TYP6N_full=zeros(s1*N_at,s1*N_at);
    H_NNN_TYP6N_full=kron(H_NNN_TYP6N,TYP6N);

    H_NNN=H_NNN_TYP1N_full+...
                     H_NNN_TYP2N_full+...
                     H_NNN_TYP3N_full+...
                     H_NNN_TYP4N_full+...
                     H_NNN_TYP5N_full+...
                     H_NNN_TYP6N_full;
    [U, V(:,i)] = eig(H_NN + H_NNN, 'vector');
end    