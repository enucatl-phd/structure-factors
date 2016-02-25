clear;
clear all;
addpath('../functions');

%% Physical parameters
re= 2.8179403267e-15;
E=25;
L=12e-2;
p2=2e-6;
R=1.86e-6/2;
fv=0.4;
%
lambda=12.4/E*1e-10;
k=2*pi/lambda;
epsilon=lambda/p2*L;
%
Vs=4*pi/3*R^3;
% Contrast
% Silica & glycerol
[delta_mat(1),beta_mat(1)]= HE_SiO2_delta_beta_2pt0(E);
[delta_mat(2),beta_mat(2)]= HE_Glycerol_delta_beta(E);
%
delta=delta_mat(1)-delta_mat(2);
beta=beta_mat(1)-beta_mat(2);
%
Delta_chi=delta+j*beta;
%Delta_rho=k^2*Delta_chi/2/pi/re;

%% Sampling

N=512*2*2;
Range=4*2*epsilon;
Delta_s=Range/N;
x=(-N/2:N/2-1)*Delta_s;
f=(-N/2:N/2-1)/(N*Delta_s);
[X,Y]=meshgrid(x);
[fx,fy]=meshgrid(f);
clear f x;
R_s=sqrt(X.^2+Y.^2);
F=sqrt(fx.^2+fy.^2);
Q=2*pi*F;

%% Form Factor

PQ=(SphereFormFactor2(Q,R)).^2;
SQ=HardSphereStructureFactor2(Q,R,fv);

%I=PQ;
%I=fv*Delta_chi*conj(Delta_chi)*Vs*PQ;
I=fv*Delta_chi*conj(Delta_chi)*Vs*PQ.*SQ;
%figure(1);imagesc(log10(I));colormap('gray');
%figure(2);plot(log10(I(N/2+1,:)));

%% Real space autocorrelation

autocorrelation=k^2*ift2(I,1/(N*Delta_s));

figure(3);imagesc(autocorrelation);
figure(4);plot(X(N/2+1,:),autocorrelation(N/2+1,:));

%% DFEC from autocorrelation

test=ifftshift(autocorrelation);
sampling=ifftshift(X);
position=sampling(1,:);
corr_values=test(1,:);
clear sampling test;

autocorrelation_1=interp1(position,corr_values,epsilon)
mu_d_exp=corr_values(1)-autocorrelation_1

%% DFEC from the S.K. Lynch calculation

mu_d_SKLynch=DFEC_sphere(2*R,lambda,epsilon,fv,1,Delta_chi)
