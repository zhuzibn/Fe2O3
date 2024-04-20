% atomistic model for alp_Fe2O3
% require nvidia GPU
clear all;clc;close all;tic
%% control parameter
constantfile;
clear gam
rk4=1;%1:rk4,0:heun Method,2:4th predictor-corrector
DMIenable=1;%enable DMI interaction?
dwcalc=1;%1:simulate dw motion 0: no domain wall
thermalenable=0;%enable thermal field?
dipolemode=1;%enable dipole field?
%% use fixed or randomatom distribution
load_fixed_atom_distrib=0;%load fixed atom distribution
save_fixed_atom_distrib=0;%save fixed atom distribution
if save_fixed_atom_distrib && load_fixed_atom_distrib
    error('only one of load_fixed_atom_distrib or save_fixed_atom_distrib can be enabled');
end
if load_fixed_atom_distrib
    display('use fixed atom distribution')
elseif save_fixed_atom_distrib
    display('save fixed atom distribution')
else
    display('use random atom distribution')
    rng('shuffle');
end
%% optional control
%gpuDevice(1)%select GPU device
%% system generation
natomW=100;natomL=100;natomH=6;
c=13.882e-10;% Phy. Rev.B.107,184426 (2023)
h=natomH*c;
systemgeneration();
%% FiM parameters
Ksim1=1.8073e-25;% second-order anisotropy
Ksim2=1.7624e-27;% forth-order anisotropy
J1=1.3740e-21;J2=-1.1696e-21;J3=4.0413e-21;J4=2.8041e-21;%[Phy. Rev.B.107,184426 (2023)]
g=2;%g-factor
gam=g*mub/(hbar*ele);%1/(s.T)
mus=4.2313*mub;
d=0.35*1e-9;%lattice constant
tz=h;
ms=mus/d^3;%[A/m], saturation magnetization
if DMIenable
    Dsim=-4.8065e-23;%[J], DMI
else
    Dsim=0;
end
alp=0.001;%Gilbert damping
%% electrical parameters
jcSOT=0e9;%[A/m2]
jcSTT=0e9;%[A/m2]
Hext=[0,0,0e-3];%
%% SOT parameters
SOT_DLT=0;%1(0),enable(disable) SOT damping torque
SOT_FLT=0;%1(0),enable(disable) SOT field-like torque
psjSHE=[0,1,0];%spin flux polarization
psjSHEx=psjSHE(1);
psjSHEy=psjSHE(2);
psjSHEz=psjSHE(3);
thetaSH=0.2;%spin hall angle
if SOT_FLT
    chi=0.2;%ratio of FLT/DLT
else
    chi=0;
end
BDSOTRE=SOT_DLT*hbar/2*thetaSH*jcSOT/(ms*tz);%[T]
BDSOTTM=SOT_DLT*hbar/2*thetaSH*jcSOT/(ms*tz);
%% STT parameters
STT_DLT=1;%1(0),enable(disable) SOT damping torque
STT_FLT=0;%1(0),enable(disable) SOT field-like torque
psjSTT=[0,0,1];%spin flux polarization
psjSTTx=psjSTT(1);
psjSTTy=psjSTT(2);
psjSTTz=psjSTT(3);
etaSTT=0.8;%spin hall angle
if STT_FLT
    chiSTT=0.2;%ratio of FLT/DLT
else
    chiSTT=0;
end
BDSTTRE=STT_DLT*hbar/2*etaSTT*jcSTT/(ms*tz);%[T]
BDSTTTM=STT_DLT*hbar/2*etaSTT*jcSTT/(ms*tz);
%% other parameters
T=100;%[K]
%% time control
gpusave=1e-12;%how often saving gpu data
gpurun_number=200;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ;
tstep=1e-16;
savetstep=1000;%this is used to reduce data size
gpusteps=round(gpusave/tstep);
runtime=gpurun_number*gpusave;%second run for dw motion
totstep=round(runtime/tstep);
t=linspace(tstep,runtime,totstep);%This need to be optimized

%% increase data (optional)
% saven=200;
% runtime=saven*gpusave;%runtime for this code

if ~mod(gpusteps,savetstep)==0
    error('gpusteps should be multiple integer times of savetstep, otherwise there might be errors')
end
tmp1=ones(1,gpusteps);
tmp2=tmp1(1:savetstep:end);
final_m_savestep=size(tmp2,2);
clear tmp1 tmp1
if (SOT_DLT || SOT_FLT) && ~(rk4==1)
    error('only rk4 is implemented for spin torque driven')
end
%% dynamic calc
integrate_llg(); toc
%% save data
clear hdipolex1 hdipoley1 hdipolez1  
clear hdipolex2 hdipoley2 hdipolez2  
clear hdipolex3 hdipoley3 hdipolez3  
clear hdipolex4 hdipoley4 hdipolez4  
clear hdipolex5 hdipoley5 hdipolez5  
clear hdipolex6 hdipoley6 hdipolez6  
clear mmx_0n1 mmy_0n1 mmz_0n1
clear mmx_0p1 mmy_0p1 mmz_0p1
clear mmx_0n2 mmy_0n2 mmz_0n2
clear mmx_0p2 mmy_0p2 mmz_0p2
clear mmx_n10 mmy_n10 mmz_n10
clear mmx_p10 mmy_p10 mmz_p10
clear mmx_n1n1 mmy_n1n1 mmz_n1n1
clear mmx_p1n1 mmy_p1n1 mmz_p1n1
clear mmx_p1p1 mmy_p1p1 mmz_p1p1
clear shift_matx_black shift_matx_blue shift_matx_gr shift_matx_or
clear shift_matx_purple shift_matx_red
clear shift_maty_black shift_maty_blue shift_maty_gr shift_maty_or
clear shift_maty_purple shift_maty_red
clear shift_matz_black shift_matz_blue shift_matz_gr shift_matz_or
clear shift_matz_purple shift_matz_red
save('finalnew.mat')