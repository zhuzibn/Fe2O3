% atomistic model for alp_Fe2O3
% require nvidia GPU
clear all;clc;close all;tic
%% control parameter
constantfile;
clear gam
rk4=1;%1:rk4,0:heun Method,2:4th predictor-corrector
DMIenable=1;%enable DMI interaction?
dipolemode=0;%enable dipole field?
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
natomW=20;natomL=20;natomH=6;
systemgeneration();
%% FiM parameters
Ksim1=-8e-22;% 2ed-order
hk6=1.76e-27;% 6th-order
J1=1.3740e-21;
J2=-1.1696e-21;
J3=4.0413e-21;
J4=2.8041e-21;%[Phy. Rev.B.107,184426 (2023)]
g=2;%g-factor
gam=g*mub/(hbar*ele);%1/(s.T)
mus=4.2313*mub;
a=5.067e-10;
c=13.882e-10;
tz=c/6*natomH;
ms=12*mus/c/a^2/sqrt(3)*2;%[A/m], saturation magnetization
% Shell-resolved z-DMI magnitudes [J]; shells 1 and 2 are omitted.
% D3/D4=1.25 is a Fig.3(b) estimate; absolute values are user-selected.
if DMIenable
    Dsim3=1.6e-23;
    Dsim4=1.28e-23;
else
    Dsim3=0;
    Dsim4=0;
end
alp=0.001;%Gilbert damping
%% electrical parameters
jcSOT=6e12;%[A/m2]
jcSTT=0e9;%[A/m2]
Hext=[0,0,0];%
%% SOT parameters
SOT_DLT=1;%1(0),enable(disable) SOT damping torque
SOT_FLT=1;%1(0),enable(disable) SOT field-like torque
psj_theta=90/180*pi;
psj_phi=210/180*pi;
%spin flux polarization
psjSHE=[cos(psj_phi)*sin(psj_theta),sin(psj_phi)*sin(psj_theta),cos(psj_theta)];
psjSHEx=psjSHE(1);
psjSHEy=psjSHE(2);
psjSHEz=psjSHE(3);
thetaSH=0.2;%spin hall angle
if SOT_FLT
    chi=5;%ratio of FLT/DLT
else
    chi=0;
end
BDSOTRE=SOT_DLT*hbar/2*thetaSH*jcSOT/(ms*tz);%[T]
BDSOTTM=SOT_DLT*hbar/2*thetaSH*jcSOT/(ms*tz);
%% STT parameters
STT_DLT=0;%1(0),enable(disable) SOT damping torque
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
%% time control
gpusave=1e-12;%how often saving gpu data
gpurun_number=150;%run time /ps                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ;
tstep=1e-16;%one step is 0.1fs
savetstep=5000;%save data every ~ loops
gpusteps=round(gpusave/tstep);
runtime=gpurun_number*gpusave;
if ~mod(gpusteps,savetstep)==0
    error('gpusteps should be multiple integer times of savetstep, otherwise there might be errors')
end
tmp1=ones(1,gpusteps);
tmp2=tmp1(1:savetstep:end);
final_m_savestep=size(tmp2,2);
clear tmp1 tmp1
t=(0:gpurun_number*final_m_savestep-1)*savetstep*tstep;%TIME-01: elapsed time of each saved frame
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
clear mmx_l1 mmy_l1 mmz_l1
clear mmx_r1 mmy_r1 mmz_r1
clear mmx_l2 mmy_l2 mmz_l2
clear mmx_r2 mmy_r2 mmz_r2
clear mmx_u1 mmy_u1 mmz_u1
clear mmx_d1 mmy_d1 mmz_d1
clear mmx_u1l1 mmy_u1l1 mmz_u1l1
clear mmx_d1l1 mmy_d1l1 mmz_d1l1
clear mmx_d1r1 mmy_d1r1 mmz_d1r1
clear shift_matx_black shift_matx_blue shift_matx_gr shift_matx_or
clear shift_matx_purple shift_matx_red
clear shift_maty_black shift_maty_blue shift_maty_gr shift_maty_or
clear shift_maty_purple shift_maty_red
clear shift_matz_black shift_matz_blue shift_matz_gr shift_matz_or
clear shift_matz_purple shift_matz_red BDSTTTM BDSTTRE cosphi costheta 
clear atomtype_layer3red atomtype_layer3black atomtype_layer2p sxx syy szz
clear mmx_u1r1 mmy_u1r1 mmz_u1r1
clear atomtype_layer1or atomtype_layer1gr atomtype_layer2blue
clear mmxtmp mmytmp mmztmp
clear hdmi_x hdmi_y hdmi_z
save('final.mat')