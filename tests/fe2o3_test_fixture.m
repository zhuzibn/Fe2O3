function p = fe2o3_test_fixture(core)
% Read actual production defaults, without integration, saving or UI effects.
root=fileparts(fileparts(mfilename('fullpath')));
oldpath=path; cleanup=onCleanup(@() path(oldpath));
addpath(root);
cfg=fe2o3_default_config;
cfg.natomW=core(1); cfg.natomL=core(2); cfg.natomH=core(3);
[cfg,d]=fe2o3_validate_config(cfg);
% A test-only script workspace; production uses explicit assignments.
for name=fieldnames(cfg)', eval([name{1} '=cfg.(name{1});']); end
for name=fieldnames(d)', eval([name{1} '=d.(name{1});']); end
psjSHEx=psjSHE(1); psjSHEy=psjSHE(2); psjSHEz=psjSHE(3);
psjSTTx=psjSTT(1); psjSTTy=psjSTT(2); psjSTTz=psjSTT(3);
systemgeneration;
names={'natomW','natomL','natomH','atomtype_',...
    'atomtype_layer1gr','atomtype_layer1or','atomtype_layer2blue',...
    'atomtype_layer2p','atomtype_layer3red','atomtype_layer3black',...
    'mx_init','my_init','mz_init','mu_0','gam','mus','alp',...
    'J1','J2','J3','J4','Ksim1','hk6','Dsim1','Dsim3','Dsim4',...
    'Hext','BDSOTRE','BDSTTRE','chi','psjSHEx','psjSHEy','psjSHEz',...
    'psjSTTx','psjSTTy','psjSTTz','rk4','dipolemode'};
p=struct;
for j=1:numel(names), p.(names{j})=eval(names{j}); end
p.root=root;
p.core=core;
end
