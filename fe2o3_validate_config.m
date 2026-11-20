function [cfg,d] = fe2o3_validate_config(cfg)
% Validate a complete input struct and derive the retained legacy arithmetic.
defaults=fe2o3_default_config;
assert(isstruct(cfg) && isscalar(cfg) && ...
    isequal(sort(fieldnames(cfg)),sort(fieldnames(defaults))),...
    'Fe2O3:InvalidConfig','Supply exactly the fields from fe2o3_default_config.');
names=fieldnames(defaults);
for index=1:numel(names)
    name=names{index}; value=cfg.(name);
    if strcmp(name,'restart_file')
        assert((ischar(value) && isrow(value)) || (isstring(value) && isscalar(value) && ...
            ~ismissing(value) && strlength(value)>0),...
            'Fe2O3:InvalidConfig','restart_file must be a nonempty text path.');
        cfg.(name)=char(value);
    else
        assert((isnumeric(value) || islogical(value)) && ~issparse(value) && ...
            ~isa(value,'gpuArray') && isreal(value) && ...
            isequal(size(value),size(defaults.(name))) && all(isfinite(value),'all'),...
            'Fe2O3:InvalidConfig','Invalid finite real input: %s.',name);
        cfg.(name)=double(value); % Integer/single inputs must not alter arithmetic precision.
    end
end
flags={'DMIenable','dipolemode','load_fixed_atom_distrib','save_fixed_atom_distrib',...
    'SOT_DLT','SOT_FLT','STT_DLT','STT_FLT'};
for name=flags
    assert(ismember(cfg.(name{1}),[0,1]),'Fe2O3:InvalidConfig','%s must be 0 or 1.',name{1});
end
assert(ismember(cfg.rk4,[0,1,2]),'Fe2O3:InvalidSolver','rk4 must be 0, 1 or 2.');
assert(~(cfg.load_fixed_atom_distrib && cfg.save_fixed_atom_distrib),...
    'Fe2O3:InvalidConfig','Distribution flags are mutually exclusive.');
for name={'natomW','natomL','natomH','gpurun_number','savetstep'}
    value=cfg.(name{1});
    assert(value>=1 && value==fix(value) && value<=flintmax,...
        'Fe2O3:InvalidConfig','%s must be a positive exactly represented integer.',name{1});
end
assert(mod(cfg.natomH,3)==0,'Fe2O3:InvalidConfig','natomH must be a positive multiple of 3.');
for name={'g','mus','a','c','gpusave','tstep'}
    assert(cfg.(name{1})>0,'Fe2O3:InvalidConfig','%s must be positive.',name{1});
end
assert(cfg.alp>=0,'Fe2O3:InvalidConfig','alp must be nonnegative.');
assert(abs(norm(cfg.psjSTT)-1)<=1e-12,'Fe2O3:InvalidConfig','psjSTT must be a unit row vector.');
constantfile;
d=struct('ele',ele,'mub',mub,'hbar',hbar,'mu_0',mu_0,'kb',kb);
d.gam=cfg.g*mub/(hbar*ele);
d.tz=cfg.c/6*cfg.natomH;
d.ms=12*cfg.mus/cfg.c/cfg.a^2/sqrt(3)*2;
d.Dsim1=0; d.Dsim3=cfg.Dsim3*cfg.DMIenable; d.Dsim4=cfg.Dsim4*cfg.DMIenable;
d.psjSHE=[cos(cfg.psj_phi)*sin(cfg.psj_theta),sin(cfg.psj_phi)*sin(cfg.psj_theta),cos(cfg.psj_theta)];
d.chi=cfg.chi*cfg.SOT_FLT;
d.chiSTT=0.2*cfg.STT_FLT; % Legacy diagnostic; BFSTT still uses chi below.
d.BDSOTRE=cfg.SOT_DLT*hbar/2*cfg.thetaSH*cfg.jcSOT/(d.ms*d.tz);
d.BDSOTTM=cfg.SOT_DLT*hbar/2*cfg.thetaSH*cfg.jcSOT/(d.ms*d.tz);
d.BDSTTRE=cfg.STT_DLT*hbar/2*cfg.etaSTT*cfg.jcSTT/(d.ms*d.tz);
d.BDSTTTM=cfg.STT_DLT*hbar/2*cfg.etaSTT*cfg.jcSTT/(d.ms*d.tz);
d.BFSOT=d.chi*d.BDSOTRE; d.BFSTT=d.chi*d.BDSTTRE;
ratio=cfg.gpusave/cfg.tstep;
d.gpusteps=round(ratio);
assert(isfinite(ratio) && ratio<=flintmax && d.gpusteps>=1 && ...
    abs(ratio-d.gpusteps)<=32*eps(max(1,abs(ratio))) && mod(d.gpusteps,cfg.savetstep)==0,...
    'Fe2O3:InvalidConfig','gpusave/tstep must be integral and divisible by savetstep.');
d.runtime=cfg.gpurun_number*cfg.gpusave;
d.final_m_savestep=d.gpusteps/cfg.savetstep;
d.sample_count=cfg.gpurun_number*d.final_m_savestep;
d.core_dimensions=[cfg.natomW,cfg.natomL,cfg.natomH];
d.padded_dimensions=d.core_dimensions+[4,4,0];
for name=fieldnames(d)'
    assert(all(isfinite(d.(name{1})),'all'),'Fe2O3:InvalidConfig','Nonfinite derived %s.',name{1});
end
assert(d.ms>0 && d.tz>0 && d.gam>0 && d.runtime>0 && d.sample_count<=flintmax && ...
    cfg.gpurun_number*d.gpusteps<=flintmax,...
    'Fe2O3:InvalidConfig','Derived scales/counts are not representable.');
assert(cfg.rk4==1 || ~any([d.BDSOTRE,d.BFSOT,d.BDSTTRE,d.BFSTT]~=0),...
    'Fe2O3:UnsupportedTorque','SOT/STT torques require rk4=1.');
d.t=(0:d.sample_count-1)*cfg.savetstep*cfg.tstep;
end
