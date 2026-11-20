function fe2o3_result = fe2o3_run(cfg,output_file,overwrite)
% Run the retained GPU scripts in local scope. Empty output_file means no save.
if nargin<2, output_file=''; end
if nargin<3, overwrite=false; end
started=tic;
[cfg,derived]=fe2o3_validate_config(cfg);
restart_input='';
if cfg.load_fixed_atom_distrib, restart_input=cfg.restart_file; end
destination=fe2o3_output_path(output_file,restart_input,overwrite);
root=fileparts(mfilename('fullpath'));
sources=source_identity(root);
rng_before=rng;
rk4=cfg.rk4;
DMIenable=cfg.DMIenable;
dipolemode=cfg.dipolemode;
load_fixed_atom_distrib=cfg.load_fixed_atom_distrib;
save_fixed_atom_distrib=cfg.save_fixed_atom_distrib;
restart_file=cfg.restart_file;
natomW=cfg.natomW;
natomL=cfg.natomL;
natomH=cfg.natomH;
Ksim1=cfg.Ksim1;
hk6=cfg.hk6;
J1=cfg.J1;
J2=cfg.J2;
J3=cfg.J3;
J4=cfg.J4;
g=cfg.g;
mus=cfg.mus;
a=cfg.a;
c=cfg.c;
alp=cfg.alp;
jcSOT=cfg.jcSOT;
jcSTT=cfg.jcSTT;
Hext=cfg.Hext;
SOT_DLT=cfg.SOT_DLT;
SOT_FLT=cfg.SOT_FLT;
psj_theta=cfg.psj_theta;
psj_phi=cfg.psj_phi;
thetaSH=cfg.thetaSH;
STT_DLT=cfg.STT_DLT;
STT_FLT=cfg.STT_FLT;
psjSTT=cfg.psjSTT;
etaSTT=cfg.etaSTT;
gpusave=cfg.gpusave;
gpurun_number=cfg.gpurun_number;
tstep=cfg.tstep;
savetstep=cfg.savetstep;
ele=derived.ele;
mub=derived.mub;
hbar=derived.hbar;
mu_0=derived.mu_0;
kb=derived.kb;
gam=derived.gam;
tz=derived.tz;
ms=derived.ms;
Dsim1=derived.Dsim1;
Dsim3=derived.Dsim3;
Dsim4=derived.Dsim4;
psjSHE=derived.psjSHE;
chi=derived.chi;
chiSTT=derived.chiSTT;
BDSOTRE=derived.BDSOTRE;
BDSOTTM=derived.BDSOTTM;
BDSTTRE=derived.BDSTTRE;
BDSTTTM=derived.BDSTTTM;
gpusteps=derived.gpusteps;
runtime=derived.runtime;
final_m_savestep=derived.final_m_savestep;
t=derived.t;
psjSHEx=psjSHE(1); psjSHEy=psjSHE(2); psjSHEz=psjSHE(3);
psjSTTx=psjSTT(1); psjSTTy=psjSTT(2); psjSTTz=psjSTT(3);
if ~load_fixed_atom_distrib && ~save_fixed_atom_distrib, rng('shuffle'); end
restart_info=struct('format','none');
systemgeneration;
rng_initial=rng;
geometry=struct('atomtype_',atomtype_,...
    'atomtype_layer1gr',atomtype_layer1gr,'atomtype_layer1or',atomtype_layer1or,...
    'atomtype_layer2blue',atomtype_layer2blue,'atomtype_layer2p',atomtype_layer2p,...
    'atomtype_layer3red',atomtype_layer3red,'atomtype_layer3black',atomtype_layer3black,...
    'array_order','W,L,H,time');
initial_state=struct('mx',mx_init,'my',my_init,'mz',mz_init);
device=gpuDevice;
environment=struct('matlab_version',version,'platform',computer,'os',system_dependent('getos'),...
    'parallel_toolbox',ver('parallel'),'gpu_name',device.Name,...
    'compute_capability',device.ComputeCapability,'driver',device.GraphicsDriverVersion,...
    'toolkit',device.ToolkitVersion,'total_memory',device.TotalMemory);
integrate_llg;
assert(isequal(sources,source_identity(root)),...
    'Fe2O3:SourceChanged','Project sources or Git state changed during the run.');
trajectory=struct('mmx',mmx,'mmy',mmy,'mmz',mmz,...
    'hdmi_allx',hdmi_allx,'hdmi_ally',hdmi_ally,'hdmi_allz',hdmi_allz,'t',t);
fe2o3_result=struct('schema_version','1.0','configuration',cfg,'derived',derived,...
    'geometry',geometry,'initial_state',initial_state,'trajectory',trajectory,...
    'units',struct('spin','dimensionless','field','T','time','s'),...
    'sampling','pre-update; final endpoint excluded');
fe2o3_result.model=struct('classification','Class 4; material-specific atomistic fixed-length LLG/LLGS',...
    'boundary','empty lateral padding and finite end layers; not verified periodic',...
    'mapping','six software site masks; crystallographic mapping unresolved; array layer direction reverses physical z',...
    'dipole_distance_m',0.5067e-9,'dipole_vector_components_m',[0.5067e-9,0.25335e-9,0.4388e-9],...
    'dipole_refresh','chunk-local updates 0,50,100,...; cached within RK4 stages',...
    'solver','1=RK4 with fresh exchange/anisotropy/DMI stages; 0=legacy Heun; 2=legacy PC4 with three RK4 bootstrap updates');
fe2o3_result.limitations={...
    'Material parameter provenance, complete crystallographic mapping and directed DMI bonds remain unresolved.',...
    'The full hk6 energy mapping and dipole truncation justification remain unresolved.',...
    'Heun/PC4 are experimental fixed-field paths; nonzero torque requires RK4.',...
    'DLT flags also gate FLT amplitudes; BFSTT uses SOT chi, not the diagnostic chiSTT.',...
    'save_fixed_atom_distrib is a retained no-op; no thermal field is implemented.',...
    'Restart uses the last saved spin, not the final integrated endpoint or solver/cache history.'};
fe2o3_result.provenance=struct(...
    'created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')),...
    'entry_point','fe2o3_run','output_file',destination,'overwrite',logical(overwrite),...
    'source',sources,'environment',environment,'rng_before',rng_before,...
    'rng_initial',rng_initial,'rng_after',rng,'restart',restart_info,'elapsed_seconds',toc(started));
if ~isempty(destination)
    fe2o3_save_result(fe2o3_result,destination,restart_input,overwrite);
end
end

function info=source_identity(root)
% Check resolution, including untracked helpers; Git identity alone is insufficient.
names={'constantfile','distrib','systemgeneration','exchangej_','dmi','dipole_',...
    'field_calc','integrate_llg','atomgpurk4','atomgpu','atomgpupc4',...
    'fe2o3_default_config','fe2o3_validate_config','fe2o3_run',...
    'fe2o3_load_restart','fe2o3_sha256','fe2o3_output_path','fe2o3_save_result'};
for index=1:numel(names)
    actual=which(names{index}); expected=fullfile(root,[names{index} '.m']);
    actual=char(java.io.File(actual).getCanonicalPath());
    expected=char(java.io.File(expected).getCanonicalPath());
    if ispc, actual=lower(actual); expected=lower(expected); end
    assert(strcmp(actual,expected),...
        'Fe2O3:ShadowedSource','Unexpected source resolution: %s -> %s',names{index},actual);
end
files=dir(fullfile(root,'*.m'));
hashes=struct('file',{},'sha256',{});
for index=1:numel(files)
    hashes(index)=struct('file',files(index).name,'sha256',fe2o3_sha256(fullfile(root,files(index).name)));
end
info=struct('project_sources',hashes,'git_head','','git_branch','','git_status','');
commands={'rev-parse HEAD','branch --show-current','status --porcelain=v1 --untracked-files=normal'};
fields={'git_head','git_branch','git_status'};
for index=1:numel(commands)
    [status,value]=system(sprintf('git -C "%s" %s',root,commands{index}));
    assert(status==0,'Fe2O3:ProvenanceUnavailable','Cannot read Git provenance: %s',commands{index});
    info.(fields{index})=strtrim(value);
end
end

