function destination = run_regression(destination)
% Run all short cases against immutable accepted inputs, never create baselines.
here=fileparts(mfilename('fullpath')); root=fileparts(here);
oldpath=path; cleanup=onCleanup(@() path(oldpath)); addpath(root); addpath(here);
accepted_file=fullfile(here,'baseline_manifest.json');
assert(isfile(accepted_file),'Fe2O3:MissingBaseline','Transfer the accepted baseline first.');
accepted=jsondecode(fileread(accepted_file));
assert(strcmp(accepted.status,'accepted'),'Fe2O3:MissingBaseline','No accepted reference.');
reference=fullfile(here,accepted.directory);
% Verify the full suite and raw invariants before using any fixed restart input.
compare_regression(reference);
if nargin<1, destination=fe2o3_candidate_directory;
else, destination=fe2o3_candidate_directory(destination); end
cases=fe2o3_regression_cases;
manifest=struct('schema_version','1.0','status','candidate',...
    'cases',struct('id',{},'record_file',{},'record_sha256',{},'artifacts',{}));
for index=1:numel(cases)
    spec=cases(index); folder=fullfile(destination,spec.id); mkdir(folder);
    cfg=fe2o3_default_config;
    cfg.natomW=spec.core(1); cfg.natomL=spec.core(2); cfg.natomH=spec.core(3);
    cfg.rk4=spec.solver; cfg.load_fixed_atom_distrib=spec.restart; cfg.dipolemode=spec.dipole;
    cfg.jcSOT=spec.jcSOT; cfg.jcSTT=spec.jcSTT; cfg.STT_DLT=spec.STT_DLT; cfg.STT_FLT=spec.STT_FLT;
    cfg.tstep=spec.tstep; cfg.gpusave=spec.steps*spec.tstep;
    cfg.gpurun_number=spec.chunks; cfg.savetstep=spec.stride;
    if spec.restart, cfg.restart_file=fullfile(reference,spec.id,'Relax_0.mat'); end
    result=fe2o3_run(cfg,fullfile(folder,'result.mat'));
    record=legacy_record(result,spec);
    save(fullfile(folder,'record.mat'),'record');
    artifact=struct('file',[spec.id '/result.mat'],'sha256',fe2o3_file_hash(fullfile(folder,'result.mat')));
    manifest.cases(index)=struct('id',spec.id,'record_file',[spec.id '/record.mat'],...
        'record_sha256',fe2o3_file_hash(fullfile(folder,'record.mat')),'artifacts',artifact);
end
fid=fopen(fullfile(destination,'manifest.json'),'w'); assert(fid~=-1);
fprintf(fid,'%s\n',jsonencode(manifest,'PrettyPrint',true)); fclose(fid);
compare_regression(destination);
fprintf('Candidate suite: %s\n',destination);
end

function record=legacy_record(result,spec)
% Only adapt field names/legacy aliases, from actual new-run inputs and outputs.
cfg=result.configuration; d=result.derived;
names={'natomW','natomL','natomH','rk4','DMIenable','dipolemode',...
    'load_fixed_atom_distrib','save_fixed_atom_distrib','Ksim1','hk6','J1','J2','J3','J4',...
    'g','gam','mus','a','c','tz','ms','Dsim1','Dsim3','Dsim4','alp','jcSOT','jcSTT','Hext',...
    'SOT_DLT','SOT_FLT','psj_theta','psj_phi','psjSHE','thetaSH','chi','BDSOTRE','BDSOTTM',...
    'STT_DLT','STT_FLT','psjSTT','etaSTT','chiSTT','BDSTTRE','BDSTTTM',...
    'gpusave','gpurun_number','tstep','savetstep','gpusteps','runtime','final_m_savestep',...
    'ele','mub','hbar','mu_0','kb'};
configuration=struct;
for index=1:numel(names)
    name=names{index};
    if isfield(d,name), configuration.(name)=d.(name);
    else, configuration.(name)=cfg.(name); end
end
configuration.natomW=d.padded_dimensions(1); configuration.natomL=d.padded_dimensions(2);
record=struct('schema_version','1.0','specification',spec,'configuration',configuration,...
    'geometry',result.geometry,'initial_state',result.initial_state,'units',result.units,...
    'sampling',result.sampling,'trajectory',result.trajectory,'provenance',result.provenance);
end
