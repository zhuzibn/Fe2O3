function run_interface_checks
% CPU input/ownership failures plus one small production GPU run and consumers.
root=fileparts(fileparts(mfilename('fullpath')));
oldpath=path; oldcwd=pwd; saved_rng=rng;
cleanup=onCleanup(@() restore_environment(oldpath,oldcwd,saved_rng));
addpath(root); addpath(fullfile(root,'benchmark'));
check_output_manifest;
scratch=fe2o3_candidate_directory; diary(fullfile(scratch,'interface.log'));
diary_cleanup=onCleanup(@() diary('off'));
cfg=fe2o3_default_config; [normalized,d]=fe2o3_validate_config(cfg);
assert(isequal(cfg,normalized) && d.sample_count==50 && d.runtime==25e-12);
assert(isequal(d.t,(0:49)*5000*1e-16));
typed=cfg; typed.natomH=int32(6); typed.natomW=single(20); typed.DMIenable=true;
[typed,td]=fe2o3_validate_config(typed); assert(isequal(typed,cfg) && isequal(td,d));
bad=cfg; bad.unknown=1; reject(@() fe2o3_run(bad),'Fe2O3:InvalidConfig');
bad=rmfield(cfg,'J1'); reject(@() fe2o3_run(bad),'Fe2O3:InvalidConfig');
bad_values={'natomW',0;'natomL',1.5;'natomH',4;'gpurun_number',0;...
    'savetstep',3;'tstep',0;'gpusave',1.25e-16;'mus',0;'g',-1;'a',0;'c',0;...
    'alp',-0.01;'J1',NaN;'hk6',Inf;'Hext',[0;0;0];'jcSOT',1i;...
    'DMIenable',2;'psjSTT',[0,0,2];'restart_file','';'restart_file',"";'chiSTT',1};
for index=1:size(bad_values,1)
    bad=cfg; bad.(bad_values{index,1})=bad_values{index,2};
    reject(@() fe2o3_run(bad),'Fe2O3:InvalidConfig');
end
bad=cfg; bad.rk4=3; reject(@() fe2o3_run(bad),'Fe2O3:InvalidSolver');
bad=cfg; bad.save_fixed_atom_distrib=1; bad.load_fixed_atom_distrib=1;
reject(@() fe2o3_run(bad),'Fe2O3:InvalidConfig');
for solver=[0,2]
    bad=cfg; bad.rk4=solver; bad.jcSOT=1;
    reject(@() fe2o3_run(bad),'Fe2O3:UnsupportedTorque');
end
off=cfg; off.DMIenable=0; off.SOT_DLT=0; off.jcSOT=1; off.rk4=0;
[~,od]=fe2o3_validate_config(off); assert(od.Dsim3==0 && od.Dsim4==0 && od.BFSOT==0);
driven=cfg; driven.jcSTT=1e9; driven.STT_DLT=1; driven.STT_FLT=1;
[~,drive]=fe2o3_validate_config(driven);
assert(drive.chiSTT==0.2 && drive.BFSTT==drive.chi*drive.BDSTTRE);
fprintf('PASS canonical defaults, numeric classes, invalid input and retained torque gating\n');

output=fullfile(scratch,'result.mat');
cfg.natomW=5; cfg.natomL=9; cfg.natomH=6;
cfg.tstep=1e-16; cfg.gpusave=4e-16; cfg.gpurun_number=2; cfg.savetstep=2;
J1=123; Hext=[99,99,99]; %#ok<NASGU> Caller variables cannot enter the runner.
result=fe2o3_run(cfg,output);
loaded=load(output); assert(isequal(fieldnames(loaded),{'fe2o3_result'}));
assert(isequaln(result,loaded.fe2o3_result));
assert(isequal(result.configuration,cfg) && result.derived.sample_count==4);
assert(~isa(result.trajectory.mmx,'gpuArray') && ~isa(result.geometry.atomtype_,'gpuArray'));
assert(all(isfinite(result.trajectory.mmx),'all'));
sources=result.provenance.source;
assert(~isempty(sources.git_head));
[git_status,actual_status]=system(sprintf('git -C "%s" status --porcelain=v1 --untracked-files=normal',root));
assert(git_status==0 && strcmp(sources.git_status,strtrim(actual_status)));
assert(any(strcmp({sources.project_sources.file},'fe2o3_run.m')));
for item=sources.project_sources
    assert(strcmp(item.sha256,fe2o3_sha256(fullfile(root,item.file))));
end
memory=fe2o3_run(cfg); assert(isequal(memory.trajectory,result.trajectory));
assert(isempty(memory.provenance.output_file));
fprintf('PASS in-memory isolation, versioned save/load, CPU arrays and actual dirty-source hashes\n');

before=fe2o3_sha256(output);
reject(@() fe2o3_run(cfg,output),'Fe2O3:OutputExists');
assert(strcmp(before,fe2o3_sha256(output)));
reject(@() fe2o3_run(cfg,scratch),'Fe2O3:InvalidOutput');
reject(@() fe2o3_run(cfg,fullfile(scratch,'missing','file.mat')),'Fe2O3:InvalidOutput');
reject(@() fe2o3_run(cfg,'',2),'Fe2O3:InvalidOutput');
baseline=fullfile(root,'benchmark','baseline','cq-20261120-v1','sentinel.mat');
reject(@() fe2o3_run(cfg,baseline,true),'Fe2O3:ProtectedBaseline');
reject(@() fe2o3_run(cfg,fullfile(root,'benchmark','candidates','..','baseline','sentinel.mat'),true),...
    'Fe2O3:ProtectedBaseline');
if ispc, reject(@() fe2o3_run(cfg,upper(baseline),true),'Fe2O3:ProtectedBaseline'); end
if ispc
    alias=fullfile(scratch,'benchmark-alias');
    [status,message]=system(sprintf('cmd /c mklink /J "%s" "%s"',alias,fullfile(root,'benchmark')));
    assert(status==0,'Fe2O3:TestFailure','Cannot verify junction guard: %s',message);
    reject(@() fe2o3_run(cfg,fullfile(alias,'baseline','sentinel.mat'),true),'Fe2O3:ProtectedBaseline');
end
% A writer arrives after successful preflight. Neither check nor publishing may clobber it.
raced=fullfile(scratch,'raced.mat'); fe2o3_output_path(raced,'',false);
fid=fopen(raced,'w'); fprintf(fid,'another writer'); fclose(fid); hash=fe2o3_sha256(raced);
reject(@() fe2o3_save_result(result,raced,'',false),'Fe2O3:OutputExists');
assert(strcmp(hash,fe2o3_sha256(raced)));
% Verify the final exclusive filesystem primitive on this filesystem too.
try
    java.nio.file.Files.createLink(java.io.File(raced).toPath(),java.io.File(output).toPath());
    error('Fe2O3:TestFailure','Exclusive creation replaced an existing path.');
catch exception
    assert(~strcmp(exception.identifier,'Fe2O3:TestFailure'));
end
assert(strcmp(hash,fe2o3_sha256(raced)));
replacement=result; replacement.provenance.test_marker='explicit replacement';
fe2o3_save_result(replacement,output,'',true);
loaded=load(output); assert(isequaln(loaded.fe2o3_result,replacement));
fprintf('PASS existing/path/baseline guards, late writer protection and explicit replacement\n');

geometry=result.geometry;
[x,y,z,info]=fe2o3_load_restart(output,geometry);
assert(isequal(x,result.trajectory.mmx(:,:,:,end)) && info.selected_frame==4);
assert(isequal(y,result.trajectory.mmy(:,:,:,end)) && isequal(z,result.trajectory.mmz(:,:,:,end)));
restart=cfg; restart.load_fixed_atom_distrib=1; restart.restart_file=output;
reject(@() fe2o3_run(restart,output,true),'Fe2O3:RestartCollision');
restarted=fe2o3_run(restart);
assert(isequal(restarted.initial_state,struct('mx',x,'my',y,'mz',z)));
assert(strcmp(restarted.provenance.restart.sha256,fe2o3_sha256(output)));
legacy=fullfile(scratch,'legacy.mat');
mmx=result.trajectory.mmx; mmy=result.trajectory.mmy; mmz=result.trajectory.mmz;
save(legacy,'mmx','mmy','mmz'); [lx,ly,lz]=fe2o3_load_restart(legacy,geometry);
assert(isequal(x,lx) && isequal(y,ly) && isequal(z,lz));
malformed=fullfile(scratch,'malformed.mat');
for index=1:5
    mmx=result.trajectory.mmx; mmy=result.trajectory.mmy; mmz=result.trajectory.mmz;
    occupied=find(geometry.atomtype_,1); empty=find(~geometry.atomtype_,1);
    if index==1, mmx(occupied)=NaN;
    elseif index==2, mmx(occupied)=0; mmy(occupied)=0; mmz(occupied)=0;
    elseif index==3, mmx(empty)=1;
    elseif index==4, mmy=mmy(:,:,:,1);
    else, mmx=complex(mmx,ones(size(mmx))); end
    save(malformed,'mmx','mmy','mmz');
    reject(@() fe2o3_load_restart(malformed,geometry),'Fe2O3:InvalidRestart');
end
for name={'atomtype_layer1gr','atomtype_layer1or','atomtype_layer2blue',...
        'atomtype_layer2p','atomtype_layer3red','atomtype_layer3black','atomtype_'}
    fe2o3_result=result; fe2o3_result.geometry.(name{1})(1)=1-fe2o3_result.geometry.(name{1})(1);
    save(malformed,'fe2o3_result');
    reject(@() fe2o3_load_restart(malformed,geometry),'Fe2O3:InvalidRestart');
end
fe2o3_result=result; fe2o3_result.schema_version='0.0'; save(malformed,'fe2o3_result');
reject(@() fe2o3_load_restart(malformed,geometry),'Fe2O3:InvalidRestart');
fprintf('PASS legacy/new restart, last saved sample, full-history corruption and all seven masks\n');

% Execute both consumer branches with disposable invisible figures.
visibility=get(groot,'defaultFigureVisible'); set(groot,'defaultFigureVisible','off');
plot_cleanup=onCleanup(@() set(groot,'defaultFigureVisible',visibility));
plot_new(result); plot_legacy(result);
shadow=fullfile(scratch,'shadow'); mkdir(shadow);
fid=fopen(fullfile(shadow,'field_calc.m'),'w'); fprintf(fid,'error(''shadow must not run'');\n'); fclose(fid);
cd(shadow); rehash;
reject(@() fe2o3_run(cfg),'Fe2O3:ShadowedSource'); cd(oldcwd);
fprintf('PASS both plotting schemas and shadowed-source rejection\nInterface checks: %s\n',scratch);
% Hash actual executed bytes even when the parent Git commit is unchanged.
copy=fullfile(scratch,'modified-source'); mkdir(copy);
for item=dir(fullfile(root,'*.m'))'
    copyfile(fullfile(root,item.name),fullfile(copy,item.name));
end
fid=fopen(fullfile(copy,'field_calc.m'),'a'); fprintf(fid,'\n%% Test-owned source identity marker.\n'); fclose(fid);
cd(copy); clear fe2o3_run field_calc; rehash;
modified=fe2o3_run(cfg);
assert(isequal(modified.trajectory,result.trajectory));
items=modified.provenance.source.project_sources;
item=items(strcmp({items.file},'field_calc.m'));
assert(strcmp(item.sha256,fe2o3_sha256(fullfile(copy,'field_calc.m'))) && ...
    ~strcmp(item.sha256,fe2o3_sha256(fullfile(root,'field_calc.m'))));
cd(oldcwd); clear fe2o3_run field_calc; rehash;
fprintf('PASS modified-source provenance with unchanged numerical results\n');
end

function plot_new(fe2o3_result)
before=findall(groot,'Type','figure');
One_data_curve;
delete(setdiff(findall(groot,'Type','figure'),before));
end

function plot_legacy(result)
mmx=result.trajectory.mmx; mmy=result.trajectory.mmy; mmz=result.trajectory.mmz;
t=result.trajectory.t; BDSOTRE=result.derived.BDSOTRE; psjSHE=result.derived.psjSHE;
before=findall(groot,'Type','figure');
One_data_curve;
delete(setdiff(findall(groot,'Type','figure'),before));
end

function reject(action,identifier)
try, action(); catch exception
    assert(strcmp(exception.identifier,identifier),'Fe2O3:TestFailure',...
        'Expected %s, received %s: %s',identifier,exception.identifier,exception.message);
    return
end
error('Fe2O3:TestFailure','Expected rejection: %s.',identifier);
end

function restore_environment(oldpath,oldcwd,saved_rng)
cd(oldcwd); path(oldpath); rng(saved_rng);
end
