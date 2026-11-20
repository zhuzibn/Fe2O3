function output_root = capture_legacy_reference(output_root)
% Prepare two isolated original-code suites for review. Never accept data.
here=fileparts(mfilename('fullpath')); root=fileparts(here);
oldpath=path; oldcwd=pwd; saved_rng=rng;
cleanup=onCleanup(@() restore_environment(oldpath,oldcwd,saved_rng));
addpath(here); addpath(root); addpath(fullfile(root,'tests'));
source=jsondecode(fileread(fullfile(here,'legacy_source_manifest.json')));
for item=source.files(:)'
    assert(strcmp(fe2o3_file_hash(fullfile(root,item.file)),item.sha256),...
        'Fe2O3:LegacySourceChanged','Original source changed: %s. Do not capture replacement baselines.',item.file);
end
device=gpuDevice;
environment=struct('matlab_version',version,'platform',computer,...
    'parallel_toolbox',ver('parallel'),'gpu_name',device.Name,...
    'compute_capability',device.ComputeCapability,'driver',device.GraphicsDriverVersion,...
    'toolkit',device.ToolkitVersion,'total_memory',device.TotalMemory);
if nargin<1, output_root=fe2o3_candidate_directory;
else, output_root=fe2o3_candidate_directory(output_root); end
output_root=fullfile(output_root,'reference-proposal'); mkdir(output_root);
write_json(fullfile(output_root,'status.json'),struct('status','incomplete-reference-proposal','source_commit',source.commit));
cases=fe2o3_regression_cases;
for repeat=1:2
    suite=fullfile(output_root,sprintf('repeat-%d',repeat)); mkdir(suite);
    original=fullfile(suite,'original-source'); mkdir(original);
    for item=source.files(:)', copyfile(fullfile(root,item.file),fullfile(original,item.file)); end
    empty_artifacts=struct('file',{},'sha256',{});
    empty_cases=struct('id',{},'record_file',{},'record_sha256',{},'artifacts',{});
    manifest=struct('schema_version','1.0','status','reference-proposal',...
        'source_commit',source.commit,'environment',environment,...
        'created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')),...
        'cases',empty_cases);
    for index=1:numel(cases)
        spec=cases(index); case_dir=fullfile(suite,spec.id); mkdir(case_dir);
        for item=source.files(:)', copyfile(fullfile(original,item.file),fullfile(case_dir,item.file)); end
        % Fixed nonuniform input is materialized before launching old main.
        p=fe2o3_test_fixture(spec.core);
        if spec.restart
            [rr,cc,kk]=ndgrid(1:spec.core(1),1:spec.core(2),1:spec.core(3));
            v=cat(4,(-1).^(kk+1),0.07*sin(rr+2*cc+kk),0.05*cos(2*rr-cc+3*kk));
            v=v./sqrt(sum(v.^2,4));
            mmx=zeros(size(p.atomtype_)); mmy=mmx; mmz=mmx;
            mmx(3:end-2,3:end-2,:)=v(:,:,:,1);
            mmy(3:end-2,3:end-2,:)=v(:,:,:,2);
            mmz(3:end-2,3:end-2,:)=v(:,:,:,3);
            mmx=mmx.*p.atomtype_; mmy=mmy.*p.atomtype_; mmz=mmz.*p.atomtype_;
            save(fullfile(case_dir,'Relax_0.mat'),'mmx','mmy','mmz');
        end
        text=fileread(fullfile(original,'main.m'));
        [patched,changes]=prepare_main(text,spec);
        write_text(fullfile(case_dir,'main.m'),patched);
        write_json(fullfile(case_dir,'substitutions.json'),changes);
        % cwd precedes path. Original script has clear all; execute in its own
        % function scope so it cannot clear orchestrator state or onCleanup.
        cd(case_dir); rehash;
        assert(strcmp(which('main'),fullfile(case_dir,'main.m')));
        for item=source.files(:)'
            if ~strcmp(item.file,'main.m')
                assert(strcmp(fe2o3_file_hash(fullfile(case_dir,item.file)),item.sha256));
            end
        end
        started=tic; execute_original_main; elapsed=toc(started);
        cd(root);
        data=load(fullfile(case_dir,'legacy_workspace.mat'),...
            'cq_configuration','cq_geometry','cq_initial_state','cq_rng_before','cq_rng_initial','cq_rng_after',...
            'mmx','mmy','mmz','hdmi_allx','hdmi_ally','hdmi_allz','t');
        record=struct('schema_version','1.0','specification',spec,...
            'configuration',data.cq_configuration,'geometry',data.cq_geometry,...
            'initial_state',data.cq_initial_state,...
            'units',struct('spin','dimensionless','field','T','time','s'),...
            'sampling','pre-update; final endpoint excluded');
        for name={'mmx','mmy','mmz','hdmi_allx','hdmi_ally','hdmi_allz','t'}
            record.trajectory.(name{1})=data.(name{1});
        end
        record.provenance=struct('source_commit',source.commit,'original_sources',source.files,...
            'environment',environment,'rng_before',data.cq_rng_before,'rng_initial',data.cq_rng_initial,...
            'rng_after',data.cq_rng_after,'elapsed_seconds',elapsed,...
            'entry_point','isolated original main.m; logged control and metadata substitutions',...
            'case_directory',case_dir,'created_utc',manifest.created_utc,...
            'git_head',git_output(root,'rev-parse HEAD'),'git_status',git_output(root,'status --short'));
        input=record.initial_state; save(fullfile(case_dir,'initial_input.mat'),'input');
        save(fullfile(case_dir,'record.mat'),'record');
        entry=struct('id',spec.id,'record_file',[spec.id '/record.mat'],...
            'record_sha256',fe2o3_file_hash(fullfile(case_dir,'record.mat')),'artifacts',empty_artifacts);
        artifacts={'legacy_workspace.mat','initial_input.mat','substitutions.json'};
        if spec.restart, artifacts{end+1}='Relax_0.mat'; end
        for item=source.files(:)'
            artifacts{end+1}=item.file;
        end
        for j=1:numel(artifacts)
            entry.artifacts(end+1)=struct('file',[spec.id '/' artifacts{j}],...
                'sha256',fe2o3_file_hash(fullfile(case_dir,artifacts{j})));
        end
        for item=source.files(:)'
            entry.artifacts(end+1)=struct('file',['original-source/' item.file],'sha256',item.sha256);
        end
        manifest.cases(index)=entry;
        write_json(fullfile(suite,'manifest.json'),manifest);
        fprintf('Captured repeat %d %s: %d samples, %.3f s\n',repeat,spec.id,numel(record.trajectory.t),elapsed);
    end
end
metrics=compare_regression(fullfile(output_root,'repeat-2'),fullfile(output_root,'repeat-1'));
proposal=struct('status','reference-proposal-awaiting-review','source_commit',source.commit,...
    'repeat_1_manifest_sha256',fe2o3_file_hash(fullfile(output_root,'repeat-1','manifest.json')),...
    'repeat_2_manifest_sha256',fe2o3_file_hash(fullfile(output_root,'repeat-2','manifest.json')),...
    'metrics',metrics,'limits',struct('spin',2e-12,'DMI_T',2e-10,'norm',2e-12));
write_json(fullfile(output_root,'status.json'),proposal);
fprintf('REFERENCE PROPOSAL ONLY: %s\n',output_root);
end

function execute_original_main
main;
end

function [text,changes]=prepare_main(text,s)
changes=struct('before',{},'after',{});
[text,changes]=replace_once(text,changes,'clear all;clc;close all;tic',...
    sprintf('clear all;clc;close all;tic\ncq_rng_before=rng;'));
pairs={'natomW=20;natomL=20;natomH=6;',sprintf('natomW=%d;natomL=%d;natomH=%d;',s.core);...
    'rk4=1;',sprintf('rk4=%d;',s.solver);...
    'load_fixed_atom_distrib=0;',sprintf('load_fixed_atom_distrib=%d;',s.restart);...
    'dipolemode=0;',sprintf('dipolemode=%d;',s.dipole);...
    'jcSOT=0e12;',sprintf('jcSOT=%.17g;',s.jcSOT);...
    'jcSTT=0e9;',sprintf('jcSTT=%.17g;',s.jcSTT);...
    'STT_DLT=0;',sprintf('STT_DLT=%d;',s.STT_DLT);...
    'STT_FLT=0;',sprintf('STT_FLT=%d;',s.STT_FLT);...
    'gpusave=1e-12;',sprintf('gpusave=%.17g;',s.steps*s.tstep);...
    'gpurun_number=25;',sprintf('gpurun_number=%d;',s.chunks);...
    'savetstep=5000;',sprintf('savetstep=%d;',s.stride)};
for j=1:size(pairs,1), [text,changes]=replace_once(text,changes,pairs{j,1},pairs{j,2}); end
% Metadata snapshots run before main's legacy cleanup removes input aliases.
names={'natomW','natomL','natomH','rk4','DMIenable','dipolemode',...
    'load_fixed_atom_distrib','save_fixed_atom_distrib','Ksim1','hk6','J1','J2','J3','J4',...
    'g','gam','mus','a','c','tz','ms','Dsim1','Dsim3','Dsim4','alp','jcSOT','jcSTT','Hext',...
    'SOT_DLT','SOT_FLT','psj_theta','psj_phi','psjSHE','thetaSH','chi','BDSOTRE','BDSOTTM',...
    'STT_DLT','STT_FLT','psjSTT','etaSTT','chiSTT','BDSTTRE','BDSTTTM',...
    'gpusave','gpurun_number','tstep','savetstep','gpusteps','runtime','final_m_savestep',...
    'ele','mub','hbar','mu_0','kb'};
geometry={'atomtype_','atomtype_layer1gr','atomtype_layer1or','atomtype_layer2blue',...
    'atomtype_layer2p','atomtype_layer3red','atomtype_layer3black'};
hook=sprintf('cq_configuration=struct;\ncq_geometry=struct;\n');
for j=1:numel(names), hook=[hook sprintf('cq_configuration.%s=%s;\n',names{j},names{j})]; end
for j=1:numel(geometry), hook=[hook sprintf('cq_geometry.%s=%s;\n',geometry{j},geometry{j})]; end
hook=[hook sprintf(['cq_geometry.array_order=''W,L,H,time'';\n'...
    'cq_initial_state=struct(''mx'',mx_init,''my'',my_init,''mz'',mz_init);\n'...
    'cq_rng_initial=rng;\n%% dynamic calc'])];
[text,changes]=replace_once(text,changes,'%% dynamic calc',hook);
[text,changes]=replace_once(text,changes,'save(''Relax.mat'')',...
    sprintf('cq_rng_after=rng;\nsave(''legacy_workspace.mat'')'));
end

function [text,changes]=replace_once(text,changes,before,after)
assert(count(text,before)==1,'Fe2O3:LegacySourceChanged','Expected exactly one substitution: %s',before);
changes(end+1)=struct('before',before,'after',after); text=strrep(text,before,after);
end

function out=git_output(root,args)
[status,value]=system(sprintf('git -C "%s" %s',root,args));
out=struct('exit_status',status,'value',strtrim(value));
end

function write_json(file,value)
write_text(file,[jsonencode(value,'PrettyPrint',true) newline]);
end

function write_text(file,value)
fid=fopen(file,'w'); assert(fid~=-1); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',value);
end

function restore_environment(oldpath,oldcwd,saved_rng)
cd(oldcwd); path(oldpath); rng(saved_rng);
end
