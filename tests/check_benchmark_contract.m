function check_benchmark_contract(proposal_root)
% Deliberate corruptions in disposable copies, never actual references.
root=fileparts(fileparts(mfilename('fullpath')));
oldpath=path; oldcwd=pwd; cleanup=onCleanup(@() restore_environment(oldpath,oldcwd));
addpath(fullfile(root,'benchmark')); cd(root);
if nargin==0
    accepted=jsondecode(fileread(fullfile(root,'benchmark','baseline_manifest.json')));
    reference=fullfile(root,'benchmark',accepted.directory); candidate=reference;
    compare_regression(reference); % Enforce the tracked acceptance pin first.
else
    reference=fullfile(proposal_root,'repeat-1'); candidate=fullfile(proposal_root,'repeat-2');
end
before={fe2o3_file_hash(fullfile(reference,'manifest.json')),fe2o3_file_hash(fullfile(candidate,'manifest.json'))};
compare_regression(candidate,reference);
work=fe2o3_candidate_directory; mutated=fullfile(work,'mutated');
copyfile(candidate,mutated);
fprintf('Mutation fixtures: %s\n',work);
manifest=jsondecode(fileread(fullfile(mutated,'manifest.json')));
entry=manifest.cases(1); data=load(fullfile(mutated,entry.record_file),'record'); original=data.record;
sites=find(original.geometry.atomtype_); at=sites(1); shape=size(original.geometry.atomtype_);
[r,c,k]=ind2sub(shape,at); rejected=0;
mutations={'intermediate-spin','intermediate-DMI','wrong-clock','wrong-config','wrong-mapping',...
    'wrong-initial','nan','inf','forged-norm-summary','wrong-shape','wrong-schema'};
identifiers={'RegressionMismatch','RegressionMismatch','WrongClock','IncompatibleInput','IncompatibleInput',...
    'IncompatibleInput','InvalidArray','InvalidArray','InvalidSpin','InvalidArray','InvalidRecord'};
for mutation=1:numel(mutations)
    name=mutations(mutation);
    record=original;
    switch name{1}
        case 'intermediate-spin'
            x=record.trajectory.mmx(r,c,k,2); y=record.trajectory.mmy(r,c,k,2);
            record.trajectory.mmx(r,c,k,2)=cos(0.1)*x-sin(0.1)*y;
            record.trajectory.mmy(r,c,k,2)=sin(0.1)*x+cos(0.1)*y;
        case 'intermediate-DMI', record.trajectory.hdmi_allx(r,c,k,2)=record.trajectory.hdmi_allx(r,c,k,2)+1e-6;
        case 'wrong-clock', record.trajectory.t(2)=record.trajectory.t(2)+1e-20;
        case 'wrong-config', record.configuration.Ksim1=record.configuration.Ksim1*2;
        case 'wrong-mapping', record.geometry.atomtype_layer1gr(r,c)=1-record.geometry.atomtype_layer1gr(r,c);
        case 'wrong-initial'
            record.initial_state.mx(r,c,k)=-record.initial_state.mx(r,c,k);
            record.trajectory.mmx(r,c,k,1)=record.initial_state.mx(r,c,k);
        case 'nan', record.trajectory.mmx(r,c,k,2)=NaN;
        case 'inf', record.trajectory.hdmi_allx(r,c,k,2)=Inf;
        case 'forged-norm-summary'
            record.trajectory.mmx(r,c,k,2)=2; record.norm_error_summary=0;
        case 'wrong-shape', record.trajectory.mmx=record.trajectory.mmx(:,:,:,1:end-1);
        case 'wrong-schema', record.schema_version='999';
    end
    save(fullfile(mutated,entry.record_file),'record');
    current=manifest; current.cases(1).record_sha256=fe2o3_file_hash(fullfile(mutated,entry.record_file));
    write_manifest(mutated,current);
    must_fail(@() compare_regression(mutated,reference),name{1},identifiers{mutation}); rejected=rejected+1;
end
% Restore only owned fixture files; original suite remains read-only.
copyfile(fullfile(candidate,entry.record_file),fullfile(mutated,entry.record_file));
current=manifest; current.cases=current.cases([1,3:end]); write_manifest(mutated,current);
must_fail(@() compare_regression(mutated,reference),'missing-case','MissingCase'); rejected=rejected+1;
write_manifest(mutated,manifest);
% Invalid baseline is rejected too; its manifest still records original hash.
record=original; record.trajectory.mmx(r,c,k,2)=NaN;
save(fullfile(mutated,entry.record_file),'record');
must_fail(@() compare_regression(candidate,mutated),'baseline-hash-mismatch','ArtifactHash'); rejected=rejected+1;
current=manifest; current.cases(1).record_sha256=fe2o3_file_hash(fullfile(mutated,entry.record_file));
write_manifest(mutated,current);
must_fail(@() compare_regression(candidate,mutated),'baseline-nonfinite','InvalidArray'); rejected=rejected+1;
% Exercise missing acceptance without moving/deleting the actual manifest.
missing=fullfile(work,'missing-acceptance'); mkdir(missing);
text=fileread(fullfile(root,'benchmark','compare_regression.m'));
text=strrep(text,'function metrics = compare_regression(', 'function metrics = missing_baseline_compare(');
fid=fopen(fullfile(missing,'missing_baseline_compare.m'),'w'); fprintf(fid,'%s',text); fclose(fid);
addpath(missing);
must_fail(@() missing_baseline_compare(candidate),'missing-accepted-baseline','MissingBaseline'); rejected=rejected+1;
must_fail(@() compare_regression(candidate,fullfile(work,'absent-raw-suite')), ...
    'missing-raw-suite','MissingSuite'); rejected=rejected+1;
sentinel=fullfile(work,'sentinel.txt'); fid=fopen(sentinel,'w'); fprintf(fid,'preserve'); fclose(fid);
sentinel_hash=fe2o3_file_hash(sentinel);
paths={work,sentinel,fullfile(root,'benchmark','baseline'),...
    fullfile(root,'benchmark','candidates','..','baseline','nested'),...
    fullfile('benchmark','baseline','nested')};
if ispc, paths{end+1}=upper(fullfile(root,'benchmark','baseline','nested')); end
for index=1:numel(paths)
    identifier='ProtectedBaseline'; if index<=2, identifier='OutputExists'; end
    must_fail(@() fe2o3_candidate_directory(paths{index}),sprintf('protected-path-%d',index),identifier);
    rejected=rejected+1;
end
% Alias an existing parent; never create/write accepted baseline storage.
if ispc
    alias=fullfile(work,'benchmark-alias');
    [status,message]=system(sprintf('cmd /c mklink /J "%s" "%s"',alias,fullfile(root,'benchmark')));
    if status==0
        must_fail(@() fe2o3_candidate_directory(fullfile(alias,'baseline','nested')), 'junction-path','ProtectedBaseline');
        rejected=rejected+1;
    else
        fprintf('UNVERIFIED junction creation (status %d): %s\n',status,message);
    end
end
assert(strcmp(fe2o3_file_hash(sentinel),sentinel_hash));
assert(strcmp(fe2o3_file_hash(fullfile(reference,'manifest.json')),before{1}));
assert(strcmp(fe2o3_file_hash(fullfile(candidate,'manifest.json')),before{2}));
% Rechecks every supporting artifact hash as well as the data.
compare_regression(candidate,reference);
fprintf('PASS %d rejection checks; original suite artifacts unchanged.\n',rejected);
end

function must_fail(action,name,identifier)
failed=false;
try, action(); catch ME
    assert(strcmp(ME.identifier,['Fe2O3:' identifier]),'Fe2O3:WrongFailure',...
        '%s failed for unexpected reason %s: %s',name,ME.identifier,ME.message);
    failed=true; fprintf('REJECTED %s: %s\n',name,ME.identifier);
end
assert(failed,'Fe2O3:FalsePass','Corruption/path unexpectedly passed: %s',name);
end

function restore_environment(oldpath,oldcwd)
cd(oldcwd); path(oldpath);
end

function write_manifest(folder,value)
fid=fopen(fullfile(folder,'manifest.json'),'w'); assert(fid~=-1);
cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
