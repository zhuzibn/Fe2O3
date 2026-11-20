function metrics = compare_regression(candidate_dir,reference_dir)
% Compare all cases/samples and recomputed invariants; failure raises.
% Explicit reference paths permit proposal repeat checks, never promotion.
if nargin<2
    path=fullfile(fileparts(mfilename('fullpath')),'baseline_manifest.json');
    assert(isfile(path),'Fe2O3:MissingBaseline','No accepted baseline manifest; do not bootstrap one.');
    accepted=jsondecode(fileread(path));
    assert(strcmp(accepted.status,'accepted'),'Fe2O3:MissingBaseline','Reference is not accepted.');
    reference_dir=fullfile(fileparts(path),accepted.directory);
    assert(strcmp(fe2o3_file_hash(fullfile(reference_dir,'manifest.json')),accepted.manifest_sha256),...
        'Fe2O3:ArtifactHash','Accepted manifest hash mismatch.');
end
expected=fe2o3_regression_cases;
candidate=read_suite(candidate_dir,expected); reference=read_suite(reference_dir,expected);
metrics=struct('case_id',{},'spin_error',{},'dmi_error_T',{},...
    'candidate_norm_error',{},'reference_norm_error',{});
for index=1:numel(expected)
    a=read_record(candidate_dir,candidate.cases(index));
    b=read_record(reference_dir,reference.cases(index));
    na=validate_record(a,expected(index)); nb=validate_record(b,expected(index));
    for name={'specification','configuration','geometry','initial_state','units','sampling'}
        assert(isequaln(a.(name{1}),b.(name{1})),...
            'Fe2O3:IncompatibleInput','%s: incompatible %s.',expected(index).id,name{1});
    end
    assert(isequal(a.trajectory.t,b.trajectory.t),...
        'Fe2O3:WrongClock','%s: timestamps differ.',expected(index).id);
    spin_error=0; dmi_error=0;
    for name={'mmx','mmy','mmz'}
        spin_error=max(spin_error,max(abs(a.trajectory.(name{1})-b.trajectory.(name{1})),[],'all'));
    end
    for name={'hdmi_allx','hdmi_ally','hdmi_allz'}
        dmi_error=max(dmi_error,max(abs(a.trajectory.(name{1})-b.trajectory.(name{1})),[],'all'));
    end
    assert(spin_error<=2e-12 && dmi_error<=2e-10,...
        'Fe2O3:RegressionMismatch','%s: spin %.12g, DMI %.12g T.',expected(index).id,spin_error,dmi_error);
    metrics(index)=struct('case_id',expected(index).id,'spin_error',spin_error,...
        'dmi_error_T',dmi_error,'candidate_norm_error',na,'reference_norm_error',nb);
    fprintf('PASS %s: spin %.12g; DMI %.12g T; norms %.12g / %.12g\n',...
        expected(index).id,spin_error,dmi_error,na,nb);
end
end

function m=read_suite(folder,expected)
filename=fullfile(folder,'manifest.json');
assert(isfile(filename),'Fe2O3:MissingSuite','Missing suite manifest: %s',filename);
m=jsondecode(fileread(filename));
assert(isstruct(m) && isscalar(m) && isfield(m,'schema_version') && strcmp(m.schema_version,'1.0'),...
    'Fe2O3:InvalidSuite','Invalid suite schema.');
assert(isfield(m,'cases') && numel(m.cases)==numel(expected) && ...
    isequal({m.cases.id},{expected.id}),'Fe2O3:MissingCase','Suite inventory/order does not match all six required cases.');
end

function r=read_record(folder,entry)
filename=fullfile(folder,entry.record_file);
assert(strcmp(fe2o3_file_hash(filename),entry.record_sha256),...
    'Fe2O3:ArtifactHash','Record hash mismatch: %s',filename);
for artifact=entry.artifacts(:)'
    assert(strcmp(fe2o3_file_hash(fullfile(folder,artifact.file)),artifact.sha256),...
        'Fe2O3:ArtifactHash','Supporting artifact hash mismatch: %s',artifact.file);
end
data=load(filename,'record');
assert(isfield(data,'record'),'Fe2O3:InvalidRecord','Missing record.'); r=data.record;
end

function norm_error=validate_record(r,spec)
required={'schema_version','specification','configuration','geometry','initial_state','units','sampling','trajectory','provenance'};
assert(isstruct(r) && isscalar(r) && all(isfield(r,required)) && strcmp(r.schema_version,'1.0'),...
    'Fe2O3:InvalidRecord','Invalid record schema.');
assert(isequal(r.specification,spec),'Fe2O3:IncompatibleInput','Unexpected case specification.');
assert(isequal(r.units,struct('spin','dimensionless','field','T','time','s')) && ...
    strcmp(r.sampling,'pre-update; final endpoint excluded'),...
    'Fe2O3:InvalidRecord','Units/sampling contract changed.');
assert(isstruct(r.configuration) && isscalar(r.configuration));
finite_numeric(r.configuration);
c=r.configuration; shape=spec.core+[4,4,0]; count=spec.steps*spec.chunks/spec.stride;
assert(isequal([c.natomW,c.natomL,c.natomH],shape) && c.rk4==spec.solver && ...
    c.gpusteps==spec.steps && c.gpurun_number==spec.chunks && c.savetstep==spec.stride && ...
    c.tstep==spec.tstep && c.dipolemode==spec.dipole && c.jcSOT==spec.jcSOT && c.jcSTT==spec.jcSTT && ...
    c.STT_DLT==spec.STT_DLT && c.STT_FLT==spec.STT_FLT && c.load_fixed_atom_distrib==spec.restart,...
    'Fe2O3:IncompatibleInput','Effective configuration does not match declared case.');
g=r.geometry;
assert(isfield(g,'atomtype_') && isequal(size(g.atomtype_),shape) && ...
    all(ismember(g.atomtype_(:),[0,1])),'Fe2O3:InvalidMapping','Malformed occupancy.');
assert(isfield(g,'array_order') && strcmp(g.array_order,'W,L,H,time'));
names={'atomtype_layer1gr','atomtype_layer1or','atomtype_layer2blue','atomtype_layer2p','atomtype_layer3red','atomtype_layer3black'};
for index=1:6
    assert(isfield(g,names{index}) && isequal(size(g.(names{index})),shape(1:2)) && ...
        all(ismember(g.(names{index})(:),[0,1])),'Fe2O3:InvalidMapping','Malformed site mask.');
end
mask=logical(g.atomtype_);
assert(any(mask,'all'),'Fe2O3:InvalidMapping','Empty occupied state.');
initial={'mx','my','mz'}; squares=zeros(shape);
for index=1:3
    value=r.initial_state.(initial{index});
    assert(isnumeric(value) && isreal(value) && isequal(size(value),shape) && all(isfinite(value),'all'));
    assert(all(value(~mask)==0),'Fe2O3:InvalidSpin','Nonzero initial empty sites.');
    squares=squares+value.^2;
end
assert(max(abs(sqrt(squares(mask))-1))<=2e-12,'Fe2O3:InvalidSpin','Invalid initial norms.');
names={'mmx','mmy','mmz','hdmi_allx','hdmi_ally','hdmi_allz'};
for index=1:numel(names)
    value=r.trajectory.(names{index});
    assert(isnumeric(value) && isreal(value) && isequal(size(value),[shape,count]) && all(isfinite(value),'all'),...
        'Fe2O3:InvalidArray','Malformed/nonfinite %s.',names{index});
    reshaped=reshape(value,[],count);
    assert(all(reshaped(~mask(:),:)==0,'all'),'Fe2O3:InvalidSpin','Nonzero empty trajectory sites.');
end
clock=(0:count-1)*spec.stride*spec.tstep;
assert(isequal(r.trajectory.t,clock),'Fe2O3:WrongClock','Clock does not match physical sample contract.');
squares=r.trajectory.mmx.^2+r.trajectory.mmy.^2+r.trajectory.mmz.^2;
norms=reshape(sqrt(squares),[],count);
norm_error=max(abs(norms(mask(:),:)-1),[],'all');
assert(norm_error<=2e-12,'Fe2O3:InvalidSpin','Invalid raw trajectory norms.');
for index=1:3
    assert(isequal(r.trajectory.(names{index})(:,:,:,1),r.initial_state.(initial{index})),...
        'Fe2O3:IncompatibleInput','Initial frame differs from supplied initial state.');
end
end

function finite_numeric(s)
for name=fieldnames(s)'
    v=s.(name{1});
    if isnumeric(v), assert(isreal(v) && all(isfinite(v),'all'),'Fe2O3:InvalidConfig','Nonfinite/complex configuration.');
    elseif isstruct(v), finite_numeric(v); end
end
end
