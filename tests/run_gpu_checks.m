function run_gpu_checks
% Port of PR #20 independent stage/analytic checks plus sample/cache checks.
root=fileparts(fileparts(mfilename('fullpath')));
oldpath=path; oldcwd=pwd; saved_rng=rng;
cleanup=onCleanup(@() restore_environment(oldpath,oldcwd,saved_rng));
addpath(root); addpath(fullfile(root,'tests'));
device=gpuDevice;
fprintf('MATLAB %s; %s; GPU %s; compute %s; driver %s\n',...
    version,computer,device.Name,device.ComputeCapability,string(device.GraphicsDriverVersion));
scratch_root=fullfile(root,'.cq-checks');
if ~isfolder(scratch_root), mkdir(scratch_root); end
scratch=tempname(scratch_root); mkdir(scratch);
wrapper=["eval(review_field_source);";...
    "if review_trace_enabled";...
    "review_entry.s=gather([mmxtmp(:),mmytmp(:),mmztmp(:)]);";...
    "review_entry.b=gather([hhx(:),hhy(:),hhz(:)]);";...
    "review_entry.exchange=gather([hex_x(:),hex_y(:),hex_z(:)]);";...
    "review_entry.anisotropy=gather([hani_x(:),hani_y(:),hani_z(:)]);";...
    "review_entry.dmi=gather([hdmi_x(:),hdmi_y(:),hdmi_z(:)]);";...
    "review_entry.dipole=zeros(size(review_entry.s));";...
    "if dipolemode, review_entry.dipole=gather([hdipo_x(:),hdipo_y(:),hdipo_z(:)]); end";...
    "review_trace{end+1}=review_entry;";"end"];
fid=fopen(fullfile(scratch,'field_calc.m'),'w');
assert(fid~=-1); fprintf(fid,'%s\n',wrapper); fclose(fid);
cd(scratch); diary('checks.log');
assert(strcmp(which('field_calc'),fullfile(scratch,'field_calc.m')));
fprintf('Disposable instrumentation/log: %s\n',scratch);
p=fe2o3_test_fixture([5,9,6]); occupied=p.atomtype_(:)~=0;
rng(20,'twister'); s=randn(numel(occupied),3); s=normalize(s,occupied);
p.BDSOTRE=0.18; p.BDSTTRE=0.04; h=1e-16;
r=fe2o3_test_integrate(p,s,h,2,1,1,true);
assert(numel(r.trace)==8);
k=cell(1,4); expected=s; stage_error=0; field_error=0;
for stage=1:4
    stage_error=max(stage_error,max(abs(r.trace{stage}.s-expected),[],'all'));
    fresh=fe2o3_test_field(p,expected);
    field_error=max(field_error,max(abs(r.trace{stage}.b-fresh.b),[],'all'));
    k{stage}=vector_rhs(expected,fresh.b,p);
    if stage<3, expected=s+h/2*k{stage};
    elseif stage==3, expected=s+h*k{stage}; end
end
expected=normalize(s+h/6*(k{1}+2*k{2}+2*k{3}+k{4}),occupied);
step_error=max(abs(frame(r,2)-expected),[],'all');
assert(stage_error<2e-12 && field_error<2e-10 && step_error<2e-12);
for term={'exchange','anisotropy','dmi'}
    assert(max(abs(r.trace{2}.(term{1})-r.trace{1}.(term{1})),[],'all')>1e-12);
end
check_invariants(r,occupied);
fprintf('PASS RK4 stages: state %.12g; field %.12g T; step %.12g\n',stage_error,field_error,step_error);
% Exact nonlinear Bz=b*sz solution, synthetic parameters, as in PR #20.
q=p; q.gam=1; q.mus=1; q.alp=0.3;
q.J1=0; q.J2=0; q.J3=0; q.J4=0; q.Ksim1=0.75; q.hk6=0;
q.Dsim1=0; q.Dsim3=0; q.Dsim4=0; q.Hext=[0,0,0];
q.BDSOTRE=0; q.BDSTTRE=0; q.chi=0;
z0=0.4; initial=repmat([sqrt(1-z0^2),0,z0],numel(occupied),1);
initial(~occupied,:)=0; b=2*q.Ksim1/q.mus; scale=q.gam/(1+q.alp^2);
z=1/sqrt(1+(1/z0^2-1)*exp(-2*scale*q.alp*b));
phi=(atanh(z)-atanh(z0))/q.alp;
exact=[sqrt(1-z^2)*cos(phi),sqrt(1-z^2)*sin(phi),z];
counts=[8,16,32,64]; errors=zeros(size(counts));
for index=1:numel(counts)
    r=fe2o3_test_integrate(q,initial,1/counts(index),counts(index)+1,1,1,false);
    check_invariants(r,occupied);
    actual=frame(r,counts(index)+1);
    errors(index)=max(abs(actual(occupied,:)-exact),[],'all');
end
orders=log2(errors(1:end-1)./errors(2:end));
assert(all(diff(errors)<0) && all(orders>3.8 & orders<4.2));
fprintf('PASS nonlinear analytic errors: %.12g %.12g %.12g %.12g; orders %.8f %.8f %.8f\n',errors,orders);
% Spin/DMI pairing and chunk identity for the dipole-disabled path.
dense=fe2o3_test_integrate(p,s,h,8,3,1,false);
for stride=[1,2,8]
    r=fe2o3_test_integrate(p,s,h,8,3,stride,false);
    single=fe2o3_test_integrate(p,s,h,24,1,stride,false);
    assert(isequal(r.x,single.x) && isequal(r.y,single.y) && isequal(r.z,single.z));
    assert(isequal(r.t,(0:24/stride-1)*stride*h));
    for component={'x','y','z','dx','dy','dz'}
        assert(isequal(r.(component{1}),dense.(component{1})(:,:,:,1:stride:end)));
    end
    assert(isequal(r.t,dense.t(1:stride:end)));
    for index=1:numel(r.t)
        f=fe2o3_test_field(p,frame(r,index));
        actual=[reshape(r.dx(:,:,:,index),[],1),reshape(r.dy(:,:,:,index),[],1),reshape(r.dz(:,:,:,index),[],1)];
        assert(max(abs(actual-f.dmi),[],'all')<2e-10);
    end
    check_invariants(r,occupied);
end
fprintf('PASS strides 1/2/8: chunk equivalence, every sampled DMI field, raw invariants\n');
% Keep existing cache behavior: refresh at global updates 0,50,52,102.
p.dipolemode=1;
r=fe2o3_test_integrate(p,s,h,52,2,1,true);
check_invariants(r,occupied);
assert(numel(r.trace)==416);
for update=0:103
    first=4*update+1;
    if ismember(update,[0,50,52,102])
        f=fe2o3_test_field(p,r.trace{first}.s); cached=f.dipole;
    end
    for stage=first:first+3
        assert(max(abs(r.trace{stage}.dipole-cached),[],'all')<2e-12);
    end
end
f=fe2o3_test_field(p,r.trace{2}.s);
gap=max(abs(f.dipole-r.trace{2}.dipole),[],'all'); assert(gap>1e-12);
fprintf('PASS retained cache refreshes 0/50/52/102; stage-2 fresh gap %.12g T\n',gap);
p.dipolemode=0; p.BDSOTRE=0; p.BDSTTRE=0;
for solver=[0,2]
    p.rk4=solver;
    r=fe2o3_test_integrate(p,s,h,8,3,1,false);
    single=fe2o3_test_integrate(p,s,h,24,1,1,false);
    check_invariants(r,occupied);
    assert(isequal(r.x,single.x) && isequal(r.y,single.y) && isequal(r.z,single.z));
    fprintf('PASS legacy mode %d: finite/norm/empty and chunk-state continuity only\n',solver);
end
fprintf('PASS GPU checks; no production output or accepted reference written.\n');
end

function k=vector_rhs(s,b,p)
v=-cross(s,b,2)-p.alp*cross(s,cross(s,b,2),2);
polarizations={[p.psjSHEx,p.psjSHEy,p.psjSHEz],[p.psjSTTx,p.psjSTTy,p.psjSTTz]};
amplitudes=[p.BDSOTRE,p.BDSTTRE];
for channel=1:2
    sigma=repmat(polarizations{channel},size(s,1),1);
    one=cross(s,sigma,2); two=cross(s,one,2); bd=amplitudes(channel); bf=p.chi*bd;
    v=v+(p.alp*bd+bf)*one+(-bd+p.alp*bf)*two;
end
k=p.gam/(1+p.alp^2)*v;
end

function s=frame(r,index)
s=[reshape(r.x(:,:,:,index),[],1),reshape(r.y(:,:,:,index),[],1),reshape(r.z(:,:,:,index),[],1)];
end

function s=normalize(s,occupied)
s(occupied,:)=s(occupied,:)./vecnorm(s(occupied,:),2,2); s(~occupied,:)=0;
end

function check_invariants(r,occupied)
for index=1:size(r.x,4)
    s=frame(r,index); assert(all(isfinite(s),'all'));
    assert(all(s(~occupied,:)==0,'all'));
    assert(max(abs(vecnorm(s(occupied,:),2,2)-1))<2e-12);
end
end

function restore_environment(oldpath,oldcwd,saved_rng)
diary off; cd(oldcwd); path(oldpath); rng(saved_rng);
end
