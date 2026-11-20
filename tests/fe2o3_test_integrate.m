function out = fe2o3_test_integrate(p,s,dt,steps,chunks,stride,trace_enabled)
% Exercise unchanged production integration with disposable field logging.
% Caller supplies a fresh wrapper directory as p.wrapper when tracing.
names=fieldnames(p);
for j=1:numel(names), eval([names{j} '=p.(names{j});']); end
mx_init=reshape(s(:,1),size(atomtype_));
my_init=reshape(s(:,2),size(atomtype_));
mz_init=reshape(s(:,3),size(atomtype_));
tstep=dt; gpurun_number=chunks; gpusave=steps*dt; savetstep=stride;
assert(mod(steps,stride)==0);
cfg=fe2o3_default_config;
cfg.tstep=dt; cfg.gpusave=gpusave; cfg.gpurun_number=chunks; cfg.savetstep=stride;
[~,clock]=fe2o3_validate_config(cfg);
gpusteps=clock.gpusteps; runtime=clock.runtime;
final_m_savestep=clock.final_m_savestep; t=clock.t;
review_field_source=fileread(fullfile(p.root,'field_calc.m'));
review_trace_enabled=trace_enabled; review_trace={};
eval(fileread(fullfile(p.root,'integrate_llg.m')));
out.x=mmx; out.y=mmy; out.z=mmz;
out.dx=hdmi_allx; out.dy=hdmi_ally; out.dz=hdmi_allz;
out.t=t;
out.trace=review_trace;
end
