function cases = fe2o3_regression_cases
% Approved small-suite inventory; overrides only, not production defaults.
ids={'small-relax','driven-restart','default-core-short','heun-legacy','pc4-legacy','dipole-cache'};
cores=[5,9,6;5,9,6;20,20,6;5,9,6;5,9,6;5,9,6];
steps=[8,8,8,8,8,52]; chunks=[3,3,2,3,3,2]; strides=[2,1,2,1,1,1];
solvers=[1,1,1,0,2,1]; restart=[0,1,0,1,1,1];
cases=struct([]);
for index=1:numel(ids)
    cases(index).id=ids{index}; cases(index).core=cores(index,:);
    cases(index).steps=steps(index); cases(index).chunks=chunks(index);
    cases(index).stride=strides(index); cases(index).solver=solvers(index);
    cases(index).restart=restart(index); cases(index).dipole=double(index==6);
    cases(index).jcSOT=6e12*double(index==2);
    cases(index).jcSTT=1e9*double(index==2);
    cases(index).STT_DLT=double(index==2); cases(index).STT_FLT=double(index==2);
    cases(index).tstep=1e-16;
end
end
