function [x,y,z,info] = fe2o3_load_restart(filename,geometry)
% Validate complete histories before selecting the last saved frame.
before=fe2o3_sha256(filename);
data=load(filename);
if isfield(data,'fe2o3_result')
    result=data.fe2o3_result;
    assert(isstruct(result) && isscalar(result) && isfield(result,'schema_version') && ...
        strcmp(result.schema_version,'1.0') && isfield(result,'geometry') && isfield(result,'trajectory'),...
        'Fe2O3:InvalidRestart','Unsupported result schema.');
    for name=fieldnames(geometry)'
        assert(isfield(result.geometry,name{1}) && isequal(result.geometry.(name{1}),geometry.(name{1})),...
            'Fe2O3:InvalidRestart','Restart mapping differs: %s.',name{1});
    end
    data=result.trajectory; format='fe2o3_result/1.0';
else
    format='legacy mmx/mmy/mmz; source provenance unavailable';
end
shape=size(geometry.atomtype_); occupied=logical(geometry.atomtype_(:));
assert(all(isfield(data,{'mmx','mmy','mmz'})),'Fe2O3:InvalidRestart','Missing spin components.');
components={'mmx','mmy','mmz'}; count=size(data.mmx,4);
squares=zeros(prod(shape),count);
for index=1:3
    value=data.(components{index});
    assert(isa(value,'double') && isreal(value) && ~issparse(value) && ...
        ndims(value)<=4 && isequal([size(value,1),size(value,2),size(value,3)],shape) && ...
        size(value,4)==count && count>=1 && all(isfinite(value),'all'),...
        'Fe2O3:InvalidRestart','Invalid restart component shape/type/values.');
    flat=reshape(value,[],count);
    assert(all(flat(~occupied,:)==0,'all'),'Fe2O3:InvalidRestart','Nonzero empty restart sites.');
    squares=squares+flat.^2;
end
assert(max(abs(sqrt(squares(occupied,:))-1),[],'all')<=2e-12,...
    'Fe2O3:InvalidRestart','Restart occupied spins must be unit length; no implicit normalization.');
assert(strcmp(before,fe2o3_sha256(filename)),'Fe2O3:SourceChanged','Restart changed while loading.');
x=data.mmx(:,:,:,end); y=data.mmy(:,:,:,end); z=data.mmz(:,:,:,end);
info=struct('path',char(java.io.File(filename).getCanonicalPath()),'sha256',before,...
    'format',format,'selected_frame',count,'meaning','last saved sample; not a full solver checkpoint');
end
