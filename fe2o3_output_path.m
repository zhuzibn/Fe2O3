function destination = fe2o3_output_path(requested,restart_file,overwrite)
% Resolve real parent aliases before enforcing output ownership.
assert((islogical(overwrite) || isnumeric(overwrite)) && isscalar(overwrite) && ...
    isreal(overwrite) && ismember(overwrite,[0,1]),...
    'Fe2O3:InvalidOutput','overwrite must be true or false.');
assert((ischar(requested) && (isrow(requested) || isempty(requested))) || ...
    (isstring(requested) && isscalar(requested) && ~ismissing(requested)),...
    'Fe2O3:InvalidOutput','output_file must be a text path.');
requested=char(requested);
if isempty(requested), destination=''; return; end
destination=resolve(requested);
root=fileparts(mfilename('fullpath'));
protected=resolve(fullfile(root,'benchmark','baseline'));
key=destination;
if ispc, key=lower(key); protected=lower(protected); end
assert(~strcmp(key,protected) && ~startsWith(key,[protected filesep]),...
    'Fe2O3:ProtectedBaseline','Output is inside accepted-reference storage.');
if ~isempty(restart_file)
    input=resolve(restart_file);
    if ispc, input=lower(input); end
    assert(~strcmp(key,input),'Fe2O3:RestartCollision','Output cannot replace its restart input.');
end
assert(isfolder(fileparts(destination)) && ~isfolder(destination),...
    'Fe2O3:InvalidOutput','Output needs an existing parent directory and a file name.');
assert(overwrite || ~isfile(destination),'Fe2O3:OutputExists','Output already exists: %s',destination);
end

function resolved=resolve(requested)
ancestor=char(java.io.File(char(requested)).getAbsolutePath()); suffix='';
while ~isfolder(ancestor) && ~isfile(ancestor)
    [parent,name,extension]=fileparts(ancestor);
    assert(~isempty(parent) && ~strcmp(parent,ancestor),'Fe2O3:InvalidOutput','Cannot resolve path.');
    suffix=fullfile([name extension],suffix); ancestor=parent;
end
realpath=java.io.File(ancestor).toPath().toRealPath(javaArray('java.nio.file.LinkOption',0));
resolved=char(java.io.File(fullfile(char(realpath.toString()),suffix)).getCanonicalPath());
end
