function destination = fe2o3_candidate_directory(requested)
% Reserve a fresh directory, resolving existing aliases before baseline guard.
here=fileparts(mfilename('fullpath'));
if nargin==0 || isempty(requested)
    parent=fullfile(here,'candidates');
    if ~isfolder(parent), mkdir(parent); end
    requested=tempname(parent);
end
destination=resolve(requested);
baseline=resolve(fullfile(here,'baseline'));
comparison=destination; protected=baseline;
if ispc, comparison=lower(comparison); protected=lower(protected); end
assert(~strcmp(comparison,protected) && ~startsWith(comparison,[protected filesep]),...
    'Fe2O3:ProtectedBaseline','Candidate destination is inside baseline storage.');
% mkdirs returns false for existing files/directories, including a raced mkdir.
assert(java.io.File(destination).mkdirs(),...
    'Fe2O3:OutputExists','Destination exists or cannot be created: %s',destination);
end

function resolved=resolve(requested)
assert((ischar(requested) && isrow(requested)) || (isstring(requested) && isscalar(requested)),...
    'Fe2O3:InvalidPath','Expected a single path.');
if ~java.io.File(char(requested)).isAbsolute(), requested=fullfile(pwd,char(requested)); end
ancestor=char(java.io.File(char(requested)).getAbsolutePath()); suffix='';
while ~isfolder(ancestor) && ~isfile(ancestor)
    [parent,name,extension]=fileparts(ancestor);
    assert(~isempty(parent) && ~strcmp(parent,ancestor),'Fe2O3:InvalidPath','Cannot resolve path.');
    suffix=fullfile([name extension],suffix); ancestor=parent;
end
realpath=java.io.File(ancestor).toPath().toRealPath(javaArray('java.nio.file.LinkOption',0));
resolved=char(java.io.File(fullfile(char(realpath.toString()),suffix)).getCanonicalPath());
end
