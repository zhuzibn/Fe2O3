function check_output_manifest
% Exercise ownership on disposable copies, never overwrite real acceptance data.
root=fileparts(fileparts(mfilename('fullpath')));
oldpath=path; oldcwd=pwd;
cleanup=onCleanup(@() restore(oldpath,oldcwd));
addpath(root,fullfile(root,'benchmark'));
scratch=fe2o3_candidate_directory;
copy=fullfile(scratch,'source'); mkdir(copy);
copyfile(fullfile(root,'fe2o3_output_path.m'),copy);
mkdir(fullfile(copy,'benchmark','baseline'));
manifest=fullfile(copy,'benchmark','baseline_manifest.json');
fid=fopen(manifest,'w'); assert(fid>=0); fprintf(fid,'{"sentinel":"preserve"}\n'); fclose(fid);
before=fileread(manifest);
cd(copy); addpath(copy,'-begin'); clear fe2o3_output_path; rehash;
alias=fullfile(scratch,'manifest-hardlink.json');
java.nio.file.Files.createLink(java.io.File(alias).toPath(),java.io.File(manifest).toPath());
paths={manifest,fullfile('benchmark','baseline_manifest.json'),...
    fullfile(copy,'benchmark','baseline','..','baseline_manifest.json'),alias};
if ispc
    junction=fullfile(scratch,'benchmark-junction');
    [status,message]=system(sprintf('cmd /c mklink /J "%s" "%s"',junction,fullfile(copy,'benchmark')));
    assert(status==0,'Fe2O3:TestFailure','Cannot verify junction: %s',message);
    paths=[paths,{upper(manifest),fullfile(junction,'baseline_manifest.json')}];
end
for overwrite=[false,true]
    for k=1:numel(paths)
        try
            fe2o3_output_path(paths{k},'',overwrite);
            error('Fe2O3:TestFailure','Accepted protected manifest destination.');
        catch exception
            assert(strcmp(exception.identifier,'Fe2O3:ProtectedBaseline'),...
                'Fe2O3:TestFailure','Wrong rejection: %s',exception.message);
        end
    end
end
assert(strcmp(before,fileread(manifest)) && strcmp(before,fileread(alias)));
fresh=fullfile(scratch,'ordinary.mat');
assert(~isempty(fe2o3_output_path(fresh,'',false)));
fid=fopen(fresh,'w'); fprintf(fid,'ordinary existing result'); fclose(fid);
assert(~isempty(fe2o3_output_path(fresh,'',true)));
assert(strcmp(fileread(fresh),'ordinary existing result'));
fprintf('PASS acceptance manifest: direct, parent, hard-link and platform aliases; %d rejections; bytes preserved\n',2*numel(paths));
end

function restore(oldpath,oldcwd)
cd(oldcwd); path(oldpath); clear fe2o3_output_path; rehash;
end
