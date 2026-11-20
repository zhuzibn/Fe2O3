function fe2o3_save_result(fe2o3_result,destination,restart_file,overwrite)
% Publish a complete MAT file. Exclusive hard-link creation prevents races.
destination=fe2o3_output_path(destination,restart_file,overwrite);
assert(~isempty(destination),'Fe2O3:InvalidOutput','Saving needs a destination.');
temporary=[tempname(fileparts(destination)) '.mat'];
assert(java.io.File(temporary).createNewFile(),'Fe2O3:InvalidOutput','Cannot reserve temporary output.');
cleanup=onCleanup(@() remove_temporary(temporary));
save(temporary,'fe2o3_result','-v7.3');
fe2o3_output_path(destination,restart_file,overwrite); % Recheck after serialization.
source=java.io.File(temporary).toPath(); target=java.io.File(destination).toPath();
if overwrite
    options=javaArray('java.nio.file.CopyOption',2);
    options(1)=java.nio.file.StandardCopyOption.REPLACE_EXISTING;
    options(2)=java.nio.file.StandardCopyOption.ATOMIC_MOVE;
    java.nio.file.Files.move(source,target,options);
else
    try
        java.nio.file.Files.createLink(target,source);
    catch exception
        if isfile(destination)
            error('Fe2O3:OutputExists','Output appeared while saving: %s',destination);
        end
        rethrow(exception); % No unsafe fallback on filesystems without hard links.
    end
end
end

function remove_temporary(filename)
% Only this invocation's reserved temporary file; moved files no longer exist.
if isfile(filename), delete(filename); end
end
