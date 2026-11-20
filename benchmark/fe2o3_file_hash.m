function value = fe2o3_file_hash(filename)
% SHA-256 of raw bytes, including source line endings.
fid=fopen(filename,'rb');
assert(fid~=-1,'Fe2O3:MissingArtifact','Cannot read %s.',filename);
cleanup=onCleanup(@() fclose(fid));
digest=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes=fread(fid,1024*1024,'*uint8'); digest.update(bytes);
end
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
