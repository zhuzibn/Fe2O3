function value = fe2o3_sha256(filename)
% SHA-256 of the exact file bytes.
fid=fopen(filename,'rb');
assert(fid~=-1,'Fe2O3:MissingArtifact','Cannot read %s.',filename);
cleanup=onCleanup(@() fclose(fid));
digest=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    digest.update(fread(fid,1024*1024,'*uint8'));
end
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
