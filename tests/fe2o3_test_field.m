function out = fe2o3_test_field(p,s)
% Execute the actual field script on supplied CPU spins; no trajectory/save.
names=fieldnames(p);
for j=1:numel(names), eval([names{j} '=p.(names{j});']); end
muigpu=mus; dipole_refresh=true;
mmxtmp=reshape(s(:,1),size(atomtype_));
mmytmp=reshape(s(:,2),size(atomtype_));
mmztmp=reshape(s(:,3),size(atomtype_));
eval(fileread(fullfile(p.root,'field_calc.m')));
out.b=[hhx(:),hhy(:),hhz(:)];
out.exchange=[hex_x(:),hex_y(:),hex_z(:)];
out.anisotropy=[hani_x(:),hani_y(:),hani_z(:)];
out.dmi=[hdmi_x(:),hdmi_y(:),hdmi_z(:)];
out.dipole=zeros(size(s));
if dipolemode, out.dipole=[hdipo_x(:),hdipo_y(:),hdipo_z(:)]; end
out.shells=cell(1,4);
for j=1:4
    out.shells{j}=[reshape(eval(sprintf('mmxtmpJ%d',j)),[],1),...
        reshape(eval(sprintf('mmytmpJ%d',j)),[],1),...
        reshape(eval(sprintf('mmztmpJ%d',j)),[],1)];
end
end
