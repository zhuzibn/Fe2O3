function run_cpu_checks
% Focused current-code checks; no production trajectory or MAT output.
root=fileparts(fileparts(mfilename('fullpath')));
oldpath=path; cleanup=onCleanup(@() path(oldpath)); addpath(root);
fprintf('CPU checks: MATLAB %s; %s\n',version,computer);
cores=[5,9,6;4,8,3;6,10,9];
expected_counts=[114,328,276,544;32,106,72,140;240,654,584,1160];
for row=1:size(cores,1)
    p=fe2o3_test_fixture(cores(row,:));
    [mask,colors,A]=fe2o3_expected_geometry(cores(row,:));
    assert(isequal(mask,p.atomtype_));
    masks={p.atomtype_layer1gr,p.atomtype_layer1or,p.atomtype_layer2blue,...
        p.atomtype_layer2p,p.atomtype_layer3red,p.atomtype_layer3black};
    for color=1:6
        k=ceil(color/2);
        assert(isequal(logical(masks{color}),colors(:,:,k)==color));
    end
    occupied=find(mask); n=numel(occupied);
    assert(isequal(cellfun(@nnz,A),expected_counts(row,:)));
    for shell=1:4
        assert(isequal(A{shell},A{shell}'));
        assert(all(diag(A{shell})==0));
    end
    for probe=1:n
        s=zeros(numel(mask),3); s(occupied(probe),:)=[1,2,-3];
        f=fe2o3_test_field(p,s);
        for shell=1:4
            assert(isequal(f.shells{shell}(occupied,:),A{shell}(:,probe)*[1,2,-3]));
        end
    end
    [~,~,k]=ind2sub(size(mask),occupied);
    W=(1-2*mod(k,2)).*(p.Dsim3*A{3}-p.Dsim4*A{4});
    assert(nnz(W+W')==0);
    s=[p.mx_init(:),p.my_init(:),p.mz_init(:)];
    assert(all(s(~mask(:),:)==0,'all'));
    assert(max(abs(s(occupied,1)-(-1).^(k+1)))<1e-14);
    assert(max(abs(s(occupied,2:3)),[],'all')<1e-14);
    fprintf('PASS topology %dx%dx%d: %d sites; J entries %d/%d/%d/%d; DMI %d\n',...
        cores(row,:),n,cellfun(@nnz,A),nnz(W));
end
p=fe2o3_test_fixture([5,9,6]);
[mask,colors,A]=fe2o3_expected_geometry(p.core); occupied=find(mask);
n=numel(occupied); s=zeros(numel(mask),3);
v=reshape(sin(1:3*n),n,3); v=v./vecnorm(v,2,2); s(occupied,:)=v;
q=p; q.mus=1; q.J1=0; q.J2=0; q.J3=0; q.J4=0;
q.Ksim1=0; q.hk6=0; q.Dsim1=0; q.Dsim3=0; q.Dsim4=0;
q.Hext=[0,0,0]; q.dipolemode=0;
for shell=1:4
    r=q; r.(sprintf('J%d',shell))=(-1)^shell*0.2;
    f=fe2o3_test_field(r,s); expected=-r.(sprintf('J%d',shell))*A{shell}*v;
    assert(max(abs(f.b(occupied,:)-expected),[],'all')<2e-12);
    energy=@(u) 0.5*r.(sprintf('J%d',shell))*sum(u.*(A{shell}*u),'all');
    check_derivative(energy,v,f.b(occupied,:));
end
r=q; r.Ksim1=-0.7; r.Hext=[0.2,-0.3,0.4];
f=fe2o3_test_field(r,s);
expected=repmat(r.Hext,n,1); expected(:,3)=expected(:,3)+2*r.Ksim1*v(:,3);
assert(max(abs(f.b(occupied,:)-expected),[],'all')<2e-12);
check_derivative(@(u) -r.Ksim1*sum(u(:,3).^2)-sum(u.*r.Hext,'all'),v,expected);
[~,~,k]=ind2sub(size(mask),occupied);
for weights=[0.7,0.4;0.7,0;0,0.4;0,0]'
    r=q; r.Dsim3=weights(1); r.Dsim4=weights(2);
    W=(1-2*mod(k,2)).*(r.Dsim3*A{3}-r.Dsim4*A{4});
    assert(nnz(W+W')==0); f=fe2o3_test_field(r,s);
    expected=[W*v(:,2),-W*v(:,1),zeros(n,1)];
    assert(max(abs(f.dmi(occupied,:)-expected),[],'all')<2e-12);
    check_derivative(@(u) -u(:,1)'*W*u(:,2),v,expected);
end
fprintf('PASS fields: isolated exchange, K1, external, D3/D4/off; Cartesian energy derivatives\n');
% hk6 arithmetic on poles, off-pole/nonunit states; no invented hk6 energy.
r=q; r.Ksim1=-0.7; r.hk6=0.03;
special=[0,0,1;0,0,-1;1,0,0;0.6,0.8,0;0.4,-0.7,1.2;0,0,0];
for test=1:size(special,1)
    u=special(test,:); state=zeros(size(s)); state(occupied,:)=repmat(u,n,1);
    f=fe2o3_test_field(r,state); expected=angular_field(u,r.Ksim1,r.hk6);
    assert(all(isfinite(f.b),'all'));
    assert(max(abs(f.b(occupied,:)-expected),[],'all')<2e-12);
    assert(all(f.b(~mask(:),:)==0,'all'));
end
% Direct retained tensor on independently indexed J2 neighbors.
r=p; r.dipolemode=1; f=fe2o3_test_field(r,s); expected=zeros(n,3);
for target=1:n
    [tr,tc,~]=ind2sub(size(mask),occupied(target));
    for source=find(A{2}(target,:))
        [sr,sc,~]=ind2sub(size(mask),occupied(source));
        if sr==tr
            rx=sign(sc-tc)*0.5067e-9; ry=0;
        else
            direction=ismember(colors(occupied(target)),[1,3,6])*2-1;
            rx=-direction*0.25335e-9; ry=-sign(sr-tr)*0.4388e-9;
        end
        bond=[rx,ry,0]; distance=0.5067e-9;
        expected(target,:)=expected(target,:)+p.mu_0*p.mus/(4*pi)*...
            (3*bond*dot(v(source,:),bond)/distance^5-v(source,:)/distance^3);
    end
end
assert(max(abs(f.dipole(occupied,:)-expected),[],'all')<2e-12);
fprintf('PASS hk6 guarded arithmetic and retained dipole tensor (not material validation)\n');
check_rhs;
fprintf('PASS all CPU checks; no simulation results or accepted references written.\n');
end

function check_derivative(energy,s,b)
h=1e-6; worst=0;
for index=1:numel(s)
    plus=s; minus=s; plus(index)=plus(index)+h; minus(index)=minus(index)-h;
    numerical=-(energy(plus)-energy(minus))/(2*h);
    worst=max(worst,abs(numerical-b(index)));
end
assert(worst<1e-8,'Energy derivative error %.12g exceeds 1e-8.',worst);
end

function b=angular_field(u,K,h6)
normu=norm(u); z=0; q=0; phi=0;
if normu~=0
    z=u(3)/normu; q=hypot(u(1),u(2))/normu; phi=atan2(u(2),u(1));
end
c=cos(phi); v=sin(phi);
b=[12*h6*q^5*c*v^2*(2*v+1)*(2*v-1)*(-3+4*v^2),...
    -12*h6*q^5*c^2*v*(2*v+1)*(2*v-1)*(-3+4*v^2),...
    -6*h6*z*q^4*(2*c^2-1)*(16*c^4-16*c^2+1)+2*K*u(3)];
end

function check_rhs
u=[0.3,-0.4,sqrt(0.75)]; b=[0.2,0.6,-0.1];
sigma=[1,2,3]/sqrt(14); polar=[-2,1,0]/sqrt(5); alpha=0.3;
channels=[zeros(1,4);diag([0.2,-0.4,0.1,0.3]);0.2,-0.4,0.1,0.3];
for amplitudes=channels'
    a=num2cell([u,sigma,polar,1/(1+alpha^2),alpha,b,amplitudes']);
    [x,y,z]=atomgpurk4(a{:}); actual=[x,y,z];
    rhs=-cross(u,b);
    for pair=1:2
        if pair==1, p=sigma; else, p=polar; end
        bd=amplitudes(2*pair-1); bf=amplitudes(2*pair);
        rhs=rhs-bd*cross(u,cross(u,p))+bf*cross(u,p);
    end
    % Solve implicit Gilbert equation independently; gamma=1.
    C=[0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];
    expected=((eye(3)-alpha*C)\rhs')';
    assert(max(abs(actual-expected))<2e-12);
end
[x,y,z]=atomgpurk4(1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,0);
assert(isequal([x,y,z],[0,1,0]));
[~,~,z]=atomgpurk4(1,0,0,0,1,0,0,0,1,1,0.2,0,0,1,0,0,0,0);
assert(z>0);
fprintf('PASS LLGS implicit-Gilbert residual, isolated/combined torques, precession/damping signs\n');
end
