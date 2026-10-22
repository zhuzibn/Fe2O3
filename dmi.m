%% Shell-resolved z-DMI neighbor sums [J].
% exchangej_ must first evaluate the current-state J3/J4 neighbor sums.
% A1/A2 DMI is omitted; their exchange interactions are unchanged.
% W = diag((-1)^layer) * (Dsim3*A3 - Dsim4*A4).
% This retains the audited directed-bond signs of the former operator.
% Array layer increases opposite to physical z. Surface/edge omissions are
% inherited from exchangej_; no uniform-sublattice approximation is used.
% Energy: E_D = -x'*W*y. Field: B_D = (W*y,-W*x,0)/mus.
for i=1:natomH
    dmi_layer_sign=1-2*mod(i,2);
    dmi_sumx(:,:,i)=dmi_layer_sign.*(Dsim3.*mmxtmpJ3(:,:,i)-Dsim4.*mmxtmpJ4(:,:,i));
    dmi_sumy(:,:,i)=dmi_layer_sign.*(Dsim3.*mmytmpJ3(:,:,i)-Dsim4.*mmytmpJ4(:,:,i));
    dmi_sumz(:,:,i)=dmi_layer_sign.*(Dsim3.*mmztmpJ3(:,:,i)-Dsim4.*mmztmpJ4(:,:,i));
end
clear dmi_layer_sign
