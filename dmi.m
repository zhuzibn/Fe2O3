%% Shell-resolved z-DMI neighbor sums [J].
% Phys. Rev. B 107, 184426(2023)
% DMI of J1, J3 and J4
% The J3/J4 neighbor lists used by mmxtmpJ3/mmxtmpJ4 are reciprocal:
% each included target-to-neighbor connection has a reverse connection
% in the same shell. For example, layer-1 green reads layer-2 purple
% through J3, and that purple atom reads the same green atom through J3.
% All included J3/J4 bonds connect adjacent layers. Reversing a bond
% therefore changes the target-layer parity and flips dmi_layer_sign,
% giving opposite directed DMI weights: W_ab = -W_ba.
for i=1:natomH
    dmi_layer_sign=1-2*mod(i,2);
    dmi_sumx(:,:,i)=dmi_layer_sign.*(Dsim1.*mmxtmpJ1(:,:,i)+Dsim3.*mmxtmpJ3(:,:,i)-Dsim4.*mmxtmpJ4(:,:,i));
    dmi_sumy(:,:,i)=dmi_layer_sign.*(Dsim1.*mmytmpJ1(:,:,i)+Dsim3.*mmytmpJ3(:,:,i)-Dsim4.*mmytmpJ4(:,:,i));
    dmi_sumz(:,:,i)=dmi_layer_sign.*(Dsim1.*mmztmpJ1(:,:,i)+Dsim3.*mmztmpJ3(:,:,i)-Dsim4.*mmztmpJ4(:,:,i));
end
clear dmi_layer_sign
