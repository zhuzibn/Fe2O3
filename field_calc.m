mmxtmp=mmxtmp.*atomtype_;
mmytmp=mmytmp.*atomtype_;
mmztmp=mmztmp.*atomtype_;
mmxtmpJ1=zeros(size(mmxtmp),'like',mmxtmp);
mmxtmpJ2=zeros(size(mmxtmp),'like',mmxtmp);
mmxtmpJ3=zeros(size(mmxtmp),'like',mmxtmp);
mmxtmpJ4=zeros(size(mmxtmp),'like',mmxtmp);
mmytmpJ1=zeros(size(mmxtmp),'like',mmxtmp);
mmytmpJ2=zeros(size(mmxtmp),'like',mmxtmp);
mmytmpJ3=zeros(size(mmxtmp),'like',mmxtmp);
mmytmpJ4=zeros(size(mmxtmp),'like',mmxtmp);
mmztmpJ1=zeros(size(mmxtmp),'like',mmxtmp);
mmztmpJ2=zeros(size(mmxtmp),'like',mmxtmp);
mmztmpJ3=zeros(size(mmxtmp),'like',mmxtmp);
mmztmpJ4=zeros(size(mmxtmp),'like',mmxtmp);
dmi_sumx=zeros(size(mmxtmp),'like',mmxtmp);
dmi_sumy=zeros(size(mmxtmp),'like',mmxtmp);
dmi_sumz=zeros(size(mmxtmp),'like',mmxtmp);
%% 10 ways to move/ combine up and down=20ways
% The matrix below is used to cauclate the DMI,J3,J4
% mmx_r1: r1 means shift the matrix [0,1]
mmx_r1=circshift(mmxtmp,[0,1]); %[0 1]  shift right one position
mmy_r1=circshift(mmytmp,[0,1]);
mmz_r1=circshift(mmztmp,[0,1]);

mmx_d1l1=circshift(mmxtmp,[1,-1]);
mmy_d1l1=circshift(mmytmp,[1,-1]);
mmz_d1l1=circshift(mmztmp,[1,-1]);

mmx_u1l1=circshift(mmxtmp,[-1,-1]);
mmy_u1l1=circshift(mmytmp,[-1,-1]);
mmz_u1l1=circshift(mmztmp,[-1,-1]);

mmx_l1=circshift(mmxtmp,[0,-1]);
mmy_l1=circshift(mmytmp,[0,-1]);
mmz_l1=circshift(mmztmp,[0,-1]);

mmx_u1r1=circshift(mmxtmp,[-1,1]);
mmy_u1r1=circshift(mmytmp,[-1,1]);
mmz_u1r1=circshift(mmztmp,[-1,1]);

mmx_l2=circshift(mmxtmp,[0,-2]);
mmy_l2=circshift(mmytmp,[0,-2]);
mmz_l2=circshift(mmztmp,[0,-2]);

mmx_d1=circshift(mmxtmp,[1,0]);%[1 0]  shift down one position
mmy_d1=circshift(mmytmp,[1,0]);
mmz_d1=circshift(mmztmp,[1,0]);

mmx_u1=circshift(mmxtmp,[-1,0]);
mmy_u1=circshift(mmytmp,[-1,0]);
mmz_u1=circshift(mmztmp,[-1,0]);

mmx_d1r1=circshift(mmxtmp,[1,1]);
mmy_d1r1=circshift(mmytmp,[1,1]);
mmz_d1r1=circshift(mmztmp,[1,1]);

mmx_r2=circshift(mmxtmp,[0,2]);
mmy_r2=circshift(mmytmp,[0,2]);
mmz_r2=circshift(mmztmp,[0,2]);
%% calcuation for exchange interaction
exchangej_()

mmxtmpJ1=mmxtmpJ1.*atomtype_;
mmytmpJ1=mmytmpJ1.*atomtype_;
mmztmpJ1=mmztmpJ1.*atomtype_;

mmxtmpJ2=mmxtmpJ2.*atomtype_;
mmytmpJ2=mmytmpJ2.*atomtype_;
mmztmpJ2=mmztmpJ2.*atomtype_;

mmxtmpJ3=mmxtmpJ3.*atomtype_;
mmytmpJ3=mmytmpJ3.*atomtype_;
mmztmpJ3=mmztmpJ3.*atomtype_;

mmxtmpJ4=mmxtmpJ4.*atomtype_;
mmytmpJ4=mmytmpJ4.*atomtype_;
mmztmpJ4=mmztmpJ4.*atomtype_;

hex_x=-(J1.*mmxtmpJ1+J2.*(mmxtmpJ2)+J3.*(mmxtmpJ3)+J4.*(mmxtmpJ4))./muigpu;%[T]
hex_y=-(J1.*mmytmpJ1+J2.*(mmytmpJ2)+J3.*(mmytmpJ3)+J4.*(mmytmpJ4))./muigpu;%[T]
hex_z=-(J1.*mmztmpJ1+J2.*(mmztmpJ2)+J3.*(mmztmpJ3)+J4.*(mmztmpJ4))./muigpu;%[T]
%% anisotropy
%in steps of rk4, mx, my, mz may exceed 1, renormalize to calculate ani
ani_norm=sqrt(mmxtmp.^2+mmytmp.^2+mmztmp.^2);
ani_denom=ani_norm+(ani_norm==0);
ani_ux=mmxtmp./ani_denom;
ani_uy=mmytmp./ani_denom;

costheta=mmztmp./ani_denom;
sintheta=hypot(ani_ux,ani_uy);
cosphi=ani_ux./(sintheta+(sintheta==0));
sinphi=ani_uy./(sintheta+(sintheta==0));

hk_x=hk6.*(6.*sintheta.^5.*2.*cosphi.*sinphi.*(2.*sinphi+1).*(2.*sinphi-1).*(-3+4.*sinphi.^2).*sinphi);
hk_y=hk6.*(-6.*sintheta.^5*2.*cosphi.*sinphi.*(2.*sinphi+1).*(2.*sinphi-1).*(-3+4.*sinphi.^2).*cosphi);
hk_z=hk6.*(-6.*costheta.*sintheta.^4.*(-1+2.*cosphi.^2).*(16.*cosphi.^4-16.*cosphi.^2+1));

mask=(sintheta~=0) & (atomtype_~=0);

hk_x2=mask.*hk_x;
hk_y2=mask.*hk_y;
hk_z2=mask.*hk_z;

hani_x=hk_x2./muigpu;
hani_y=hk_y2./muigpu;
hani_z=2*Ksim1./muigpu.*mmztmp+hk_z2./muigpu;
%% DMI
dmi()
dmi_sumx=dmi_sumx.*atomtype_;
dmi_sumy=dmi_sumy.*atomtype_;
dmi_sumz=dmi_sumz.*atomtype_;

hdmi_x=dmi_sumy./muigpu;
hdmi_y=-dmi_sumx./muigpu;
hdmi_z=zeros(size(hex_x),'like',hex_x);
%% dipole
if dipolemode
    if dipole_refresh
        hdipolex1=zeros(size(mmxtmp),'like',mmxtmp);
        hdipolex2=hdipolex1;hdipolex3=hdipolex1;
        hdipolex4=hdipolex1;hdipolex5=hdipolex1;hdipolex6=hdipolex1;
        hdipoley1=hdipolex1;hdipoley2=hdipolex1;hdipoley3=hdipolex1;
        hdipoley4=hdipolex1;hdipoley5=hdipolex1;hdipoley6=hdipolex1;
        hdipolez1=hdipolex1;hdipolez2=hdipolex1;hdipolez3=hdipolex1;
        hdipolez4=hdipolex1;hdipolez5=hdipolex1;hdipolez6=hdipolex1;
        dipole_();
        hdipo_x=hdipo_x.*atomtype_;
        hdipo_y=hdipo_y.*atomtype_;
        hdipo_z=hdipo_z.*atomtype_;
        dipole_refresh=false;
    end
else
    hdipo_x=0;hdipo_y=0;hdipo_z=0;
end
%% Total Heff
hhx=(hex_x+hani_x+hdmi_x+Hext(1)+hdipo_x).*atomtype_;
hhy=(hex_y+hani_y+hdmi_y+Hext(2)+hdipo_y).*atomtype_;
hhz=(hex_z+hani_z+hdmi_z+Hext(3)+hdipo_z).*atomtype_;
