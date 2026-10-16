%reduce data size （optional)
% mmx_show=zeros(natomW,natomL,saven);
% mmy_show=zeros(natomW,natomL,saven);
% mmz_show=zeros(natomW,natomL,saven);
mmx_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
mmy_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
mmz_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
% J1 matrix
mmxtmpJ1=zeros(natomW,natomL,natomH);
mmytmpJ1=zeros(natomW,natomL,natomH);
mmztmpJ1=zeros(natomW,natomL,natomH);
% J2 matrix
mmxtmpJ2=zeros(natomW,natomL,natomH);
mmytmpJ2=zeros(natomW,natomL,natomH);
mmztmpJ2=zeros(natomW,natomL,natomH);
% J3 matrix
mmxtmpJ3=zeros(natomW,natomL,natomH);
mmytmpJ3=zeros(natomW,natomL,natomH);
mmztmpJ3=zeros(natomW,natomL,natomH);
% J4 matrix
mmxtmpJ4=zeros(natomW,natomL,natomH);
mmytmpJ4=zeros(natomW,natomL,natomH);
mmztmpJ4=zeros(natomW,natomL,natomH);
% DMI net layer matrix
mmxtmpd_nex=zeros(natomW,natomL,natomH);
mmytmpd_nex=zeros(natomW,natomL,natomH);
mmztmpd_nex=zeros(natomW,natomL,natomH);
mmxtmpd_pre=zeros(natomW,natomL,natomH);
mmytmpd_pre=zeros(natomW,natomL,natomH);
mmztmpd_pre=zeros(natomW,natomL,natomH);
%
hdmi_x=zeros(natomW,natomL,natomH);
hdmi_y=zeros(natomW,natomL,natomH);
hdmi_z=zeros(natomW,natomL,natomH);
hdmi_allx=zeros(natomW,natomL,natomH,gpurun_number*final_m_savestep);
hdmi_ally=zeros(natomW,natomL,natomH,gpurun_number*final_m_savestep);
hdmi_allz=zeros(natomW,natomL,natomH,gpurun_number*final_m_savestep);
if dipolemode
hdipolex1=zeros(natomW,natomL,natomH);
hdipolex2=zeros(natomW,natomL,natomH);
hdipolex3=zeros(natomW,natomL,natomH);
hdipolex4=zeros(natomW,natomL,natomH);
hdipolex5=zeros(natomW,natomL,natomH);
hdipolex6=zeros(natomW,natomL,natomH);
hdipoley1=zeros(natomW,natomL,natomH);
hdipoley2=zeros(natomW,natomL,natomH);
hdipoley3=zeros(natomW,natomL,natomH);
hdipoley4=zeros(natomW,natomL,natomH);
hdipoley5=zeros(natomW,natomL,natomH);
hdipoley6=zeros(natomW,natomL,natomH);
hdipolez1=zeros(natomW,natomL,natomH);
hdipolez2=zeros(natomW,natomL,natomH);
hdipolez3=zeros(natomW,natomL,natomH);
hdipolez4=zeros(natomW,natomL,natomH);
hdipolez5=zeros(natomW,natomL,natomH);
hdipolez6=zeros(natomW,natomL,natomH);
end
hdipo_x=zeros(natomW,natomL,natomH);
hdipo_y=zeros(natomW,natomL,natomH);
hdipo_z=zeros(natomW,natomL,natomH);
%%parameter
scalgpu=gam/(1+alp^2);%scale parameter
BDSOT=BDSOTRE;
BDSTT=BDSTTRE;
gamatom=scalgpu*(1+alp^2);
BFSOT=chi*BDSOT;
BFSTT=chi*BDSTT;
muigpu=mus;
ct3run=round((runtime)/gpusave);
ct3=1;
dipole_tstep=50;%decrease the dipole filed step
while ~(ct3>ct3run)

    mmx=zeros(natomW,natomL,natomH,gpusteps);
    mmy=zeros(natomW,natomL,natomH,gpusteps);
    mmz=zeros(natomW,natomL,natomH,gpusteps);
    hdmi_tmpx=zeros(natomW,natomL,natomH,gpusteps);
    hdmi_tmpy=zeros(natomW,natomL,natomH,gpusteps);
    hdmi_tmpz=zeros(natomW,natomL,natomH,gpusteps);

    if ~(ct3==1)
        mmx(:,:,:,1)=tmp2xn0;mmy(:,:,:,1)=tmp2yn0;mmz(:,:,:,1)=tmp2zn0;
    else
        mmx(:,:,:,1)=mx_init;mmy(:,:,:,1)=my_init;mmz(:,:,:,1)=mz_init;
    end
    clear tmpx tmpy tmpz
    ct1=1; %count 1
    while ct1<gpusteps+1
        mmxtmp=mmx(:,:,:,ct1);
        mmytmp=mmy(:,:,:,ct1);
        mmztmp=mmz(:,:,:,ct1);
        mmxtmp=mmxtmp.*atomtype_;
        mmytmp=mmytmp.*atomtype_;
        mmztmp=mmztmp.*atomtype_;
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
        %% calcuation for anisotropy field

      
        sintheta=sin(acos(mmztmp));
        cosphi=mmxtmp./(sintheta+(sintheta==0));%ANI-01: avoid 0/0 at the +-z poles
        sinphi=mmytmp./(sintheta+(sintheta==0));
        costheta=mmztmp;

        hk_x=hk6.*(6.*sintheta.^5.*2.*cosphi.*sinphi.*(2.*sinphi+1).*(2.*sinphi-1).*(-3+4.*sinphi.^2).*sinphi);
        hk_y=hk6.*(-6.*sintheta.^5*2.*cosphi.*sinphi.*(2.*sinphi+1).*(2.*sinphi-1).*(-3+4.*sinphi.^2).*cosphi);
        hk_z=hk6.*(-6.*costheta.*sintheta.^4.*(-1+2.*cosphi.^2).*(16.*cosphi.^4-16.*cosphi.^2+1));
        mask=~(abs(mmztmp)==1);
        hk_x2=mask.*hk_x;
        hk_y2=mask.*hk_y;
        hk_z2=mask.*hk_z;
        hani_x=hk_x2./muigpu;
        hani_y=hk_y2./muigpu;
        hani_z=2*Ksim1./muigpu.*mmztmp+hk_z2./muigpu;%+4*Ksim2./muigpu.*mmztmp.^3;%[T]
       %% calcuation for DMI field
        dmi()
        mmxtmpd_nex=mmxtmpd_nex.*atomtype_;
        mmytmpd_nex=mmytmpd_nex.*atomtype_;
        mmztmpd_nex=mmztmpd_nex.*atomtype_;
        mmxtmpd_pre=mmxtmpd_pre.*atomtype_;
        mmytmpd_pre=mmytmpd_pre.*atomtype_;
        mmztmpd_pre=mmztmpd_pre.*atomtype_;

        hdmi_x=Dsim./muigpu.*(mmytmpd_nex+mmytmpd_pre);
        hdmi_y=-Dsim./muigpu.*(mmxtmpd_nex+mmxtmpd_pre);
        hdmi_z=zeros(size(hex_x,1),size(hex_x,2),size(hex_x,3));
        

        if dipolemode
            if mod(ct1,dipole_tstep)==1
                dipole_();
                hdipo_x=hdipo_x.*atomtype_;
                hdipo_y=hdipo_y.*atomtype_;
                hdipo_z=hdipo_z.*atomtype_;
            end
        end 

        hhx=hex_x+hani_x+hdmi_x+Hext(1)+hdipo_x;
        hhy=hex_y+hani_y+hdmi_y+Hext(2)+hdipo_y;
        hhz=hex_z+hani_z+hdmi_z+Hext(3)+hdipo_z;

        mmxtmp=mmxtmp+ato_s;
        mmytmp=mmytmp+ato_s;
        mmztmp=mmztmp+ato_s;
        mmxtmp=gpuArray(mmxtmp);
        mmytmp=gpuArray(mmytmp);
        mmztmp=gpuArray(mmztmp);
        if rk4==2%4th predictor-corrector
            if ct3==1 && ~(ct1>3)
                [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                    tstep,hhx,hhy,hhz);
            else
                [sxx,syy,szz]=arrayfun(@atomgpupc4,tmpxn0,tmpyn0,tmpzn0,...
                    tmpxn1,tmpyn1,tmpzn1,tmpxn2,tmpyn2,tmpzn2,tmpxn3,tmpyn3,tmpzn3,...
                    scalgpu,alp,tstep,hhx,hhy,hhz);
            end
        elseif rk4==1 %rk4
            [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
                psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
                tstep,hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);

        else%heun
            [sxx,syy,szz]=arrayfun(@atomgpu,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                tstep,hhx,hhy,hhz);%
        end

        mmx(:,:,:,ct1+1)=gather(sxx); mmy(:,:,:,ct1+1)=gather(syy); mmz(:,:,:,ct1+1)=gather(szz);
        hdmi_tmpx(:,:,:,ct1)=hdmi_x;
        hdmi_tmpy(:,:,:,ct1)=hdmi_y;
        hdmi_tmpz(:,:,:,ct1)=hdmi_z;
        ct1=ct1+1;
        if ~(ct3==1 && ~(ct1>3)) && ct1>3
            tmpxn0=mmx(:,:,:,ct1);tmpyn0=mmy(:,:,:,ct1);tmpzn0=mmz(:,:,:,ct1);
            tmpxn1=mmx(:,:,:,ct1-1);tmpyn1=mmy(:,:,:,ct1-1);tmpzn1=mmz(:,:,:,ct1-1);
            tmpxn2=mmx(:,:,:,ct1-2);tmpyn2=mmy(:,:,:,ct1-2);tmpzn2=mmz(:,:,:,ct1-2);
            tmpxn3=mmx(:,:,:,ct1-3);tmpyn3=mmy(:,:,:,ct1-3);tmpzn3=mmz(:,:,:,ct1-3);
        elseif ~(ct3==1 && ~(ct1>3)) && ct1==2
            tmpxn0=mmx(:,:,:,ct1);tmpyn0=mmy(:,:,:,ct1);tmpzn0=mmz(:,:,:,ct1);
            tmpxn1=tmp2xn0;tmpyn1=tmp2yn0;tmpzn1=tmp2zn0;
            tmpxn2=tmp2xn1;tmpyn2=tmp2yn1;tmpzn2=tmp2zn1;
            tmpxn3=tmp2xn2;tmpyn3=tmp2yn2;tmpzn3=tmp2zn2;
        elseif ~(ct3==1 && ~(ct1>3)) && ct1==3
            tmpxn0=mmx(:,:,:,ct1);tmpyn0=mmy(:,:,:,ct1);tmpzn0=mmz(:,:,:,ct1);
            tmpxn1=mmx(:,:,:,ct1-1);tmpyn1=mmy(:,:,:,ct1-1);tmpzn1=mmz(:,:,:,ct1-1);
            tmpxn2=tmp2xn0;tmpyn2=tmp2yn0;tmpzn2=tmp2zn0;
            tmpxn3=tmp2xn1;tmpyn3=tmp2yn1;tmpzn3=tmp2zn1;
        end
    end
    tmp2xn0=mmx(:,:,:,end);tmp2yn0=mmy(:,:,:,end);tmp2zn0=mmz(:,:,:,end);
    tmp2xn1=mmx(:,:,:,end-1);tmp2yn1=mmy(:,:,:,end-1);tmp2zn1=mmz(:,:,:,end-1);
    tmp2xn2=mmx(:,:,:,end-2);tmp2yn2=mmy(:,:,:,end-2);tmp2zn2=mmz(:,:,:,end-2);
    mmx_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmx(:,:,:,1:savetstep:gpusteps));
    mmy_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmy(:,:,:,1:savetstep:gpusteps));
    mmz_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmz(:,:,:,1:savetstep:gpusteps));
    hdmi_allx(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=hdmi_tmpx(:,:,:,1:savetstep:gpusteps);
    hdmi_ally(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=hdmi_tmpy(:,:,:,1:savetstep:gpusteps);
    hdmi_allz(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=hdmi_tmpz(:,:,:,1:savetstep:gpusteps);
    ct3=ct3+1;
end

clear mmx mmy mmz tmp2xn0 tmp2yn0 tmp2zn0 tmp2xn1 tmp2yn1 tmp2zn1
clear tmp2xn2 tmp2yn2 tmp2zn2
clear tmpxn0 tmpxn1 tmpxn2 tmpxn3
clear tmpyn0 tmpyn1 tmpyn2 tmpyn3
clear tmpzn0 tmpzn1 tmpzn2 tmpzn3
clear mmxtmpJ1 mmxtmpJ2 mmxtmpJ3 mmxtmpJ4
clear mmytmpJ1 mmytmpJ2 mmytmpJ3 mmytmpJ4
clear mmztmpJ1 mmztmpJ2 mmztmpJ3 mmztmpJ4
clear mmxtmpd_nex mmytmpd_nex mmztmpd_nex
clear mmxtmpd_pre mmytmpd_pre mmztmpd_pre
clear hdmi_tmpx hdmi_tmpy hdmi_tmpz

mmx=mmx_.*atomtype_;
mmy=mmy_.*atomtype_;
mmz=mmz_.*atomtype_;
hdmi_allx=hdmi_allx.*atomtype_;
hdmi_ally=hdmi_ally.*atomtype_;
hdmi_allz=hdmi_allz.*atomtype_;
clear mmx_ mmy_ mmz_
