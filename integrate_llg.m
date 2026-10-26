%reduce data size （optional)
% mmx_show=zeros(natomW,natomL,saven);
% mmy_show=zeros(natomW,natomL,saven);
% mmz_show=zeros(natomW,natomL,saven);
mmx_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
mmy_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
mmz_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
hdmi_allx=zeros(natomW,natomL,natomH,gpurun_number*final_m_savestep);
hdmi_ally=zeros(natomW,natomL,natomH,gpurun_number*final_m_savestep);
hdmi_allz=zeros(natomW,natomL,natomH,gpurun_number*final_m_savestep);
%%parameter
scalgpu=gam/(1+alp^2);%scale parameter
BDSOT=BDSOTRE;
BDSTT=BDSTTRE;
BFSOT=chi*BDSOT;
BFSTT=chi*BDSTT;
muigpu=mus;
ct3run=round((runtime)/gpusave);
ct3=1;
while ~(ct3>ct3run)

    mmx=zeros(natomW,natomL,natomH,gpusteps+1);
    mmy=zeros(natomW,natomL,natomH,gpusteps+1);
    mmz=zeros(natomW,natomL,natomH,gpusteps+1);
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
        mmxtmp0=gpuArray(mmx(:,:,:,ct1).*atomtype_);
        mmytmp0=gpuArray(mmy(:,:,:,ct1).*atomtype_);
        mmztmp0=gpuArray(mmz(:,:,:,ct1).*atomtype_);

        %Step 1
        mmxtmp=mmxtmp0;mmytmp=mmytmp0;mmztmp=mmztmp0;
        field_calc();
        hdmi_tmpx(:,:,:,ct1)=gather(hdmi_x);
        hdmi_tmpy(:,:,:,ct1)=gather(hdmi_y);
        hdmi_tmpz(:,:,:,ct1)=gather(hdmi_z);
        [kk1x,kk1y,kk1z]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
            psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
            hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);

        %Step 2
        mmxtmp=mmxtmp0+0.5*tstep*kk1x;
        mmytmp=mmytmp0+0.5*tstep*kk1y;
        mmztmp=mmztmp0+0.5*tstep*kk1z;
        field_calc();
        [kk2x,kk2y,kk2z]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
            psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
            hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);

        %Step 3
        mmxtmp=mmxtmp0+0.5*tstep*kk2x;
        mmytmp=mmytmp0+0.5*tstep*kk2y;
        mmztmp=mmztmp0+0.5*tstep*kk2z;
        field_calc();
        [kk3x,kk3y,kk3z]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
            psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
            hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);

        %Step 4
        mmxtmp=mmxtmp0+tstep*kk3x;
        mmytmp=mmytmp0+tstep*kk3y;
        mmztmp=mmztmp0+tstep*kk3z;
        field_calc();
        [kk4x,kk4y,kk4z]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
            psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
            hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);

        %Total
        sxx=mmxtmp0+tstep/6*(kk1x+2*kk2x+2*kk3x+kk4x);
        syy=mmytmp0+tstep/6*(kk1y+2*kk2y+2*kk3y+kk4y);
        szz=mmztmp0+tstep/6*(kk1z+2*kk2z+2*kk3z+kk4z);
        normsn=sqrt(sxx.^2+syy.^2+szz.^2);
        if gather(any((~isfinite(normsn) | normsn==0) & atomtype_~=0,'all'))
            error('Fe2O3:InvalidSpin','RK4 produced a nonfinite or zero occupied spin.');
        end
        normsn=normsn+(atomtype_==0);
        sxx=sxx./normsn.*atomtype_;
        syy=syy./normsn.*atomtype_;
        szz=szz./normsn.*atomtype_;
        mmx(:,:,:,ct1+1)=gather(sxx);mmy(:,:,:,ct1+1)=gather(syy);mmz(:,:,:,ct1+1)=gather(szz);
        ct1=ct1+1;

    end
    tmp2xn0=mmx(:,:,:,end);tmp2yn0=mmy(:,:,:,end);tmp2zn0=mmz(:,:,:,end);
    mmx_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmx(:,:,:,1:savetstep:gpusteps));
    mmy_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmy(:,:,:,1:savetstep:gpusteps));
    mmz_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmz(:,:,:,1:savetstep:gpusteps));
    hdmi_allx(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=hdmi_tmpx(:,:,:,1:savetstep:gpusteps);
    hdmi_ally(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=hdmi_tmpy(:,:,:,1:savetstep:gpusteps);
    hdmi_allz(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=hdmi_tmpz(:,:,:,1:savetstep:gpusteps);
    ct3=ct3+1;
end

clear mmx mmy mmz tmp2xn0 tmp2yn0 tmp2zn0
clear mmxtmp0 mmytmp0 mmztmp0 normsn
clear kk1x kk1y kk1z kk2x kk2y kk2z kk3x kk3y kk3z kk4x kk4y kk4z
clear ani_norm ani_denom ani_ux ani_uy
clear mmxtmpJ1 mmxtmpJ2 mmxtmpJ3 mmxtmpJ4
clear mmytmpJ1 mmytmpJ2 mmytmpJ3 mmytmpJ4
clear mmztmpJ1 mmztmpJ2 mmztmpJ3 mmztmpJ4
clear dmi_sumx dmi_sumy dmi_sumz
clear hdmi_tmpx hdmi_tmpy hdmi_tmpz

mmx=mmx_.*atomtype_;
mmy=mmy_.*atomtype_;
mmz=mmz_.*atomtype_;
hdmi_allx=hdmi_allx.*atomtype_;
hdmi_ally=hdmi_ally.*atomtype_;
hdmi_allz=hdmi_allz.*atomtype_;
clear mmx_ mmy_ mmz_
