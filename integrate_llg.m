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

    if ~(ct3==1)
        mmx(:,:,:,1)=tmp2xn0;mmy(:,:,:,1)=tmp2yn0;mmz(:,:,:,1)=tmp2zn0;
    else
        mmx(:,:,:,1)=mx_init;mmy(:,:,:,1)=my_init;mmz(:,:,:,1)=mz_init;
    end
    clear tmpx tmpy tmpz
    ct1=1; %count 1
    while ct1<gpusteps

        mmxtmp=mmx(:,:,:,ct1);
        mmytmp=mmy(:,:,:,ct1);
        mmztmp=mmz(:,:,:,ct1);
        mmxtmp=mmxtmp.*atomtype_;
        mmytmp=mmytmp.*atomtype_;
        mmztmp=mmztmp.*atomtype_;
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
          field_calc();
            mmxtmp=mmxtmp+ato_s;
            mmytmp=mmytmp+ato_s;
            mmztmp=mmztmp+ato_s;
            mmxtmp=gpuArray(mmxtmp);
            mmytmp=gpuArray(mmytmp);
            mmztmp=gpuArray(mmztmp);
            mmxtmp0=mmxtmp;
            mmytmp0=mmytmp;
            mmztmp0=mmztmp;
            [kk1x,kk1y,kk1z]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
                psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
                tstep,hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);
            mmxtmp=mmxtmp0+1/2*tstep*kk1x;
            mmytmp=mmytmp0+1/2*tstep*kk1y;
            mmztmp=mmztmp0+1/2*tstep*kk1z;
            mmxtmp=mmxtmp.*atomtype_;
            mmytmp=mmytmp.*atomtype_;
            mmztmp=mmztmp.*atomtype_;
            field_calc();
            mmytmp=mmytmp+ato_s;
            mmztmp=mmztmp+ato_s;
            mmxtmp=gpuArray(mmxtmp);
            mmytmp=gpuArray(mmytmp);
            mmztmp=gpuArray(mmztmp);
            [kk2x,kk2y,kk2z]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
                psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
                tstep,hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);
            mmxtmp=mmxtmp0+1/2*tstep*kk2x;
            mmytmp=mmytmp0+1/2*tstep*kk2y;
            mmztmp=mmztmp0+1/2*tstep*kk2z;

            mmxtmp=mmxtmp.*atomtype_;
            mmytmp=mmytmp.*atomtype_;
            mmztmp=mmztmp.*atomtype_;
            field_calc();
            mmytmp=mmytmp+ato_s;
            mmztmp=mmztmp+ato_s;
            mmxtmp=gpuArray(mmxtmp);
            mmytmp=gpuArray(mmytmp);
            mmztmp=gpuArray(mmztmp);
            [kk3x,kk3y,kk3z]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
                psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
                tstep,hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);
            mmxtmp=mmxtmp0+1/2*tstep*kk3x;
            mmytmp=mmytmp0+1/2*tstep*kk3y;
            mmztmp=mmztmp0+1/2*tstep*kk3z;

            mmxtmp=mmxtmp.*atomtype_;
            mmytmp=mmytmp.*atomtype_;
            mmztmp=mmztmp.*atomtype_;
            field_calc();
            mmytmp=mmytmp+ato_s;
            mmztmp=mmztmp+ato_s;
            mmxtmp=gpuArray(mmxtmp);
            mmytmp=gpuArray(mmytmp);
            mmztmp=gpuArray(mmztmp);
            [kk4x,kk4y,kk4z]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,...
                psjSHEy,psjSHEz,psjSTTx,psjSTTy,psjSTTz,scalgpu,alp,...
                tstep,hhx,hhy,hhz,BDSOT,BFSOT,BDSTT,BFSTT);
            snx=mmxtmp0+1/6*tstep*(kk1x+2*kk2x+2*kk3x+kk4x);
            sny=mmytmp0+1/6*tstep*(kk1y+2*kk2y+2*kk3y+kk4y);
            snz=mmztmp0+1/6*tstep*(kk1z+2*kk2z+2*kk3z+kk4z);
            normsn=sqrt(snx.^2+sny.^2+snz.^2);
            sxx=snx./normsn;
            syy=sny./normsn;
            szz=snz./normsn;


        else%heun
            [sxx,syy,szz]=arrayfun(@atomgpu,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                tstep,hhx,hhy,hhz);%
        end

        mmx(:,:,:,ct1+1)=sxx; mmy(:,:,:,ct1+1)=syy; mmz(:,:,:,ct1+1)=szz;

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
    mmx_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmx(:,:,:,1:savetstep:end));
    mmy_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmy(:,:,:,1:savetstep:end));
    mmz_(:,:,:,(ct3-1)*final_m_savestep+1:ct3*final_m_savestep)=gather(mmz(:,:,:,1:savetstep:end));
    %reduce data size (optional)
    % mmx_show(:,:,ct3)=gather(tmp2xn0);
    % mmy_show(:,:,ct3)=gather(tmp2yn0);
    % mmz_show(:,:,ct3)=gather(tmp2zn0);
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

mmx=mmx_.*atomtype_;
mmy=mmy_.*atomtype_;
mmz=mmz_.*atomtype_;
clear mmx_ mmy_ mmz_
t=t(1:savetstep:end);
