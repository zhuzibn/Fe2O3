mmx_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
mmy_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
mmz_=zeros(natomW,natomL,natomH, gpurun_number*final_m_savestep);
%% J1, J2, J3, J4 only shift the mx,,my,mz and set the poper boundary condition.
scalgpu=gam/(1+alp^2);%scale parameter
BDSOT=BDSOTRE;
BDSTT=BDSTTRE;        
gamatom=scalgpu*(1+alp^2);
clear ctW ctL
BFSOT=chi*BDSOT;
BFSTT=chi*BDSTT;
muigpu=mus;
ct3run=round((runtime)/gpusave);
ct3=1;
while ~(ct3>ct3run)
    
    mmx=zeros(natomW,natomL,natomH,gpusteps,'gpuArray');
    mmy=zeros(natomW,natomL,natomH,gpusteps,'gpuArray');
    mmz=zeros(natomW,natomL,natomH,gpusteps,'gpuArray');

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
        mmxtmp(atomtype_==2)=0;
        mmytmp(atomtype_==2)=0;
        mmztmp(atomtype_==2)=0;
        mmx(:,:,:,ct1)=mmxtmp;
        mmy(:,:,:,ct1)=mmytmp;
        mmz(:,:,:,ct1)=mmztmp;
        mmxtmpJ1=mmxtmp;mmytmpJ1=mmytmp;mmztmpJ1=mmztmp;
        mmxtmpJ2=mmxtmp;mmytmpJ2=mmytmp;mmztmpJ2=mmztmp;
        mmxtmpJ3=mmxtmp; mmytmpJ3=mmytmp;mmztmpJ3=mmztmp;
        mmxtmpJ4=mmxtmp; mmytmpJ4=mmytmp; mmztmpJ4=mmztmp;
        mmxtmpJ_1=mmxtmp;mmytmpJ_1=mmytmp;mmztmpJ_1=mmztmp;
        mmxtmpJ_2=mmxtmp;mmytmpJ_2=mmytmp;mmztmpJ_2=mmztmp;
        mmxtmpJ_3=mmxtmp; mmytmpJ_3=mmytmp;mmztmpJ_3=mmztmp;
        mmxtmpJ_4=mmxtmp; mmytmpJ_4=mmytmp; mmztmpJ_4=mmztmp;
       
        exchangej_()

       mmxtmpJ_1(atomtype_==2)=0;mmxtmpJ_2(atomtype_==2)=0;
       mmxtmpJ_3(atomtype_==2)=0;mmxtmpJ_4(atomtype_==2)=0;
       mmytmpJ_1(atomtype_==2)=0;mmytmpJ_2(atomtype_==2)=0;
       mmytmpJ_3(atomtype_==2)=0;mmytmpJ_4(atomtype_==2)=0;
       mmztmpJ_1(atomtype_==2)=0;mmztmpJ_2(atomtype_==2)=0;
       mmztmpJ_3(atomtype_==2)=0;mmztmpJ_4(atomtype_==2)=0;
   
        
       hex_x=-(J1.*mmxtmpJ_1+J2.*(mmxtmpJ_2)+J3.*(mmxtmpJ_3)+J4.*(mmxtmpJ_4))./muigpu;%[T]
       hex_y=-(J1.*mmytmpJ_1+J2.*(mmytmpJ_2)+J3.*(mmytmpJ_3)+J4.*(mmytmpJ_4))./muigpu;%[T]
       hex_z=-(J1.*mmztmpJ_1+J2.*(mmztmpJ_2)+J3.*(mmztmpJ_3)+J4.*(mmztmpJ_4))./muigpu;%[T]

       clear mmxtmpJ_2 mmxtmpJ_1  mmxtmpJ_3 mmxtmpJ_4
       clear mmytmpJ_2 mmytmpJ_1  mmytmpJ_3 mmytmpJ_4
       clear mmztmpJ_2 mmztmpJ_1  mmztmpJ_3 mmztmpJ_4
       
       hani_x=2*Ksim1./muigpu.*mmxtmp;
       hani_y=zeros(size(hex_x,1),size(hex_x,2),size(hex_x,3),'gpuArray');
       hani_z=zeros(size(hex_x,1),size(hex_x,2),size(hex_x,3),'gpuArray');%[T]

        mmxtmpd_nex=mmxtmp;
        mmytmpd_nex=mmytmp;
        mmztmpd_nex=mmztmp;
        mmxtmpd_pre=mmxtmp;
        mmytmpd_pre=mmytmp;
        mmztmpd_pre=mmztmp;

        mmxtmpdmi_nex=mmxtmp;
        mmytmpdmi_nex=mmytmp;
        mmztmpdmi_nex=mmztmp;
        mmxtmpdmi_pre=mmxtmp;
        mmytmpdmi_pre=mmytmp;
        mmztmpdmi_pre=mmztmp;
       dmi()
       mmxtmpdmi_nex(atomtype_==2)=0;mmytmpdmi_nex(atomtype_==2)=0;mmztmpdmi_nex(atomtype_==2)=0;
       mmxtmpdmi_pre(atomtype_==2)=0;mmytmpdmi_pre(atomtype_==2)=0;mmztmpdmi_pre(atomtype_==2)=0;

       hdmi_x=Dsim./muigpu.*(mmytmpdmi_nex-mmytmpdmi_pre);%[T]
       hdmi_y=Dsim./muigpu.*(-mmxtmpdmi_nex+mmxtmpdmi_pre);
       hdmi_z=zeros(size(hex_x,1),size(hex_x,2),size(hex_x,3),'gpuArray');

       clear mmxtmpd_nex mmytmpd_nex  mmztmpd_nex
       clear mmxtmpd_pre mmytmpd_pre  mmztmpd_pre
        hhx=hex_x+hani_x+hdmi_x+Hext(1);
        hhy=hex_y+hani_y+hdmi_y+Hext(2);
        hhz=hex_z+hani_z+hdmi_z+Hext(3);
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
        
         mmx(:,:,:,ct1+1)=sxx; mmy(:,:,:,ct1+1)=syy; mmz(:,:,:,ct1+1)=szz;
        if enablefixedge
            mmx(:,1:fixededgeL,ct1+1)=mxleft;
            mmy(:,1:fixededgeL,ct1+1)=myleft;
            mmz(:,1:fixededgeL,ct1+1)=mzleft;
            
            mmx(:,natomL-fixededgeL:end,ct1+1)=mxright;
            mmy(:,natomL-fixededgeL:end,ct1+1)=myright;
            mmz(:,natomL-fixededgeL:end,ct1+1)=mzright;
        end
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
    ct3=ct3+1;
end

clear mmx mmy mmz tmp2xn0 tmp2yn0 tmp2zn0 tmp2xn1 tmp2yn1 tmp2zn1
clear tmp2xn2 tmp2yn2 tmp2zn2
mmx=mmx_;
mmy=mmy_;
mmz=mmz_;
clear mmx_ mmy_ mmz_
t=t(1:savetstep:end);
