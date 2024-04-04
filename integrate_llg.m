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
        mmxtmp=mmxtmp.*atomtype_s;
        mmytmp=mmytmp.*atomtype_s;
        mmztmp=mmztmp.*atomtype_s;
        mmx(:,:,:,ct1)=mmxtmp;
        mmy(:,:,:,ct1)=mmytmp;
        mmz(:,:,:,ct1)=mmztmp;
        mmxtmpJ1=mmxtmp;mmytmpJ1=mmytmp;mmztmpJ1=mmztmp;
        mmxtmpJ2=mmxtmp;mmytmpJ2=mmytmp;mmztmpJ2=mmztmp;
        mmxtmpJ3=mmxtmp; mmytmpJ3=mmytmp;mmztmpJ3=mmztmp;
        mmxtmpJ4=mmxtmp; mmytmpJ4=mmytmp; mmztmpJ4=mmztmp;
        
     %%directly move in z direction
     %%formation mmx_(z(direction)n(negative)1)_pre
 

     %% 10 ways to move/ combine up and down=20ways
     %%The matrix below is used to cauclate the DMI,J3,J4
     % mmx_0p1_pre: 0p1 means shift the matrix [0,1]
     % pre means shift the matrix [0,0,-1] 
     mmx_0p1=circshift(mmxtmp,[0,1]);
     mmy_0p1=circshift(mmytmp,[0,1]);
     mmz_0p1=circshift(mmztmp,[0,1]);

     mmx_p1n1=circshift(mmxtmp,[1,-1]);
     mmy_p1n1=circshift(mmytmp,[1,-1]);
     mmz_p1n1=circshift(mmztmp,[1,-1]);
 
     mmx_n1n1=circshift(mmxtmp,[-1,-1]);
     mmy_n1n1=circshift(mmytmp,[-1,-1]);
     mmz_n1n1=circshift(mmztmp,[-1,-1]);

     mmx_0n1=circshift(mmxtmp,[0,-1]);
     mmy_0n1=circshift(mmytmp,[0,-1]);
     mmz_0n1=circshift(mmztmp,[0,-1]);

     mmx_n1p1=circshift(mmxtmp,[-1,1]);
     mmy_n1p1=circshift(mmytmp,[-1,1]);
     mmz_n1p1=circshift(mmztmp,[-1,1]); 

     mmx_0n2=circshift(mmxtmp,[0,-2]);
     mmy_0n2=circshift(mmytmp,[0,-2]);
     mmz_0n2=circshift(mmztmp,[0,-2]);

     mmx_p10=circshift(mmxtmp,[1,0]);
     mmy_p10=circshift(mmytmp,[1,0]);
     mmz_p10=circshift(mmztmp,[1,0]);

     mmx_n10=circshift(mmxtmp,[-1,0]);
     mmy_n10=circshift(mmytmp,[-1,0]);
     mmz_n10=circshift(mmztmp,[-1,0]);

     mmx_p1p1=circshift(mmxtmp,[1,1]);
     mmy_p1p1=circshift(mmytmp,[1,1]);
     mmz_p1p1=circshift(mmztmp,[1,1]);

     mmx_0p2=circshift(mmxtmp,[0,2]);
     mmy_0p2=circshift(mmytmp,[0,2]);
     mmz_0p2=circshift(mmztmp,[0,2]);
       
       exchangej_()

       mmxtmpJ1=mmxtmpJ1.*atomtype_s;
       mmxtmpJ2=mmxtmpJ2.*atomtype_s;
       mmxtmpJ3=mmxtmpJ3.*atomtype_s;
       mmxtmpJ4=mmxtmpJ4.*atomtype_s;
       mmytmpJ1=mmytmpJ1.*atomtype_s;
       mmytmpJ2=mmytmpJ2.*atomtype_s;
       mmytmpJ3=mmytmpJ3.*atomtype_s;
       mmytmpJ4=mmytmpJ4.*atomtype_s;
       mmztmpJ1=mmztmpJ1.*atomtype_s;
       mmztmpJ2=mmztmpJ2.*atomtype_s;
       mmztmpJ3=mmztmpJ3.*atomtype_s;
       mmztmpJ4=mmztmpJ4.*atomtype_s;

       hex_x=-(J1.*mmxtmpJ1+J2.*(mmxtmpJ2)+J3.*(mmxtmpJ3)+J4.*(mmxtmpJ4))./muigpu;%[T]
       hex_y=-(J1.*mmytmpJ1+J2.*(mmytmpJ2)+J3.*(mmytmpJ3)+J4.*(mmytmpJ4))./muigpu;%[T]
       hex_z=-(J1.*mmztmpJ1+J2.*(mmztmpJ2)+J3.*(mmztmpJ3)+J4.*(mmztmpJ4))./muigpu;%[T]


       
       hani_x=zeros(size(hex_x,1),size(hex_x,2),size(hex_x,3),'gpuArray');%anisotropy
       hani_y=zeros(size(hex_x,1),size(hex_x,2),size(hex_x,3),'gpuArray');
       hani_z=2*Ksim1./muigpu.*mmztmp+4*Ksim2./muigpu.*mmztmp.^3;%[T]

        mmxtmpd_nex=mmxtmp;
        mmytmpd_nex=mmytmp;
        mmztmpd_nex=mmztmp;
        mmxtmpd_pre=mmxtmp;
        mmytmpd_pre=mmytmp;
        mmztmpd_pre=mmztmp;

      
       dmi()
       mmxtmpd_nex=mmxtmpd_nex.*atomtype_s;
       mmytmpd_nex=mmytmpd_nex.*atomtype_s;
       mmztmpd_nex=mmztmpd_nex.*atomtype_s;
       mmxtmpd_pre=mmxtmpd_pre.*atomtype_s;
       mmytmpd_pre=mmytmpd_pre.*atomtype_s;
       mmztmpd_pre=mmztmpd_pre.*atomtype_s;
  
       hdmi_x=Dsim./muigpu.*(mmytmpd_nex-mmytmpd_pre);%[T]
       hdmi_y=Dsim./muigpu.*(-mmxtmpd_nex+mmxtmpd_pre);
       hdmi_z=zeros(size(hex_x,1),size(hex_x,2),size(hex_x,3),'gpuArray');


      if ct1==1||fix(ct1/50)==ct1/50
        hdipo_x=zeros(natomW,natomL,natomH,'gpuArray');
        hdipo_y=zeros(natomW,natomL,natomH,'gpuArray');
        hdipo_z=zeros(natomW,natomL,natomH,'gpuArray');

        dipole_();

        hdipo_x=hdipo_x.*atomtype_s;
        hdipo_y=hdipo_y.*atomtype_s;
        hdipo_z=hdipo_z.*atomtype_s;
      else
      end

        hhx=hex_x+hani_x+hdmi_x+Hext(1)+hdipo_x;
        hhy=hex_y+hani_y+hdmi_y+Hext(2)+hdipo_y;
        hhz=hex_z+hani_z+hdmi_z+Hext(3)+hdipo_z;

        mmxtmp=mmxtmp+ato_s;
        mmytmp=mmytmp+ato_s;
        mmztmp=mmztmp+ato_s;
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
mmx=mmx_.*atomtype_s;
mmy=mmy_.*atomtype_s;
mmz=mmz_.*atomtype_s;
clear mmx_ mmy_ mmz_
t=t(1:savetstep:end);
