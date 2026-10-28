%reduce data size （optional)
% mmx_show=zeros(natomW,natomL,saven);
% mmy_show=zeros(natomW,natomL,saven);
% mmz_show=zeros(natomW,natomL,saven);
if ~isscalar(rk4) || ~ismember(rk4,[0,1,2])
    error('Fe2O3:InvalidSolver','rk4 must be 0 (Heun), 1 (RK4), or 2 (PC4).');
end
if rk4~=1 && any([BDSOTRE,chi*BDSOTRE,BDSTTRE,chi*BDSTTRE]~=0)
    error('Fe2O3:UnsupportedTorque','SOT/STT torques require rk4=1.');
end
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
dipole_tstep=50;
if rk4==2,pc_history_count=0;end
mx_current=gpuArray(mx_init);
my_current=gpuArray(my_init);
mz_current=gpuArray(mz_init);
while ~(ct3>ct3run)

    ct1=1; %count 1
    while ct1<gpusteps+1
        mmxtmp0=mx_current.*atomtype_;
        mmytmp0=my_current.*atomtype_;
        mmztmp0=mz_current.*atomtype_;
        save_frame=mod(ct1-1,savetstep)==0;
        if save_frame
            sample_index=(ct3-1)*final_m_savestep+(ct1-1)/savetstep+1;
            mmx_(:,:,:,sample_index)=gather(mmxtmp0);
            mmy_(:,:,:,sample_index)=gather(mmytmp0);
            mmz_(:,:,:,sample_index)=gather(mmztmp0);
        end

        dipole_refresh=dipolemode && mod(ct1,dipole_tstep)==1;

        mmxtmp=mmxtmp0;
        mmytmp=mmytmp0;
        mmztmp=mmztmp0;

        field_calc();
        if save_frame
            hdmi_allx(:,:,:,sample_index)=gather(hdmi_x);
            hdmi_ally(:,:,:,sample_index)=gather(hdmi_y);
            hdmi_allz(:,:,:,sample_index)=gather(hdmi_z);
        end

        if rk4==1 || (rk4==2 && pc_history_count<3)
            %Step 1
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
        %Heun
        elseif rk4==0
            legacy_empty=1-atomtype_;
            [sxx,syy,szz]=arrayfun(@atomgpu,mmxtmp0+legacy_empty,...
                mmytmp0+legacy_empty,mmztmp0+legacy_empty,...
                scalgpu,alp,tstep,hhx,hhy,hhz);
        %pc4
        else
            legacy_empty=1-atomtype_;
            [sxx,syy,szz]=arrayfun(@atomgpupc4,...
                mmxtmp0+legacy_empty,mmytmp0+legacy_empty,mmztmp0+legacy_empty,...
                pc_x1+legacy_empty,pc_y1+legacy_empty,pc_z1+legacy_empty,...
                pc_x2+legacy_empty,pc_y2+legacy_empty,pc_z2+legacy_empty,...
                pc_x3+legacy_empty,pc_y3+legacy_empty,pc_z3+legacy_empty,...
                scalgpu,alp,tstep,hhx,hhy,hhz);
        end

        sxx=sxx.*atomtype_;syy=syy.*atomtype_;szz=szz.*atomtype_;
        normsn=sqrt(sxx.^2+syy.^2+szz.^2);

        if gather(any((~isfinite(normsn) | normsn==0) & atomtype_~=0,'all'))
            error('Fe2O3:InvalidSpin','Solver produced a nonfinite or zero occupied spin.');
        end
        normsn=normsn+(atomtype_==0);
        sxx=sxx./normsn.*atomtype_;
        syy=syy./normsn.*atomtype_;
        szz=szz./normsn.*atomtype_;

        %pc4 stages
        if rk4==2
            if pc_history_count>=2
                pc_x3=pc_x2;pc_y3=pc_y2;pc_z3=pc_z2;
            end
            if pc_history_count>=1
                pc_x2=pc_x1;pc_y2=pc_y1;pc_z2=pc_z1;
            end
            pc_x1=mmxtmp0;pc_y1=mmytmp0;pc_z1=mmztmp0;
            pc_history_count=min(pc_history_count+1,3);
        end
        mx_current=sxx;my_current=syy;mz_current=szz;
        ct1=ct1+1;

    end
    ct3=ct3+1;
end

clear mx_current my_current mz_current sample_index save_frame
clear mmxtmp0 mmytmp0 mmztmp0 normsn
clear kk1x kk1y kk1z kk2x kk2y kk2z kk3x kk3y kk3z kk4x kk4y kk4z
clear ani_norm ani_denom ani_ux ani_uy
clear pc_x1 pc_y1 pc_z1 pc_x2 pc_y2 pc_z2 pc_x3 pc_y3 pc_z3
clear pc_history_count legacy_empty dipole_refresh dipole_tstep
clear mmxtmpJ1 mmxtmpJ2 mmxtmpJ3 mmxtmpJ4
clear mmytmpJ1 mmytmpJ2 mmytmpJ3 mmytmpJ4
clear mmztmpJ1 mmztmpJ2 mmztmpJ3 mmztmpJ4
clear dmi_sumx dmi_sumy dmi_sumzs

mmx=mmx_.*atomtype_;
mmy=mmy_.*atomtype_;
mmz=mmz_.*atomtype_;
hdmi_allx=hdmi_allx.*atomtype_;
hdmi_ally=hdmi_ally.*atomtype_;
hdmi_allz=hdmi_allz.*atomtype_;
clear mmx_ mmy_ mmz_
