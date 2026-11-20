distrib();
mx_init=zeros(natomW,natomL,natomH);%initial magnetization
my_init=zeros(natomW,natomL,natomH);
mz_init=zeros(natomW,natomL,natomH);
if load_fixed_atom_distrib
    restart_geometry=struct('atomtype_',atomtype_,...
        'atomtype_layer1gr',atomtype_layer1gr,'atomtype_layer1or',atomtype_layer1or,...
        'atomtype_layer2blue',atomtype_layer2blue,'atomtype_layer2p',atomtype_layer2p,...
        'atomtype_layer3red',atomtype_layer3red,'atomtype_layer3black',atomtype_layer3black,...
        'array_order','W,L,H,time');
    [mx_init,my_init,mz_init,restart_info]=fe2o3_load_restart(restart_file,restart_geometry);
else
    thet_=90/180*pi;
    for ctL=1:natomL
        for ctW=1:natomW
            for ctH=1:natomH
                if mod (ctH,2)==1
                    phi_=0/180*pi;
                else
                    phi_=180/180*pi;
                end
                mx_init(ctW,ctL,ctH)=sin(thet_)*cos(phi_);
                my_init(ctW,ctL,ctH)=sin(thet_)*sin(phi_);
                mz_init(ctW,ctL,ctH)=cos(thet_);
                if atomtype_(ctW, ctL, ctH)==0
                    mx_init(ctW,ctL,ctH)=0;
                    my_init(ctW,ctL,ctH)=0;
                    mz_init(ctW,ctL,ctH)=0;
                end
            end
        end
    end
end


clear ctL ctW ctH


