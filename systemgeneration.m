distrib();
mx_init=zeros(natomW,natomL,natomH);%initial magnetization
my_init=zeros(natomW,natomL,natomH);
mz_init=zeros(natomW,natomL,natomH);


  thet_=70/180*pi;
    for ctL=1:natomL
        for ctW=1:natomW
         for ctH=1:natomH
             if mod (ctH,2)==1 
                phi_=30/180*pi;
             else 
                phi_=210/180*pi;
            end
            mx_init(ctW,ctL,ctH)=sin(thet_)*cos(phi_);
            my_init(ctW,ctL,ctH)=sin(thet_)*sin(phi_);
            mz_init(ctW,ctL,ctH)=cos(thet_);
            if atomtype_(ctW, ctL, ctH)==2
            mx_init(ctW,ctL,ctH)=0;
            my_init(ctW,ctL,ctH)=0;
            mz_init(ctW,ctL,ctH)=0;
            end
        end
        end
    end

clear ctL ctW ctH

if loadstartm
    clear mx_init my_init mz_init
    load(startmname);
    if natomW~=natomWcheck || natomL~=natomLcheck
       error('system not consistent') 
    end
    clear natomxcheck natomycheck
    mx_init=mmxstart;
    my_init=mmystart;
    mz_init=mmzstart;
end
