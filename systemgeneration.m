distrib();
mx_init=zeros(natomW,natomL,natomH);%initial magnetization
my_init=zeros(natomW,natomL,natomH);
mz_init=zeros(natomW,natomL,natomH);

if dwcalc
    phi_=0;
    for ctL=1:natomL
        for ctW=1:natomW
            for ctH=1:natomH
                if ctL<round(natomL/2)
                    if atomtype_(ctW,ctL,ctH)==1%TM
                        thet_=5/180*pi;
                    else
                        thet_=0;
                    end
                    mx_init(ctW,ctL)=sin(thet_)*cos(phi_);
                    my_init(ctW,ctL)=sin(thet_)*sin(phi_);
                    mz_init(ctW,ctL)=cos(thet_);
                else
                    if atomtype_(ctW,ctL)==1%TM
                        thet_=(5+180)/180*pi;
                    else
                        thet_=0;
                    end
                    mx_init(ctW,ctL)=sin(thet_)*cos(phi_);
                    my_init(ctW,ctL)=sin(thet_)*sin(phi_);
                    mz_init(ctW,ctL)=cos(thet_);
                end
            end
        end
    end
else
    thet_=75/180*pi;
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
end


clear ctL ctW ctH


