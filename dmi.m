for i=1:h
    if i==1
        [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layerfirst(natomW, natomL,atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,mmxtmp, mmytmp, mmztmp);
    elseif i==h
         [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layerend(h,natomW, natomL,atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,mmxtmp, mmytmp, mmztmp);
    else      
    switch mod(i-1,3)
        case 0
        [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layer1(atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,i, mmxtmp, mmytmp, mmztmp);
        case 1
        [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layer2(atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,i, mmxtmp, mmytmp, mmztmp);
        case 2
        [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layer3(atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,i, mmxtmp, mmytmp, mmztmp);
    
        otherwise
            error('unexpected index');
    end
    end
        mmxtmpdmi_pre(:,:,i)=mmxtmpd_pre(:,:,i);mmytmpdmi_pre(:,:,i)=mmytmpd_pre(:,:,i);mmztmpdmi_pre(:,:,i)=mmztmpd_pre(:,:,i);
        mmxtmpdmi_nex(:,:,i)=mmxtmpd_nex(:,:,i);mmytmpdmi_nex(:,:,i)=mmytmpd_nex(:,:,i);mmztmpdmi_nex(:,:,i)=mmztmpd_nex(:,:,i);
end

function [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layerfirst(natomW, natomL,atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,mmxtmp, mmytmp, mmztmp)
        
       mmxtmpd_pre(:,:,1)=circshift(atomtype_layer2p.*mmxtmp(:,:,2),[0,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,2),[-1,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,2),[1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,2),[0,-2])+circshift(atomtype_layer1or.*mmxtmp(:,:,2),[-1,0])+circshift(atomtype_layer1or.*mmxtmp(:,:,2),[1,0])...
           +circshift(atomtype_layer2p.*mmxtmp(:,:,2),[0,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,2),[1,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,2),[-1,1]);   
       mmytmpd_pre(:,:,1)=circshift(atomtype_layer2p.*mmytmp(:,:,2),[0,1])+circshift(atomtype_layer2p.*mmytmp(:,:,2),[-1,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,2),[1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,2),[0,-2])+circshift(atomtype_layer1or.*mmytmp(:,:,2),[-1,0])+circshift(atomtype_layer1or.*mmytmp(:,:,2),[1,0])...
           +circshift(atomtype_layer2p.*mmytmp(:,:,2),[0,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,2),[1,1])+circshift(atomtype_layer2p.*mmytmp(:,:,2),[-1,1]);
       mmztmpd_pre(:,:,1)=circshift(atomtype_layer2p.*mmztmp(:,:,2),[0,1])+circshift(atomtype_layer2p.*mmztmp(:,:,2),[-1,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,2),[1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,2),[0,-2])+circshift(atomtype_layer1or.*mmztmp(:,:,2),[-1,0])+circshift(atomtype_layer1or.*mmztmp(:,:,2),[1,0])...
           +circshift(atomtype_layer2p.*mmztmp(:,:,2),[0,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,2),[1,1])+circshift(atomtype_layer2p.*mmztmp(:,:,2),[-1,1]);     
       mmxtmpd_nex(:,:,1)=zeros(natomW, natomL,'gpuArray');
       mmytmpd_nex(:,:,1)=zeros(natomW, natomL,'gpuArray');
       mmztmpd_nex(:,:,1)=zeros(natomW, natomL,'gpuArray');

end
function [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layerend(h,natomW, natomL,atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,mmxtmp, mmytmp, mmztmp)
      
       mmxtmpd_nex(:,:,h)=circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[0,1])+circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[-1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[1,0])+circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[0,-2])+circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[-1,0])...
           +circshift(atomtype_layer2p.*mmxtmp(:,:,h-1),[0,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,h-1),[-1,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,h-1),[1,-1]);   
       mmytmpd_nex(:,:,h)=circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[0,1])+circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[-1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[1,0])+circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[0,-2])+circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[-1,0])...
           +circshift(atomtype_layer2p.*mmytmp(:,:,h-1),[0,1])+circshift(atomtype_layer2p.*mmytmp(:,:,h-1),[-1,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,h-1),[1,-1]);   
       mmztmpd_nex(:,:,h)=circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[0,1])+circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[-1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[1,0])+circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[0,-2])+circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[-1,0])...
           +circshift(atomtype_layer2p.*mmztmp(:,:,h-1),[0,1])+circshift(atomtype_layer2p.*mmztmp(:,:,h-1),[-1,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,h-1),[1,-1]);   
       mmxtmpd_pre(:,:,h)=zeros(natomW, natomL,'gpuArray');
       mmytmpd_pre(:,:,h)=zeros(natomW, natomL,'gpuArray');
       mmztmpd_pre(:,:,h)=zeros(natomW, natomL,'gpuArray');

end
function  [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layer1(atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,i,mmxtmp, mmytmp, mmztmp)
       mmxtmpd_pre(:,:,i)=circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[0,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[-1,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[0,-2])+circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[1,0])...
           +circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[0,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[1,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[-1,1]);   
       mmytmpd_pre(:,:,i)=circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[0,1])+circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[-1,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[0,-2])+circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[1,0])...
           +circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[0,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[1,1])+circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[-1,1]);
       mmztmpd_pre(:,:,i)=circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[0,1])+circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[-1,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[0,-2])+circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[1,0])...
           +circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[0,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[1,1])+circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[-1,1]);     
      
       mmxtmpd_nex(:,:,i)=circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[0,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[0,2])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[-1,0])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[1,0])...
           +circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[0,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[1,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[-1,1]);   
       mmytmpd_nex(:,:,i)=circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[0,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[0,2])+circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[-1,0])+circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[1,0])...
           +circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[0,1])+circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[1,1])+circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[-1,1]);   
       mmztmpd_nex(:,:,i)=circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[0,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[0,2])+circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[-1,0])+circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[1,0])...
           +circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[0,1])+circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[1,1])+circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[-1,1]);   
end

function  [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layer2(atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,i,mmxtmp, mmytmp, mmztmp)
        
       mmxtmpd_nex(:,:,i)=circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[0,1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[1,1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[-1,1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[0,-1])...
           +circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[0,2])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[-1,0])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i-1),[1,0]);   
       mmytmpd_nex(:,:,i)=circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[0,1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[1,1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[-1,1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[0,-1])...
           +circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[0,2])+circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[-1,0])+circshift(atomtype_layer1gr.*mmytmp(:,:,i-1),[1,0]);   
       mmztmpd_nex(:,:,i)=circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[0,1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[1,1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[-1,1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[0,-1])...
           +circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[0,2])+circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[-1,0])+circshift(atomtype_layer1gr.*mmztmp(:,:,i-1),[1,0]);   
       
       mmxtmpd_pre(:,:,i)=circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[0,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[-1,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i+1),[1,1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[0,-1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[-1,1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[1,1])...
           +circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[0,2])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[1,0]);   
       mmytmpd_pre(:,:,i)=circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[0,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[-1,1])+circshift(atomtype_layer2p.*mmytmp(:,:,i+1),[1,1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[0,-1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[-1,1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[1,1])...
           +circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[0,2])+circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[1,0]);   
       mmztmpd_pre(:,:,i)=circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[0,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[-1,1])+circshift(atomtype_layer2p.*mmztmp(:,:,i+1),[1,1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[0,-1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[-1,1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[1,1])...
           +circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[0,2])+circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[1,0]);   
 
end
function  [mmxtmpd_pre, mmytmpd_pre,mmztmpd_pre, mmxtmpd_nex, mmytmpd_nex,mmztmpd_nex]=layer3(atomtype_layer1gr,atomtype_layer1or,atomtype_layer2p,i,mmxtmp, mmytmp, mmztmp)
      
       mmxtmpd_nex(:,:,i)=circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[0,1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[1,0])+circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[0,-2])+circshift(atomtype_layer1or.*mmxtmp(:,:,i-1),[-1,0])...
           +circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[0,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,i-1),[1,-1]);   
       mmytmpd_nex(:,:,i)=circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[0,1])+circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[1,0])+circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[0,-2])+circshift(atomtype_layer1or.*mmytmp(:,:,i-1),[-1,0])...
           +circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[0,1])+circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,i-1),[1,-1]);   
       mmztmpd_nex(:,:,i)=circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[0,1])+circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[1,0])+circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[0,-2])+circshift(atomtype_layer1or.*mmztmp(:,:,i-1),[-1,0])...
           +circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[0,1])+circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[-1,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,i-1),[1,-1]);  
       
       mmxtmpd_pre(:,:,i)=circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[0,1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[-1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[1,-1])+circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[0,-2])+circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1or.*mmxtmp(:,:,i+1),[1,0])...
           +circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[0,-1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[-1,1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,i+1),[1,1]);   
       mmytmpd_pre(:,:,i)=circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[0,1])+circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[-1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[1,-1])+circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[0,-2])+circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1or.*mmytmp(:,:,i+1),[1,0])...
           +circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[0,-1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[-1,1])+circshift(atomtype_layer1gr.*mmytmp(:,:,i+1),[1,1]);   
       mmztmpd_pre(:,:,i)=circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[0,1])+circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[-1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[1,-1])+circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[0,-2])+circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[-1,0])+circshift(atomtype_layer1or.*mmztmp(:,:,i+1),[1,0])...
           +circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[0,-1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[-1,1])+circshift(atomtype_layer1gr.*mmztmp(:,:,i+1),[1,1]);   
 end
              
