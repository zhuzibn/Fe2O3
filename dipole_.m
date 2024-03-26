 
for i=1:h
   
    switch mod(i-1,3)
        case 0
        [hdipolex_, hdipoley_,hdipolez_]=layer1(muigpu, mu_0,i,atomtype_, atomtype_layer1or,atomtype_layer1gr,atomtype_layer2p,mmxtmp, mmytmp, mmztmp);
        case 1
        [hdipolex_, hdipoley_,hdipolez_]=layer2(muigpu, mu_0,i,atomtype_,atomtype_layer1or,atomtype_layer1gr,atomtype_layer2p, mmxtmp, mmytmp, mmztmp);
        case 2
        [hdipolex_, hdipoley_,hdipolez_]=layer3(muigpu, mu_0,i,atomtype_, atomtype_layer1or,atomtype_layer1gr,atomtype_layer2p,mmxtmp, mmytmp, mmztmp);
    
        otherwise
            error('unexpected index');
   
    end
    hdipo_x(:,:,i)= hdipolex_(:,:,i); hdipo_y(:,:,i)= hdipoley_(:,:,i); hdipo_z(:,:,i)= hdipolez_(:,:,i);
end
clear hdipolex_ hdipoley_ hdipolez_

function  [hdipolex_, hdipoley_,hdipolez_]=layer1(muigpu, mu_0,i,atomtype_,atomtype_layer1or,atomtype_layer1gr,atomtype_layer2p,mmxtmp, mmytmp, mmztmp)
        
  
          rijx_=0.5067*1e-9*1;rijy_=0;
          dot_sr=muigpu.*((circshift(atomtype_layer1or.*mmxtmp(:,:,i),[0,-2])).*rijx_+(circshift(atomtype_layer1or.*mmytmp(:,:,i),[0,-2])).*rijy_);
          dist_=0.5067*1e-9;
          hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[0,-2]))./dist_.^3);%[T]
       
          hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmytmp(:,:,i),[0,-2]))./dist_.^3);
        
          hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1or.*mmztmp(:,:,i),[0,-2]))./dist_.^3);
          
         
          rijx_=0.5067*1e-9*-1;rijy_=0;
          dot_sr=muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[0,2]).*rijx_+circshift(atomtype_layer1gr.*mmytmp(:,:,i),[0,2]).*rijy_);
      
          hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[0,2]))./dist_.^3);%[T]
          
          hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[0,2]))./dist_.^3);
        
          hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1gr.*mmztmp(:,:,i),[0,2]))./dist_.^3);
        
          %%green-orange (A-1,i) 
          rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
          dot_sr=muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[-1,0]).*rijx_+circshift(atomtype_layer1or.*mmytmp(:,:,i),[-1,0]).*rijy_);
          hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[-1,0]))./dist_.^3);%[T]
         
          hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmytmp(:,:,i),[-1,0]))./dist_.^3);
         
          hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1or.*mmztmp(:,:,i),[-1,0]))./dist_.^3);
         
          rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
          dot_sr=muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[-1,0]).*rijx_+circshift(atomtype_layer1gr.*mmytmp(:,:,i),[-1,0]).*rijy_);
          hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[-1,0]))./dist_.^3);%[T]
         
          hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[-1,0]))./dist_.^3);
        
          hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1gr.*mmztmp(:,:,i),[-1,0]))./dist_.^3);

         
          rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
          dot_sr=muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[1,0]).*rijx_+circshift(atomtype_layer1or.*mmytmp(:,:,i),[1,0]).*rijy_);
          hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[1,0]))./dist_.^3);%[T]
         
          hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmytmp(:,:,i),[1,0]))./dist_.^3);
        
          hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1or.*mmztmp(:,:,i),[1,0]))./dist_.^3);
          
          rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
          dot_sr=muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[1,0]).*rijx_+circshift(atomtype_layer1gr.*mmytmp(:,:,i),[1,0]).*rijy_);
          hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[1,0]))./dist_.^3);%[T]
         
          hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[1,0]))./dist_.^3);
        
          hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1gr.*mmztmp(:,:,i),[1,0]))./dist_.^3);

          hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
          hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
          hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);
          
end

function  [hdipolex_, hdipoley_,hdipolez_]=layer2(muigpu, mu_0,i,atomtype_,atomtype_layer1or,atomtype_layer1gr,atomtype_layer2p,mmxtmp, mmytmp, mmztmp)

      
          
          rijx_=0.5067*1e-9;rijy_=0;
          dot_sr=muigpu.*(circshift(atomtype_layer2p.*mmxtmp(:,:,i),[0,-1]).*rijx_+circshift(atomtype_layer2p.*mmytmp(:,:,i),[0,-1]).*rijy_);
          dist_=0.5067*1e-9;
          hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmxtmp(:,:,i),[0,-1]))./dist_.^3);%[T]
       
          hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmytmp(:,:,i),[0,-1]))./dist_.^3);
        
          hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer2p.*mmztmp(:,:,i),[0,-1]))./dist_.^3);
          
      
          rijx_=-0.5067*1e-9;rijy_=0;
          dot_sr=muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[0,1]).*rijx_+circshift(atomtype_layer1or.*mmytmp(:,:,i),[0,1]).*rijy_);
          hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[0,1]))./dist_.^3);%[T]
          
          hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmytmp(:,:,i),[0,1]))./dist_.^3);
        
          hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1or.*mmztmp(:,:,i),[0,1]))./dist_.^3);
        
          
          rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
          dot_sr=muigpu.*(circshift(atomtype_layer2p.*mmxtmp(:,:,i),[-1,1]).*rijx_+circshift(atomtype_layer2p.*mmytmp(:,:,i),[-1,1]).*rijy_);
          hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmxtmp(:,:,i),[-1,1]))./dist_.^3);%[T]
         
          hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmytmp(:,:,i),[-1,1]))./dist_.^3);
         
          hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer2p.*mmztmp(:,:,i),[-1,1]))./dist_.^3);
    

          rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
          dot_sr=muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[-1,-1]).*rijx_+circshift(atomtype_layer1or.*mmytmp(:,:,i),[-1,-1]).*rijy_);
          hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[-1,-1]))./dist_.^3);%[T]
         
          hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmytmp(:,:,i),[-1,-1]))./dist_.^3);
        
          hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1or.*mmztmp(:,:,i),[-1,-1]))./dist_.^3);

          rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
          dot_sr=muigpu.*((circshift(atomtype_layer2p.*mmxtmp(:,:,i),[1,1])).*rijx_+(circshift(atomtype_layer2p.*mmytmp(:,:,i),[-1,1])).*rijy_);
          hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmxtmp(:,:,i),[1,1]))./dist_.^3);%[T]
         
          hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmytmp(:,:,i),[1,1]))./dist_.^3);
        
          hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer2p.*mmztmp(:,:,i),[1,1]))./dist_.^3);
          
          rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
          dot_sr=muigpu.*((circshift(atomtype_layer1or.*mmxtmp(:,:,i),[1,-1])).*rijx_+(circshift(atomtype_layer1or.*mmytmp(:,:,i),[1,-1])).*rijy_);
          hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmxtmp(:,:,i),[1,-1]))./dist_.^3);%[T]
         
          hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1or.*mmytmp(:,:,i),[1,-1]))./dist_.^3);
        
          hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1or.*mmztmp(:,:,i),[1,-1]))./dist_.^3);

         
          hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
          hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
          hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);
        
       
end
function  [hdipolex_, hdipoley_,hdipolez_]=layer3(muigpu, mu_0,i,atomtype_,atomtype_layer1or,atomtype_layer1gr,atomtype_layer2p,mmxtmp, mmytmp, mmztmp)
 
          
          rijx_=0.5067*1e-9;rijy_=0;
          dot_sr=muigpu.*((circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[0,-1])).*rijx_+(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[0,-1])).*rijy_);
          dist_=0.5067*1e-9;
          hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[0,-1]))./dist_.^3);%[T]
       
          hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[0,-1]))./dist_.^3);
        
          hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1gr.*mmztmp(:,:,i),[0,-1]))./dist_.^3);
          
          rijx_=-0.5067*1e-9;rijy_=0;
          dot_sr=muigpu.*((circshift(atomtype_layer2p.*mmxtmp(:,:,i),[0,1])).*rijx_+(circshift(atomtype_layer2p.*mmytmp(:,:,i),[0,1])).*rijy_);
          hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmxtmp(:,:,i),[0,1]))./dist_.^3);%[T]
          
          hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmytmp(:,:,i),[0,1]))./dist_.^3);
        
          hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer2p.*mmztmp(:,:,i),[0,1]))./dist_.^3);
        
          rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
          dot_sr=muigpu.*((circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[-1,1])).*rijx_+(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[-1,1])).*rijy_);
          hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[-1,1]))./dist_.^3);%[T]
         
          hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[-1,1]))./dist_.^3);
         
          hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1gr.*mmztmp(:,:,i),[-1,1]))./dist_.^3);
         
          rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
          dot_sr=muigpu.*((circshift(atomtype_layer2p.*mmxtmp(:,:,i),[-1,-1])).*rijx_+(circshift(atomtype_layer2p.*mmytmp(:,:,i),[-1,-1])).*rijy_);
          hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmxtmp(:,:,i),[-1,-1]))./dist_.^3);%[T]
         
          hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmytmp(:,:,i),[-1,-1]))./dist_.^3);
        
          hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer2p.*mmztmp(:,:,i),[-1,-1]))./dist_.^3);

          rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
          dot_sr=muigpu.*((circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[1,1])).*rijx_+(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[1,1])).*rijy_);
          hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmxtmp(:,:,i),[1,1]))./dist_.^3);%[T]
         
          hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer1gr.*mmytmp(:,:,i),[1,1]))./dist_.^3);
        
          hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer1gr.*mmztmp(:,:,i),[1,1]))./dist_.^3);
          
          rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
          dot_sr=muigpu.*((circshift(atomtype_layer2p.*mmxtmp(:,:,i),[1,-1])).*rijx_+(circshift(atomtype_layer2p.*mmytmp(:,:,i),[1,-1])).*rijy_);
          hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmxtmp(:,:,i),[1,-1]))./dist_.^3);%[T]
         
          hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(circshift(atomtype_layer2p.*mmytmp(:,:,i),[1,-1]))./dist_.^3);
        
          hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(circshift(atomtype_layer2p.*mmztmp(:,:,i),[1,-1]))./dist_.^3);

         
          hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
          hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
          hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);
      
end

