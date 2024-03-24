        
for i=1:h
       

     if i==1

       mmxtmpJ1(:,:,1)=atomtype_layer1or.*mmxtmp(:,:,2);
       mmytmpJ1(:,:,1)=atomtype_layer1or.*mmytmp(:,:,2);
       mmztmpJ1(:,:,1)=atomtype_layer1or.*mmztmp(:,:,2);

       shift_mat=circshift(atomtype_layer1or.*mmxtmp(:,:,1),[0,-2,0])+circshift(atomtype_layer1gr.*mmxtmp(:,:,1),[0,2,0])+circshift(atomtype_layer1or.*mmxtmp(:,:,1),[1,0,0])+circshift(atomtype_layer1gr.*mmxtmp(:,:,1),[1,0,0])...
       +circshift(atomtype_layer1or.*mmxtmp(:,:,1),[-1,0,0])+circshift(atomtype_layer1gr.*mmxtmp(:,:,1),[-1,0,0]);%length direction -2 related to +2(right)
       mmxtmpJ2(:,:,1)=shift_mat;
       shift_mat=circshift(atomtype_layer1or.*mmytmp(:,:,1),[0,-2,0])+circshift(atomtype_layer1gr.*mmytmp(:,:,1),[0,2,0])+circshift(atomtype_layer1or.*mmytmp(:,:,1),[1,0,0])+circshift(atomtype_layer1gr.*mmytmp(:,:,1),[1,0,0])...
       +circshift(atomtype_layer1or.*mmytmp(:,:,1),[-1,0,0])+circshift(atomtype_layer1gr.*mmytmp(:,:,1),[-1,0,0]);%length direction -2 related to +2(right)
       mmytmpJ2(:,:,1)=shift_mat;
       shift_mat=circshift(atomtype_layer1or.*mmztmp(:,:,1),[0,-2,0])+circshift(atomtype_layer1gr.*mmztmp(:,:,1),[0,2,0])+circshift(atomtype_layer1or.*mmztmp(:,:,1),[1,0,0])+circshift(atomtype_layer1gr.*mmztmp(:,:,1),[1,0,0])...
       +circshift(atomtype_layer1or.*mmztmp(:,:,1),[-1,0,0])+circshift(atomtype_layer1gr.*mmztmp(:,:,1),[-1,0,0]);%length direction -2 related to +2(right)
       mmztmpJ2(:,:,1)=shift_mat;

       shift_mat(:,:,1)=circshift(atomtype_layer2p.*mmxtmp(:,:,2),[0,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,2),[-1,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,2),[1,-1]);
       mmxtmpJ3(:,:,1)=shift_mat;
       shift_mat(:,:,1)=circshift(atomtype_layer2p.*mmytmp(:,:,2),[0,1])+circshift(atomtype_layer2p.*mmytmp(:,:,2),[-1,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,2),[1,-1]);
       mmytmpJ3(:,:,1)=shift_mat;
       shift_mat(:,:,1)=circshift(atomtype_layer2p.*mmztmp(:,:,2),[0,1])+circshift(atomtype_layer2p.*mmztmp(:,:,2),[-1,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,2),[1,-1]);
       mmztmpJ3(:,:,1)=shift_mat;

       shift_mat=circshift(atomtype_layer1or.*mmxtmp(:,:,2),[0,-2])+circshift(atomtype_layer1or.*mmxtmp(:,:,2),[-1,0])+circshift(atomtype_layer1or.*mmxtmp(:,:,2),[1,0])...
            +circshift(atomtype_layer2p.*mmxtmp(:,:,2),[0,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,2),[1,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,2),[-1,1]);
       mmxtmpJ4(:,:,1)=shift_mat;
       shift_mat=circshift(atomtype_layer1or.*mmytmp(:,:,2),[0,-2])+circshift(atomtype_layer1or.*mmytmp(:,:,2),[-1,0])+circshift(atomtype_layer1or.*mmytmp(:,:,2),[1,0])...
            +circshift(atomtype_layer2p.*mmytmp(:,:,2),[0,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,2),[1,1])+circshift(atomtype_layer2p.*mmytmp(:,:,2),[-1,1]);
       mmytmpJ4(:,:,1)=shift_mat;
       shift_mat=circshift(atomtype_layer1or.*mmztmp(:,:,2),[0,-2])+circshift(atomtype_layer1or.*mmztmp(:,:,2),[-1,0])+circshift(atomtype_layer1or.*mmztmp(:,:,2),[1,0])...
            +circshift(atomtype_layer2p.*mmztmp(:,:,2),[0,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,2),[1,1])+circshift(atomtype_layer2p.*mmztmp(:,:,2),[-1,1]);
       mmztmpJ4(:,:,1)=shift_mat;
     elseif i==h
         mmxtmpJ1(:,:,h)=atomtype_layer2p.*mmxtmp(:,:,h-1);
         mmytmpJ1(:,:,h)=atomtype_layer2p.*mmytmp(:,:,h-1);
         mmztmpJ1(:,:,h)=atomtype_layer2p.*mmztmp(:,:,h-1);

         shift_mat=circshift(atomtype_layer1gr.*mmxtmp(:,:,h),[0,-1])+circshift(atomtype_layer2p.*mmxtmp(:,:,h),[0,1])+circshift(atomtype_layer1gr.*mmxtmp(:,:,h),[1,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,h),[1,-1])...%length direction 
         +circshift(atomtype_layer1gr.*mmxtmp(:,:,h),[-1,1])+circshift(atomtype_layer2p.*mmxtmp(:,:,h),[-1,-1]);
         mmxtmpJ2(:,:,h)=shift_mat;
         shift_mat=circshift(atomtype_layer1gr.*mmytmp(:,:,h),[0,-1])+circshift(atomtype_layer2p.*mmytmp(:,:,h),[0,1])+circshift(atomtype_layer1gr.*mmytmp(:,:,h),[1,1])+circshift(atomtype_layer2p.*mmytmp(:,:,h),[1,-1])...%length direction 
         +circshift(atomtype_layer1gr.*mmytmp(:,:,h),[-1,1])+circshift(atomtype_layer2p.*mmytmp(:,:,h),[-1,-1]);
         mmytmpJ2(:,:,h)=shift_mat;
         shift_mat=circshift(atomtype_layer1gr.*mmztmp(:,:,h),[0,-1])+circshift(atomtype_layer2p.*mmztmp(:,:,h),[0,1])+circshift(atomtype_layer1gr.*mmztmp(:,:,h),[1,1])+circshift(atomtype_layer2p.*mmztmp(:,:,h),[1,-1])...%length direction 
         +circshift(atomtype_layer1gr.*mmztmp(:,:,h),[-1,1])+circshift(atomtype_layer2p.*mmztmp(:,:,h),[-1,-1]);
         mmztmpJ2(:,:,h)=shift_mat;

         shift_mat=circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[0,-2])+circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[-1,0])+circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[1,0]);
         mmxtmpJ3(:,:,h)=shift_mat;
         shift_mat=circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[0,-2])+circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[-1,0])+circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[1,0]);
         mmytmpJ3(:,:,h)=shift_mat;
         shift_mat=circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[0,-2])+circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[-1,0])+circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[1,0]);
         mmztmpJ3(:,:,h)=shift_mat;

         shift_mat=circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[0,1]) +circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[-1,-1]) +circshift(atomtype_layer1or.*mmxtmp(:,:,h-1),[1,-1])...
             +circshift(atomtype_layer2p.*mmxtmp(:,:,h-1),[0,1]) +circshift(atomtype_layer2p.*mmxtmp(:,:,h-1),[-1,-1]) +circshift(atomtype_layer2p.*mmxtmp(:,:,h-1),[1,-1]);
         mmxtmpJ4(:,:,h)=shift_mat;
         shift_mat=circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[0,1]) +circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[-1,-1]) +circshift(atomtype_layer1or.*mmytmp(:,:,h-1),[1,-1])...
             +circshift(atomtype_layer2p.*mmytmp(:,:,h-1),[0,1]) +circshift(atomtype_layer2p.*mmytmp(:,:,h-1),[-1,-1]) +circshift(atomtype_layer2p.*mmytmp(:,:,h-1),[1,-1]);
         mmytmpJ4(:,:,h)=shift_mat;
         shift_mat=circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[0,1]) +circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[-1,-1]) +circshift(atomtype_layer1or.*mmztmp(:,:,h-1),[1,-1])...
             +circshift(atomtype_layer2p.*mmztmp(:,:,h-1),[0,1]) +circshift(atomtype_layer2p.*mmztmp(:,:,h-1),[-1,-1]) +circshift(atomtype_layer2p.*mmztmp(:,:,h-1),[1,-1]);
         mmztmpJ4(:,:,h)=shift_mat;


     else
  tmpgr_add_x= atomtype_layer1gr.*mmxtmp(:,:,i+1);
         tmpor_add_x= atomtype_layer1or.*mmxtmp(:,:,i+1);
         tmpp_add_x= atomtype_layer2p.*mmxtmp(:,:,i+1);
         tmpgr_min_x= atomtype_layer1gr.*mmxtmp(:,:,i-1);
         tmpor_min_x= atomtype_layer1or.*mmxtmp(:,:,i-1);
         tmpp_min_x= atomtype_layer2p.*mmxtmp(:,:,i-1);

         tmpgr_add_y= atomtype_layer1gr.*mmytmp(:,:,i+1);
         tmpor_add_y= atomtype_layer1or.*mmytmp(:,:,i+1);
         tmpp_add_y= atomtype_layer2p.*mmytmp(:,:,i+1);
         tmpgr_min_y= atomtype_layer1gr.*mmytmp(:,:,i-1);
         tmpor_min_y= atomtype_layer1or.*mmytmp(:,:,i-1);
         tmpp_min_y= atomtype_layer2p.*mmytmp(:,:,i-1);

         tmpgr_add_z= atomtype_layer1gr.*mmztmp(:,:,i+1);
         tmpor_add_z= atomtype_layer1or.*mmztmp(:,:,i+1);
         tmpp_add_z= atomtype_layer2p.*mmztmp(:,:,i+1);
         tmpgr_min_z= atomtype_layer1gr.*mmztmp(:,:,i-1);
         tmpor_min_z= atomtype_layer1or.*mmztmp(:,:,i-1);
         tmpp_min_z= atomtype_layer2p.*mmztmp(:,:,i-1);

         tmpgr_x= atomtype_layer1gr.*mmxtmp(:,:,i);
         tmpor_x= atomtype_layer1or.*mmxtmp(:,:,i);
         tmpp_x= atomtype_layer2p.*mmxtmp(:,:,i);

         tmpgr_y= atomtype_layer1gr.*mmytmp(:,:,i);
         tmpor_y= atomtype_layer1or.*mmytmp(:,:,i);
         tmpp_y= atomtype_layer2p.*mmytmp(:,:,i);

         tmpgr_z= atomtype_layer1gr.*mmztmp(:,:,i);
         tmpor_z= atomtype_layer1or.*mmztmp(:,:,i);
         tmpp_z= atomtype_layer2p.*mmztmp(:,:,i);
    switch mod(i-1,3)
        case 0
         
        
         mmxtmpJ1(:,:,i)=tmpgr_add_x+tmpor_min_x;
         mmytmpJ1(:,:,i)=tmpgr_add_y+tmpor_min_y;
         mmztmpJ1(:,:,i)=tmpgr_add_z+tmpor_min_z;

         shift_mat=circshift(tmpor_x,[0,-2,0])+circshift(tmpgr_x,[0,2,0])+circshift((tmpor_x+tmpgr_x),[1,0,0])...%length direction 
         +circshift((tmpor_x+tmpgr_x),[-1,0,0]);%length direction 
         mmxtmpJ2(:,:,i)=shift_mat;
         shift_mat=circshift(tmpor_y,[0,-2,0])+circshift(tmpgr_y,[0,2,0])+circshift((tmpor_y+tmpgr_y),[1,0,0])...%length direction 
         +circshift((tmpor_y+tmpgr_y),[-1,0,0]);%length direction 
         mmytmpJ2(:,:,i)=shift_mat;
         shift_mat=circshift(tmpor_z,[0,-2,0])+circshift(tmpgr_z,[0,2,0])+circshift((tmpor_z+tmpgr_z),[1,0,0])...%length direction 
         +circshift((tmpor_z+tmpgr_z),[-1,0,0]);%length direction 
         mmztmpJ2(:,:,i)=shift_mat;

         shift_mat=circshift(tmpp_min_x,[0,-1])+circshift(tmpp_add_x,[0,1])+circshift(tmpp_min_x,[-1,1])+circshift(tmpp_add_x,[-1,-1])...
             +circshift(tmpp_min_x,[1,1])+circshift(tmpp_add_x,[1,-1]);
         mmxtmpJ3(:,:,i)=shift_mat;
         shift_mat=circshift(tmpp_min_y,[0,-1])+circshift(tmpp_add_y,[0,1])+circshift(tmpp_min_y,[-1,1])+circshift(tmpp_add_y,[-1,-1])...
          +circshift(tmpp_min_y,[1,1])+circshift(tmpp_add_y,[1,-1]);
         mmytmpJ3(:,:,i)=shift_mat;
         shift_mat=circshift(tmpp_min_z,[0,-1])+circshift(tmpp_add_z,[0,1])+circshift(tmpp_min_z,[-1,1])+circshift(tmpp_add_z,[-1,-1])...
             + circshift(tmpp_min_z,[1,1])+circshift(tmpp_add_z,[1,-1]);
         mmztmpJ3(:,:,i)=shift_mat;
        
         shift_mat=circshift(tmpor_add_x,[0,-2])+circshift(tmpor_add_x,[-1,0])+circshift(tmpor_add_x,[1,0])...
            +circshift(tmpp_add_x,[0,-1])+circshift(tmpp_add_x,[1,1])+circshift(tmpp_add_x,[-1,1])...
            +circshift(tmpp_min_x,[0,1])+circshift(tmpp_min_x,[-1,-1])+circshift(tmpp_min_x,[1,-1])...
            +circshift(tmpgr_min_x,[0,2])+circshift(tmpgr_min_x,[1,0])+circshift(tmpgr_min_x,[-1,0]);
         mmxtmpJ4(:,:,i)=shift_mat;
         shift_mat=circshift(tmpor_add_y,[0,-2])+circshift(tmpor_add_y,[-1,0])+circshift(tmpor_add_y,[1,0])...
            +circshift(tmpp_add_y,[0,-1])+circshift(tmpp_add_y,[1,1])+circshift(tmpp_add_y,[-1,1])...
            +circshift(tmpp_min_y,[0,1])+circshift(tmpp_min_y,[-1,-1])+circshift(tmpp_min_y,[1,-1])...
            +circshift(tmpgr_min_y,[0,2])+circshift(tmpgr_min_y,[1,0])+circshift(tmpgr_min_y,[-1,0]);
         mmytmpJ4(:,:,i)=shift_mat;
         shift_mat=circshift(tmpor_add_z,[0,-2])+circshift(tmpor_add_z,[-1,0])+circshift(tmpor_add_z,[1,0])...
            +circshift(tmpp_add_z,[0,-1])+circshift(tmpp_add_z,[1,1])+circshift(tmpp_add_z,[-1,1])...
            +circshift(tmpp_min_z,[0,1])+circshift(tmpp_min_z,[-1,-1])+circshift(tmpp_min_z,[1,-1])...
            +circshift(tmpgr_min_z,[0,2])+circshift(tmpgr_min_z,[1,0])+circshift(tmpgr_min_z,[-1,0]);
         mmztmpJ4(:,:,i)=shift_mat;

        case 1
         mmxtmpJ1(:,:,i)=tmpp_add_x+tmpor_min_x;
         mmytmpJ1(:,:,i)=tmpp_add_y+tmpor_min_y;
         mmztmpJ1(:,:,i)=tmpp_add_z+tmpor_min_z;
 
         shift_mat=circshift(tmpp_x,[0,-1,0])+circshift(tmpor_x,[0,1,0])+circshift(tmpp_x,[1,1,0])+circshift(tmpor_x,[1,-1,0])...%length direction -2 related to +2(right)
         +circshift(tmpp_x,[-1,1,0])+circshift(tmpor_x,[-1,-1,0]);
         mmxtmpJ2(:,:,i)=shift_mat;
         shift_mat=circshift(tmpp_y,[0,-1,0])+circshift(tmpor_y,[0,1,0])+circshift(tmpp_y,[1,1,0])+circshift(tmpor_y,[1,-1,0])...%length direction -2 related to +2(right)
         +circshift(tmpp_y,[-1,1,0])+circshift(tmpor_y,[-1,-1,0]);
         mmytmpJ2(:,:,i)=shift_mat;
         shift_mat=circshift(tmpp_z,[0,-1,0])+circshift(tmpor_z,[0,1,0])+circshift(tmpp_z,[1,1,0])+circshift(tmpor_z,[1,-1,0])...%length direction -2 related to +2(right)
         +circshift(tmpp_z,[-1,1,0])+circshift(tmpor_z,[-1,-1,0]);
         mmztmpJ2(:,:,i)=shift_mat;

         shift_mat=circshift(tmpgr_min_x,[0,-1])+circshift(tmpgr_add_x,[0,2])+circshift(tmpgr_min_x,[-1,1])+circshift(tmpgr_add_x,[-1,0])...
             +circshift(tmpgr_min_x,[1,1])+circshift(tmpgr_add_x,[1,0]);
         mmxtmpJ3(:,:,i)=shift_mat;
         shift_mat=circshift(tmpgr_min_y,[0,-1])+circshift(tmpgr_add_y,[0,2])+circshift(tmpgr_min_y,[-1,1])+circshift(tmpgr_add_y,[-1,0])...
            +circshift(tmpgr_min_y,[1,1])+circshift(tmpgr_add_y,[1,0]);
         mmytmpJ3(:,:,i)=shift_mat;
         shift_mat=circshift(tmpgr_min_z,[0,-1])+circshift(tmpgr_add_z,[0,2])+circshift(tmpgr_min_z,[-1,1])+circshift(tmpgr_add_z,[-1,0])...
            +circshift(tmpgr_min_z,[1,1])+circshift(tmpgr_add_z,[1,0]);
         mmztmpJ3(:,:,i)=shift_mat;

         shift_mat=circshift(tmpp_add_x,[0,-1])+circshift(tmpp_add_x,[-1,1])+circshift(tmpp_add_x,[1,1])...
            +circshift(tmpgr_add_x,[0,-1])+circshift(tmpgr_add_x,[-1,1])+circshift(tmpgr_add_x,[1,1])...
            +circshift(tmpgr_min_x,[0,2]) +circshift(tmpgr_min_x,[-1,0]) +circshift(tmpgr_min_x,[1,0])...
             +circshift(tmpor_min_x,[0,1]) +circshift(tmpor_min_x,[-1,-1]) +circshift(tmpor_min_x,[1,-1]);
         mmxtmpJ4(:,:,i)=shift_mat;
         shift_mat=circshift(tmpp_add_y,[0,-1])+circshift(tmpp_add_y,[-1,1])+circshift(tmpp_add_y,[1,1])...
             +circshift(tmpgr_add_y,[0,-1])+circshift(tmpgr_add_y,[-1,1])+circshift(tmpgr_add_y,[1,1])...
            +circshift(tmpgr_min_y,[0,2]) +circshift(tmpgr_min_y,[-1,0]) +circshift(tmpgr_min_y,[1,0])...
             +circshift(tmpor_min_y,[0,1]) +circshift(tmpor_min_y,[-1,-1]) +circshift(tmpor_min_y,[1,-1]);
         mmytmpJ4(:,:,i)=shift_mat;
         shift_mat=circshift(tmpp_add_z,[0,-1])+circshift(tmpp_add_z,[-1,1])+circshift(tmpp_add_z,[1,1])...
            +circshift(tmpgr_add_z,[0,-1])+circshift(tmpgr_add_z,[-1,1])+circshift(tmpgr_add_z,[1,1])...
            +circshift(tmpgr_min_z,[0,2]) +circshift(tmpgr_min_z,[-1,0]) +circshift(tmpgr_min_z,[1,0])...
             +circshift(tmpor_min_z,[0,1]) +circshift(tmpor_min_z,[-1,-1]) +circshift(tmpor_min_z,[1,-1]);
         mmztmpJ4(:,:,i)=shift_mat;

        case 2
         mmxtmpJ1(:,:,i)=tmpgr_add_x+tmpp_min_x;
         mmytmpJ1(:,:,i)=tmpgr_add_y+tmpp_min_y;
         mmztmpJ1(:,:,i)=tmpgr_add_z+tmpp_min_z;

         shift_mat=circshift(tmpgr_x,[0,-1])+circshift(tmpp_x,[0,1])+circshift(tmpgr_x,[1,1])+circshift(tmpp_x,[1,-1])...%length direction -2 related to +2(right)
         +circshift(tmpgr_x,[-1,1])+circshift(tmpp_x,[-1,-1]);%length direction -2 related to +2(right)
         mmxtmpJ2(:,:,i)=shift_mat;
         shift_mat=circshift(tmpgr_y,[0,-1])+circshift(tmpp_y,[0,1])+circshift(tmpgr_y,[1,1])+circshift(tmpp_y,[1,-1])...%length direction -2 related to +2(right)
         +circshift(tmpgr_y,[-1,1])+circshift(tmpp_y,[-1,-1]);%length direction -2 related to +2(right)
         mmytmpJ2(:,:,i)=shift_mat;
         shift_mat=circshift(tmpgr_z,[0,-1])+circshift(tmpp_z,[0,1])+circshift(tmpgr_z,[1,1])+circshift(tmpp_z,[1,-1])...%length direction -2 related to +2(right)
         +circshift(tmpgr_z,[-1,1])+circshift(tmpp_z,[-1,-1]);%length direction -2 related to +2(right)
         mmztmpJ2(:,:,i)=shift_mat;

         shift_mat=circshift(tmpor_min_x,[0,-2])+circshift(tmpor_add_x,[0,1])+circshift(tmpor_min_x,[-1,0])+circshift(tmpor_add_x,[-1,-1])...
             +circshift(tmpor_min_x,[1,0])+circshift(tmpor_add_x,[1,-1]);
         mmxtmpJ3(:,:,i)=shift_mat;
         shift_mat=circshift(tmpor_min_y,[0,-2])+circshift(tmpor_add_y,[0,1])+circshift(tmpor_min_y,[-1,0])+circshift(tmpor_add_y,[-1,-1])...
             +circshift(tmpor_min_y,[1,0])+circshift(tmpor_add_y,[1,-1]);
         mmytmpJ3(:,:,i)=shift_mat;
         shift_mat=circshift(tmpor_min_z,[0,-2])+circshift(tmpor_add_z,[0,1])+circshift(tmpor_min_z,[-1,0])+circshift(tmpor_add_z,[-1,-1])...
          +circshift(tmpor_min_z,[1,0])+circshift(tmpor_add_z,[1,-1]);
         mmztmpJ3(:,:,i)=shift_mat;

         shift_mat=circshift(tmpgr_add_x,[0,-1])+circshift(tmpgr_add_x,[-1,1])+circshift(tmpgr_add_x,[1,1])...
            +circshift(tmpor_add_x,[1,0])+circshift(tmpor_add_x,[-1,0])+circshift(tmpor_add_x,[0,-2])...
            +circshift(tmpor_min_x,[0,1]) +circshift(tmpor_min_x,[-1,-1]) +circshift(tmpor_min_x,[1,-1])...
             +circshift(tmpp_min_x,[0,1]) +circshift(tmpp_min_x,[-1,-1]) +circshift(tmpp_min_x,[1,-1]);
         mmxtmpJ4(:,:,i)=shift_mat;
         shift_mat=circshift(tmpgr_add_y,[0,-1])+circshift(tmpgr_add_y,[-1,1])+circshift(tmpgr_add_y,[1,1])...
            +circshift(tmpor_add_y,[1,0])+circshift(tmpor_add_y,[-1,0])+circshift(tmpor_add_y,[0,-2])...
            +circshift(tmpor_min_y,[0,1]) +circshift(tmpor_min_y,[-1,-1]) +circshift(tmpor_min_y,[1,-1])...
             +circshift(tmpp_min_y,[0,1]) +circshift(tmpp_min_y,[-1,-1]) +circshift(tmpp_min_y,[1,-1]);       
         mmytmpJ4(:,:,i)=shift_mat;
         shift_mat=circshift(tmpgr_add_z,[0,-1])+circshift(tmpgr_add_z,[-1,1])+circshift(tmpgr_add_z,[1,1])...
            +circshift(tmpor_add_z,[1,0])+circshift(tmpor_add_z,[-1,0])+circshift(tmpor_add_z,[0,-2])...
            +circshift(tmpor_min_z,[0,1]) +circshift(tmpor_min_z,[-1,-1]) +circshift(tmpor_min_z,[1,-1])...
             +circshift(tmpp_min_z,[0,1]) +circshift(tmpp_min_z,[-1,-1]) +circshift(tmpp_min_z,[1,-1]);
         mmztmpJ4(:,:,i)=shift_mat;
        otherwise
            error('unexpected index');
    end
    end
         mmxtmpJzx_1(:,:,i)=mmxtmpJ1(:,:,i);
        mmytmpJzx_1(:,:,i)=mmytmpJ1(:,:,i);
        mmztmpJzx_1(:,:,i)=mmztmpJ1(:,:,i);
        mmxtmpJzx_2(:,:,i)=mmxtmpJ2(:,:,i);
        mmytmpJzx_2(:,:,i)=mmytmpJ2(:,:,i);
        mmztmpJzx_2(:,:,i)=mmztmpJ2(:,:,i);
        mmxtmpJzx_3(:,:,i)=mmxtmpJ3(:,:,i);
        mmytmpJzx_3(:,:,i)=mmytmpJ3(:,:,i);
        mmztmpJzx_3(:,:,i)=mmztmpJ3(:,:,i);
        mmxtmpJzx_4(:,:,i)=mmxtmpJ4(:,:,i);
        mmytmpJzx_4(:,:,i)=mmytmpJ4(:,:,i);
        mmztmpJzx_4(:,:,i)=mmztmpJ4(:,:,i);
end
1