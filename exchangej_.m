%% zf modify

for i=1:h 
     if i==1 || i==h
        b=0;
     else 
        b=1;
     end    
        %blue in the 2nd layer is identical to orange in 1st layer
        %green in the 1st layer is identical to red in 3rd(last) layer
        %red in the 3nd layer is identical to green in 1st layer

        %purple in the 2nd layer is identical to black in 3rd(last) layer
        %black in the 3nd layer is identical to purple in 2nd layer
        switch mod(i-1,3)
            case 0
                bb=0;
                %% calculate J1
               if i==1
                   bb=i+1;
               else 
                   bb=i;
               end
                mmxtmpJ1(:,:,i)=b*(atomtype_layer1gr.*mmxtmp(:,:,bb-1))+atomtype_layer1or.*mmxtmp(:,:,i+1);
                mmytmpJ1(:,:,i)=b*(atomtype_layer1gr.*mmytmp(:,:,bb-1))+atomtype_layer1or.*mmytmp(:,:,i+1);
                mmztmpJ1(:,:,i)=b*(atomtype_layer1gr.*mmztmp(:,:,bb-1))+atomtype_layer1or.*mmztmp(:,:,i+1);
                %% calculate J2
                %J2 experienced by green atoms

                shift_matx_gr=atomtype_layer1gr.*(mmx_0n2(:,:,i)+mmx_p10(:,:,i)...
                    +mmx_n10(:,:,i));
                shift_maty_gr=atomtype_layer1gr.*(mmy_0n2(:,:,i)+mmy_p10(:,:,i)...
                    +mmy_n10(:,:,i));
                shift_matz_gr=atomtype_layer1gr.*(mmz_0n2(:,:,i)+mmz_p10(:,:,i)...
                    +mmz_n10(:,:,i));
                %J2 experienced by orange atoms
                shift_matx_or=atomtype_layer1or.*(mmx_0p2(:,:,i)+mmx_p10(:,:,i)+...
                    mmx_n10(:,:,i));
                shift_maty_or=atomtype_layer1or.*(mmy_0p2(:,:,i)+mmy_p10(:,:,i)+...
                    mmy_n10(:,:,i));
                shift_matz_or=atomtype_layer1or.*(mmz_0p2(:,:,i)+mmz_p10(:,:,i)+...
                    mmz_n10(:,:,i));

                mmxtmpJ2(:,:,i)=shift_matx_gr+shift_matx_or;
                mmytmpJ2(:,:,i)=shift_maty_gr+shift_maty_or;
                mmztmpJ2(:,:,i)=shift_matz_gr+shift_matz_or;

                %% calculate J3
                %J3 experienced by green atoms
                shift_matx_gr=atomtype_layer1gr.*(mmx_0p1(:,:,i+1)+mmx_n1n1(:,:,i+1)+...
                mmx_p1n1(:,:,i+1));
                shift_maty_gr=atomtype_layer1gr.*(mmy_0p1(:,:,i+1)+mmy_n1n1(:,:,i+1)+...
                mmy_p1n1(:,:,i+1));
                shift_matz_gr=atomtype_layer1gr.*(mmz_0p1(:,:,i+1)+mmz_n1n1(:,:,i+1)+...
                mmz_p1n1(:,:,i+1));
                %J3 experienced by orange atoms
                %black in the 3nd layer is identical to purple in 2nd layer
                shift_matx_or=b*(atomtype_layer1or.*(mmx_0n1(:,:,bb-1)+mmx_n1p1(:,:,bb-1)+...
                mmx_p1p1(:,:,bb-1)));
                shift_maty_or=b*(atomtype_layer1or.*(mmy_0n1(:,:,bb-1)+mmy_n1p1(:,:,bb-1)+...
                mmy_p1p1(:,:,bb-1)));
                shift_matz_or=b*(atomtype_layer1or.*(mmz_0n1(:,:,bb-1)+mmz_n1p1(:,:,bb-1)+...
                mmz_p1p1(:,:,bb-1)));

                mmxtmpJ3(:,:,i)=shift_matx_gr+shift_matx_or;
                mmytmpJ3(:,:,i)=shift_maty_gr+shift_maty_or;
                mmztmpJ3(:,:,i)=shift_matz_gr+shift_matz_or;

                %% calculate J4
                %J4 experienced by green atoms
                %blue in the 2nd layer is identical to orange in 1st layer
                %black in the 3nd layer is identical to purple in 2nd layer
                shift_matx_gr=atomtype_layer1gr.*(mmx_0n2(:,:,i+1)+...
                    mmx_p10(:,:,i+1)+mmx_n10(:,:,i+1)+b*(mmx_0p1(:,:,bb-1)+...
                    mmx_n1n1(:,:,bb-1)+mmx_p1n1(:,:,bb-1)));
                shift_maty_gr=atomtype_layer1gr.*(mmy_0n2(:,:,i+1)+...
                    mmy_p10(:,:,i+1)+mmy_n10(:,:,i+1)+b*(mmy_0p1(:,:,bb-1)+...
                    mmy_n1n1(:,:,bb-1)+mmy_p1n1(:,:,bb-1)));
                shift_matz_gr=atomtype_layer1gr.*(mmz_0n2(:,:,i+1)+...
                    mmz_p10(:,:,i+1)+mmz_n10(:,:,i+1)+b*(mmz_0p1(:,:,bb-1)+...
                    mmz_n1n1(:,:,bb-1)+mmz_p1n1(:,:,bb-1)));
                %J4 experienced by orange atoms
                %red in the 3nd layer is identical to green in 1st layer
                shift_matx_or=atomtype_layer1or.*(mmx_0n1(:,:,i+1)+...
                    mmx_n1p1(:,:,i+1)+mmx_p1p1(:,:,i+1)+b*(mmx_0p2(:,:,bb-1)+...
                    mmx_n10(:,:,bb-1)+mmx_p10(:,:,bb-1)));
                shift_maty_or=atomtype_layer1or.*(mmy_0n1(:,:,i+1)+...
                    mmy_n1p1(:,:,i+1)+mmy_p1p1(:,:,i+1)+ b*(mmy_0p2(:,:,bb-1)+...
                    mmy_n10(:,:,bb-1)+mmy_p10(:,:,bb-1)));
                shift_matz_or=atomtype_layer1or.*(mmz_0n1(:,:,i+1)+...
                    mmz_n1p1(:,:,i+1)+mmz_p1p1(:,:,i+1)+b*(mmz_0p2(:,:,bb-1)+...
                    mmz_n10(:,:,bb-1)+mmz_p10(:,:,bb-1)));

                mmxtmpJ4(:,:,i)=shift_matx_gr+shift_matx_or;
                mmytmpJ4(:,:,i)=shift_maty_gr+shift_maty_or;
                mmztmpJ4(:,:,i)=shift_matz_gr+shift_matz_or;

            case 1
                %% calculate J1
                mmxtmpJ1(:,:,i)=atomtype_layer2p.*mmxtmp(:,:,i+1)+atomtype_layer2blue.*mmxtmp(:,:,i-1);
                mmytmpJ1(:,:,i)=atomtype_layer2p.*mmytmp(:,:,i+1)+atomtype_layer2blue.*mmytmp(:,:,i-1);
                mmztmpJ1(:,:,i)=atomtype_layer2p.*mmztmp(:,:,i+1)+atomtype_layer2blue.*mmztmp(:,:,i-1);

                %% calculate J2
                %J2 experienced by blue atoms
                shift_matx_blue=atomtype_layer2blue.*(mmx_0n1(:,:,i)+mmx_n1p1(:,:,i)+...
                    mmx_p1p1(:,:,i));
                shift_maty_blue=atomtype_layer2blue.*(mmy_0n1(:,:,i)+mmy_n1p1(:,:,i)+...
                    mmy_p1p1(:,:,i));
                shift_matz_blue=atomtype_layer2blue.*(mmz_0n1(:,:,i)+mmz_n1p1(:,:,i)+...
                    mmz_p1p1(:,:,i));

                %J2 experienced by purple atoms
                %blue in the 2nd layer is identical to orange in 1st layer
                shift_matx_purple=atomtype_layer2p.*(mmx_0p1(:,:,i)+mmx_n1n1(:,:,i)+...
                     mmx_p1n1(:,:,i));
                shift_maty_purple=atomtype_layer2p.*(mmy_0p1(:,:,i)+mmy_n1n1(:,:,i)+...
                     mmy_p1n1(:,:,i));
                shift_matz_purple=atomtype_layer2p.*(mmz_0p1(:,:,i)+mmz_n1n1(:,:,i)+...
                     mmz_p1n1(:,:,i));

                mmxtmpJ2(:,:,i)=shift_matx_blue+shift_matx_purple;
                mmytmpJ2(:,:,i)=shift_maty_blue+shift_maty_purple;
                mmztmpJ2(:,:,i)=shift_matz_blue+shift_matz_purple;

                %% calculate J3
                %J3 experienced by blue atoms
                %red in the 3rd layer is identical to green in 1st layer
                shift_matx_blue=atomtype_layer2blue.*(mmx_0p2(:,:,i+1)+mmx_n10(:,:,i+1)+...
                mmx_p10(:,:,i+1));
                shift_maty_blue=atomtype_layer2blue.*(mmy_0p2(:,:,i+1)+mmy_n10(:,:,i+1)+...
                mmy_p10(:,:,i+1));
                shift_matz_blue=atomtype_layer2blue.*(mmz_0p2(:,:,i+1)+mmz_n10(:,:,i+1)+...
                mmz_p10(:,:,i+1));
                %J3 experienced by purple atoms
                shift_matx_purple=atomtype_layer2p.*(mmx_0n1(:,:,i-1)+mmx_n1p1(:,:,i-1)+...
                mmx_p1p1(:,:,i-1));
                shift_maty_purple=atomtype_layer2p.*(mmy_0n1(:,:,i-1)+mmy_n1p1(:,:,i-1)+...
                mmy_p1p1(:,:,i-1));
                shift_matz_purple=atomtype_layer2p.*(mmz_0n1(:,:,i-1)+mmz_n1p1(:,:,i-1)+...
                mmz_p1p1(:,:,i-1));

                mmxtmpJ3(:,:,i)=shift_matx_blue+shift_matx_purple;
                mmytmpJ3(:,:,i)=shift_maty_blue+shift_maty_purple;
                mmztmpJ3(:,:,i)=shift_matz_blue+shift_matz_purple;

                %% calculate J4
                %J4 experienced by blue atoms
                %black in the 3nd layer is identical to purple in 2nd layer
                shift_matx_blue=atomtype_layer2blue.*(mmx_0n1(:,:,i+1)+...
                mmx_n1p1(:,:,i+1)+mmx_p1p1(:,:,i+1)+mmx_0p2(:,:,i-1)+...
                mmx_n10(:,:,i-1)+mmx_p10(:,:,i-1));
                shift_maty_black=atomtype_layer2blue.*(mmy_0n1(:,:,i+1)+...
                mmy_n1p1(:,:,i+1)+mmy_p1p1(:,:,i+1)+mmy_0p2(:,:,i-1)+...
                mmy_n10(:,:,i-1)+mmy_p10(:,:,i-1));
                shift_matz_black=atomtype_layer2blue.*(mmz_0n1(:,:,i+1)+...
                mmz_n1p1(:,:,i+1)+mmz_p1p1(:,:,i+1)+mmz_0p2(:,:,i-1)+...
                mmz_n10(:,:,i-1)+mmz_p10(:,:,i-1));
                %J4 experienced by purple atoms
                %red in the 3nd layer is identical to green in 1st layer
                shift_matx_purple=atomtype_layer2p.*(mmx_0n1(:,:,i+1)+...
                mmx_n1p1(:,:,i+1)+mmx_p1p1(:,:,i+1)+mmx_0p1(:,:,i-1)+...
                mmx_n1n1(:,:,i-1)+mmx_p1n1(:,:,i-1));
                shift_maty_purple=atomtype_layer2p.*(mmy_0n1(:,:,i+1)+...
                mmy_n1p1(:,:,i+1)+mmy_p1p1(:,:,i+1)+mmy_0p1(:,:,i-1)+...
                mmy_n1n1(:,:,i-1)+mmy_p1n1(:,:,i-1));
                shift_matz_purple=atomtype_layer2p.*(mmz_0n1(:,:,i+1)+...
                mmz_n1p1(:,:,i+1)+mmz_p1p1(:,:,i+1)+mmz_0p1(:,:,i-1)+...
                mmz_n1n1(:,:,i-1)+mmz_p1n1(:,:,i-1));

                mmxtmpJ4(:,:,i)=shift_matx_blue+shift_matx_purple;
                mmytmpJ4(:,:,i)=shift_maty_black+shift_maty_purple;
                mmztmpJ4(:,:,i)=shift_matz_black+shift_matz_purple;

            case 2
                hh=0;
               if i==h
                   hh=i-1;
               else 
                   hh=i;
               end
                %% calculate J1
                mmxtmpJ1(:,:,i)=b*(atomtype_layer3red.*mmxtmp(:,:,hh+1))+atomtype_layer3black.*mmxtmp(:,:,i-1);  
                mmytmpJ1(:,:,i)=b*(atomtype_layer3red.*mmytmp(:,:,hh+1))+atomtype_layer3black.*mmytmp(:,:,i-1);  
                mmztmpJ1(:,:,i)=b*(atomtype_layer3red.*mmztmp(:,:,hh+1))+atomtype_layer3black.*mmztmp(:,:,i-1);   

                %% calculate J2
                %J2 experienced by black atoms
                %green in the 1st layer is identical to red in 3rd(last) layer
                shift_matx_black=atomtype_layer3black.*(mmx_0n1(:,:,i)+mmx_n1p1(:,:,i)+...
                    mmx_p1p1(:,:,i));
                shift_maty_black=atomtype_layer3black.*(mmy_0n1(:,:,i)+mmy_n1p1(:,:,i)+...
                    mmy_p1p1(:,:,i));
                shift_matz_black=atomtype_layer3black.*(mmz_0n1(:,:,i)+mmz_n1p1(:,:,i)+...
                    mmz_p1p1(:,:,i));
                %J2 experienced by red atoms
                %purple in the 2nd layer is identical to black in 3rd(last) layer
                shift_matx_red=atomtype_layer3red.*(mmx_0p1(:,:,i)+mmx_n1n1(:,:,i)+...
                    mmx_p1n1(:,:,i));
                shift_maty_red=atomtype_layer3red.*(mmy_0p1(:,:,i)+mmy_n1n1(:,:,i)+...
                    mmy_p1n1(:,:,i));
                shift_matz_red=atomtype_layer3red.*(mmz_0p1(:,:,i)+mmz_n1n1(:,:,i)+...
                    mmz_p1n1(:,:,i));

                mmxtmpJ2(:,:,i)=shift_matx_black+shift_matx_red;
                mmytmpJ2(:,:,i)=shift_maty_black+shift_maty_red;
                mmztmpJ2(:,:,i)=shift_matz_black+shift_matz_red;

                %% calculate J3
                %J3 for red atoms 
                shift_matx_red=atomtype_layer3red.*(mmx_0n2(:,:,i-1)+mmx_n10(:,:,i-1)+...
                mmx_p10(:,:,i-1));
                shift_maty_red=atomtype_layer3red.*(mmy_0n2(:,:,i-1)+mmy_n10(:,:,i-1)+...
                mmy_p10(:,:,i-1));
                shift_matz_red=atomtype_layer3red.*(mmz_0n2(:,:,i-1)+mmz_n10(:,:,i-1)+...
                mmz_p10(:,:,i-1));
                %J3 for black atoms
                shift_matx_black=atomtype_layer3black.*(mmx_0p1(:,:,hh+1)+mmx_n1n1(:,:,hh+1)+...
                mmx_p1n1(:,:,hh+1));
                shift_maty_black=atomtype_layer3black.*(mmy_0p1(:,:,hh+1)+mmy_n1n1(:,:,hh+1)+...
                mmy_p1n1(:,:,hh+1));
                shift_matz_black=atomtype_layer3black.*(mmz_0p1(:,:,hh+1)+mmz_n1n1(:,:,hh+1)+...
                mmz_p1n1(:,:,hh+1));

                mmxtmpJ3(:,:,i)=shift_matx_red+b*shift_matx_black;
                mmytmpJ3(:,:,i)=shift_maty_red+b*shift_maty_black;
                mmztmpJ3(:,:,i)=shift_matz_red+b*shift_matz_black;

                %% calculate J4
                %J4 experienced by black atoms
                %blue in the 2nd layer is identical to orange in 1st layer
                shift_matx_black=atomtype_layer3black.*(b*(mmx_0n1(:,:,hh+1)+...
                    mmx_n1p1(:,:,hh+1)+mmx_p1p1(:,:,hh+1))+mmx_0p1(:,:,i-1)+...
                    mmx_n1n1(:,:,i-1)+mmx_p1n1(:,:,i-1));
                shift_maty_black=atomtype_layer3black.*(b*(mmy_0n1(:,:,hh+1)+...
                    mmy_n1p1(:,:,hh+1)+mmy_p1p1(:,:,hh+1))+mmy_0p1(:,:,i-1)+...
                    mmy_n1n1(:,:,i-1)+mmy_p1n1(:,:,i-1));
                shift_matz_black=atomtype_layer3black.*(b*(mmz_0n1(:,:,hh+1)+...
                    mmz_n1p1(:,:,hh+1)+mmz_p1p1(:,:,hh+1))+mmz_0p1(:,:,i-1)+...
                    mmz_n1n1(:,:,i-1)+mmz_p1n1(:,:,i-1));
                %J4 experienced by red atoms
                %blue in the 2nd layer is identical to orange in 1st layer
                shift_matx_red=atomtype_layer3red.*(b*(mmx_p10(:,:,hh+1)+...
                    mmx_n10(:,:,hh+1)+mmx_0n2(:,:,hh+1))+mmx_0p1(:,:,i-1)+...
                    mmx_n1n1(:,:,i-1)+mmx_p1n1(:,:,i-1));
                shift_maty_red=atomtype_layer3red.*(b*(mmy_p10(:,:,hh+1)+...
                    mmy_n10(:,:,hh+1)+mmy_0n2(:,:,hh+1))+mmy_0p1(:,:,i-1)+mmy_n1n1(:,:,i-1)+...
                    mmy_p1n1(:,:,i-1));
                shift_matz_red=atomtype_layer3red.*(b*(mmz_p10(:,:,hh+1)+...
                    mmz_n10(:,:,hh+1)+ mmz_0n2(:,:,hh+1))+mmz_0p1(:,:,i-1)+...
                    mmz_n1n1(:,:,i-1)+mmz_p1n1(:,:,i-1));

                mmxtmpJ4(:,:,i)=shift_matx_black+shift_matx_red;
                mmytmpJ4(:,:,i)=shift_maty_black+shift_maty_red;
                mmztmpJ4(:,:,i)=shift_matz_black+shift_matz_red;
            otherwise
                error('unexpected index');
        end
  end

