%%%dipole part  2024_3_26
%since the lattice constant is a=0.5067e-9, c= 13.88e-10.
% we only consider the dipole in a distance a=0.5067e-9.
%three different cases as the exchange J2.
for i=1:natomH
    dist_=0.5067*1e-9;%% the distance between two nearst atoms in the lattice plane is a=0.5067e-9.
    switch mod(i-1,3)
        case 0
            %type layer1, we have two different atoms, i.e., green and orange
            %dipole field for green atoms

            rijx_=0.5067*1e-9*1;rijy_=0;% rijx_ and rijy_ is the true distance in the lattice
            dot_sr=muigpu.*(atomtype_layer1gr.*mmx_l2(:,:,i).*rijx_+atomtype_layer1gr.*mmy_l2(:,:,i).*rijy_);
            hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmx_l2(:,:,i))./dist_.^3);%[T]
            hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmy_l2(:,:,i))./dist_.^3);
            hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1gr.*mmz_l2(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer1gr.*mmx_u1(:,:,i).*rijx_+atomtype_layer1gr.*mmy_u1(:,:,i).*rijy_);
            hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmx_u1(:,:,i))./dist_.^3);%[T]
            hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmy_u1(:,:,i))./dist_.^3);
            hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1gr.*mmz_u1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer1gr.*mmx_d1(:,:,i).*rijx_+atomtype_layer1gr.*mmy_d1(:,:,i).*rijy_);
            hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmx_d1(:,:,i))./dist_.^3);%[T]
            hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmy_d1(:,:,i))./dist_.^3);
            hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1gr.*mmz_d1(:,:,i))./dist_.^3);

            %%the dipole field for the orange atoms
            rijx_=0.5067*1e-9*-1;rijy_=0;
            dot_sr=muigpu.*(atomtype_layer1or.*mmx_r2(:,:,i).*rijx_+atomtype_layer1or.*mmy_r2(:,:,i).*rijy_);
            hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmx_r2(:,:,i))./dist_.^3);%[T]
            hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmy_r2(:,:,i))./dist_.^3);
            hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1or.*mmz_r2(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer1or.*mmx_u1(:,:,i).*rijx_+mmy_u1(:,:,i).*atomtype_layer1or.*rijy_);
            hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmx_u1(:,:,i))./dist_.^3);%[T]
            hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmy_u1(:,:,i))./dist_.^3);
            hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1or.*mmz_u1(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer1or.*mmx_d1(:,:,i).*rijx_+atomtype_layer1or.*mmy_d1(:,:,i).*rijy_);
            hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmx_d1(:,:,i))./dist_.^3);%[T]
            hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmy_d1(:,:,i))./dist_.^3);
            hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1or.*mmz_d1(:,:,i))./dist_.^3);

            hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
            hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
            hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);

        case 1
            %%type layer2 bule and purple atoms
            %%dipole field of blue atoms
            rijx_=0.5067*1e-9;rijy_=0;
            dot_sr=muigpu.*(atomtype_layer2blue.*mmx_l1(:,:,i).*rijx_+atomtype_layer2blue.*mmy_l1(:,:,i).*rijy_);
            hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmx_l1(:,:,i))./dist_.^3);%[T]
            hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmy_l1(:,:,i))./dist_.^3);
            hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2blue.*mmz_l1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer2blue.*mmx_u1r1(:,:,i).*rijx_+atomtype_layer2blue.*mmy_u1r1(:,:,i).*rijy_);
            hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmx_u1r1(:,:,i))./dist_.^3);%[T]
            hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmy_u1r1(:,:,i))./dist_.^3);
            hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2blue.*mmz_u1r1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer2blue.*mmx_d1r1(:,:,i)).*rijx_+(atomtype_layer2blue.*mmy_d1r1(:,:,i)).*rijy_);
            hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmx_d1r1(:,:,i))./dist_.^3);%[T]
            hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmy_d1r1(:,:,i))./dist_.^3);
            hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2blue.*mmz_d1r1(:,:,i))./dist_.^3);

            %% dipole field for purple atoms
            rijx_=-0.5067*1e-9;rijy_=0;
            dot_sr=muigpu.*(atomtype_layer2p.*mmx_r1(:,:,i).*rijx_+atomtype_layer2p.*mmy_r1(:,:,i).*rijy_);
            hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmx_r1(:,:,i))./dist_.^3);%[T]
            hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmy_r1(:,:,i))./dist_.^3);
            hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2p.*mmz_r1(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer2p.*mmx_u1l1(:,:,i).*rijx_+atomtype_layer2p.*mmy_u1l1(:,:,i).*rijy_);
            hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmx_u1l1(:,:,i))./dist_.^3);%[T]
            hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmy_u1l1(:,:,i))./dist_.^3);
            hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2p.*mmz_u1l1(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer2p.*mmx_d1l1(:,:,i)).*rijx_+(atomtype_layer2p.*mmy_d1l1(:,:,i)).*rijy_);
            hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmx_d1l1(:,:,i))./dist_.^3);%[T]
            hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmy_d1l1(:,:,i))./dist_.^3);
            hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2p.*mmz_d1l1(:,:,i))./dist_.^3);


            hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
            hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
            hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);

        case 2
            %%type layer3 red and black atoms
            %%dipole field of black atoms
            rijx_=0.5067*1e-9;rijy_=0;
            dot_sr=muigpu.*((atomtype_layer3black.*mmx_l1(:,:,i)).*rijx_+(atomtype_layer3black.*mmy_l1(:,:,i)).*rijy_);
            hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmx_l1(:,:,i))./dist_.^3);%[T]
            hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmy_l1(:,:,i))./dist_.^3);
            hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3black.*mmz_l1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer3black.*mmx_u1r1(:,:,i)).*rijx_+(atomtype_layer3black.*mmy_u1r1(:,:,i)).*rijy_);
            hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmx_u1r1(:,:,i))./dist_.^3);%[T]
            hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmy_u1r1(:,:,i))./dist_.^3);
            hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3black.*mmz_u1r1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer3black.*mmx_d1r1(:,:,i)).*rijx_+(atomtype_layer3black.*mmy_d1r1(:,:,i)).*rijy_);
            hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmx_d1r1(:,:,i))./dist_.^3);%[T]
            hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmy_d1r1(:,:,i))./dist_.^3);
            hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3black.*mmz_d1r1(:,:,i))./dist_.^3);

            %% dipole field for red atoms
            rijx_=-0.5067*1e-9;rijy_=0;
            dot_sr=muigpu.*((atomtype_layer3red.*mmx_r1(:,:,i)).*rijx_+(atomtype_layer3red.*mmy_r1(:,:,i)).*rijy_);
            hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmx_r1(:,:,i))./dist_.^3);%[T]
            hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmy_r1(:,:,i))./dist_.^3);
            hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3red.*mmz_r1(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer3red.*mmx_u1l1(:,:,i)).*rijx_+(atomtype_layer3red.*mmy_u1l1(:,:,i)).*rijy_);
            hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmx_u1l1(:,:,i))./dist_.^3);%[T]
            hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmy_u1l1(:,:,i))./dist_.^3);
            hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3red.*mmz_u1l1(:,:,i))./dist_.^3);


            rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer3red.*mmx_d1l1(:,:,i)).*rijx_+(atomtype_layer3red.*mmy_d1l1(:,:,i)).*rijy_);
            hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmx_d1l1(:,:,i))./dist_.^3);%[T]
            hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmy_d1l1(:,:,i))./dist_.^3);
            hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3red.*mmz_d1l1(:,:,i))./dist_.^3);


            hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
            hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
            hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);

        otherwise
            error('unexpected index');

    end
    hdipo_x(:,:,i)= hdipolex_(:,:,i); hdipo_y(:,:,i)= hdipoley_(:,:,i); hdipo_z(:,:,i)= hdipolez_(:,:,i);
end
clear hdipolex_ hdipoley_ hdipolez_

