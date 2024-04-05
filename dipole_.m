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
            dot_sr=muigpu.*(atomtype_layer1gr.*mmx_0n2(:,:,i).*rijx_+atomtype_layer1gr.*mmy_0n2(:,:,i).*rijy_);
            hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmx_0n2(:,:,i))./dist_.^3);%[T]
            hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmy_0n2(:,:,i))./dist_.^3);
            hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1gr.*mmz_0n2(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer1gr.*mmx_n10(:,:,i).*rijx_+atomtype_layer1gr.*mmy_n10(:,:,i).*rijy_);
            hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmx_n10(:,:,i))./dist_.^3);%[T]
            hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmy_n10(:,:,i))./dist_.^3);
            hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1gr.*mmz_n10(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer1gr.*mmx_p10(:,:,i).*rijx_+atomtype_layer1gr.*mmy_p10(:,:,i).*rijy_);
            hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmx_p10(:,:,i))./dist_.^3);%[T]
            hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1gr.*mmy_p10(:,:,i))./dist_.^3);
            hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1gr.*mmz_p10(:,:,i))./dist_.^3);

            %%the dipole field for the orange atoms
            rijx_=0.5067*1e-9*-1;rijy_=0;
            dot_sr=muigpu.*(atomtype_layer1or.*mmx_0p2(:,:,i).*rijx_+atomtype_layer1or.*mmy_0p2(:,:,i).*rijy_);
            hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmx_0p2(:,:,i))./dist_.^3);%[T]
            hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmy_0p2(:,:,i))./dist_.^3);
            hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1or.*mmz_0p2(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer1or.*mmx_n10(:,:,i).*rijx_+mmy_n10(:,:,i).*atomtype_layer1or.*rijy_);
            hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmx_n10(:,:,i))./dist_.^3);%[T]
            hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmy_n10(:,:,i))./dist_.^3);
            hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1or.*mmz_n10(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer1or.*mmx_p10(:,:,i).*rijx_+atomtype_layer1or.*mmy_p10(:,:,i).*rijy_);
            hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmx_p10(:,:,i))./dist_.^3);%[T]
            hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer1or.*mmy_p10(:,:,i))./dist_.^3);
            hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer1or.*mmz_p10(:,:,i))./dist_.^3);

            hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
            hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
            hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);

        case 1
            %%type layer2 bule and purple atoms
            %%dipole field of blue atoms
            rijx_=0.5067*1e-9;rijy_=0;
            dot_sr=muigpu.*(atomtype_layer2blue.*mmx_0n1(:,:,i).*rijx_+atomtype_layer2blue.*mmy_0n1(:,:,i).*rijy_);
            hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmx_0n1(:,:,i))./dist_.^3);%[T]
            hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmy_0n1(:,:,i))./dist_.^3);
            hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2blue.*mmz_0n1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer2blue.*mmx_n1p1(:,:,i).*rijx_+atomtype_layer2blue.*mmy_n1p1(:,:,i).*rijy_);
            hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmx_n1p1(:,:,i))./dist_.^3);%[T]
            hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmy_n1p1(:,:,i))./dist_.^3);
            hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2blue.*mmz_n1p1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer2blue.*mmx_p1p1(:,:,i)).*rijx_+(atomtype_layer2blue.*mmy_p1p1(:,:,i)).*rijy_);
            hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmx_p1p1(:,:,i))./dist_.^3);%[T]
            hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2blue.*mmy_p1p1(:,:,i))./dist_.^3);
            hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2blue.*mmz_p1p1(:,:,i))./dist_.^3);

            %% dipole field for purple atoms
            rijx_=-0.5067*1e-9;rijy_=0;
            dot_sr=muigpu.*(atomtype_layer2p.*mmx_0p1(:,:,i).*rijx_+atomtype_layer2p.*mmy_0p1(:,:,i).*rijy_);
            hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmx_0p1(:,:,i))./dist_.^3);%[T]
            hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmy_0p1(:,:,i))./dist_.^3);
            hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2p.*mmz_0p1(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*(atomtype_layer2p.*mmx_n1n1(:,:,i).*rijx_+atomtype_layer2p.*mmy_n1n1(:,:,i).*rijy_);
            hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmx_n1n1(:,:,i))./dist_.^3);%[T]
            hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmy_n1n1(:,:,i))./dist_.^3);
            hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2p.*mmz_n1n1(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer2p.*mmx_p1n1(:,:,i)).*rijx_+(atomtype_layer2p.*mmy_p1n1(:,:,i)).*rijy_);
            hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmx_p1n1(:,:,i))./dist_.^3);%[T]
            hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer2p.*mmy_p1n1(:,:,i))./dist_.^3);
            hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer2p.*mmz_p1n1(:,:,i))./dist_.^3);


            hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
            hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
            hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);

        case 2
            %%type layer3 red and black atoms
            %%dipole field of black atoms
            rijx_=0.5067*1e-9;rijy_=0;
            dot_sr=muigpu.*((atomtype_layer3black.*mmx_0n1(:,:,i)).*rijx_+(atomtype_layer3black.*mmy_0n1(:,:,i)).*rijy_);
            hdipolex1(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmx_0n1(:,:,i))./dist_.^3);%[T]
            hdipoley1(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmy_0n1(:,:,i))./dist_.^3);
            hdipolez1(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3black.*mmz_0n1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer3black.*mmx_n1p1(:,:,i)).*rijx_+(atomtype_layer3black.*mmy_n1p1(:,:,i)).*rijy_);
            hdipolex2(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmx_n1p1(:,:,i))./dist_.^3);%[T]
            hdipoley2(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmy_n1p1(:,:,i))./dist_.^3);
            hdipolez2(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3black.*mmz_n1p1(:,:,i))./dist_.^3);

            rijx_=-0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer3black.*mmx_p1p1(:,:,i)).*rijx_+(atomtype_layer3black.*mmy_p1p1(:,:,i)).*rijy_);
            hdipolex3(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmx_p1p1(:,:,i))./dist_.^3);%[T]
            hdipoley3(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3black.*mmy_p1p1(:,:,i))./dist_.^3);
            hdipolez3(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3black.*mmz_p1p1(:,:,i))./dist_.^3);

            %% dipole field for red atoms
            rijx_=-0.5067*1e-9;rijy_=0;
            dot_sr=muigpu.*((atomtype_layer3red.*mmx_0p1(:,:,i)).*rijx_+(atomtype_layer3red.*mmy_0p1(:,:,i)).*rijy_);
            hdipolex4(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmx_0p1(:,:,i))./dist_.^3);%[T]
            hdipoley4(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmy_0p1(:,:,i))./dist_.^3);
            hdipolez4(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3red.*mmz_0p1(:,:,i))./dist_.^3);

            rijx_=0.25335*1e-9;rijy_=-0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer3red.*mmx_n1n1(:,:,i)).*rijx_+(atomtype_layer3red.*mmy_n1n1(:,:,i)).*rijy_);
            hdipolex5(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmx_n1n1(:,:,i))./dist_.^3);%[T]
            hdipoley5(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmy_n1n1(:,:,i))./dist_.^3);
            hdipolez5(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3red.*mmz_n1n1(:,:,i))./dist_.^3);


            rijx_=0.25335*1e-9;rijy_=0.4388*1e-9;
            dot_sr=muigpu.*((atomtype_layer3red.*mmx_p1n1(:,:,i)).*rijx_+(atomtype_layer3red.*mmy_p1n1(:,:,i)).*rijy_);
            hdipolex6(:,:,i)=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmx_p1n1(:,:,i))./dist_.^3);%[T]
            hdipoley6(:,:,i)=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*(atomtype_layer3red.*mmy_p1n1(:,:,i))./dist_.^3);
            hdipolez6(:,:,i)=mu_0/(4*pi)*(-muigpu.*(atomtype_layer3red.*mmz_p1n1(:,:,i))./dist_.^3);


            hdipolez_(:,:,i)= hdipolez6(:,:,i)+hdipolez5(:,:,i)+hdipolez4(:,:,i)+hdipolez3(:,:,i)+hdipolez2(:,:,i)+hdipolez1(:,:,i);
            hdipoley_(:,:,i)= hdipoley6(:,:,i)+hdipoley5(:,:,i)+hdipoley4(:,:,i)+hdipoley3(:,:,i)+hdipoley2(:,:,i)+hdipoley1(:,:,i);
            hdipolex_(:,:,i)= hdipolex6(:,:,i)+hdipolex5(:,:,i)+hdipolex4(:,:,i)+hdipolex3(:,:,i)+hdipolex2(:,:,i)+hdipolex1(:,:,i);

        otherwise
            error('unexpected index');

    end
    hdipo_x(:,:,i)= hdipolex_(:,:,i); hdipo_y(:,:,i)= hdipoley_(:,:,i); hdipo_z(:,:,i)= hdipolez_(:,:,i);
end
clear hdipolex_ hdipoley_ hdipolez_

