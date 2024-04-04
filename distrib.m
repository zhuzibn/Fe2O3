%% distribute atoms
%1 RE,0 TM

atomtype_ = 2*ones(natomW, natomL, natomH);
for ctW=1:natomW
for ctL = 1:natomL
    for ctH = 1:natomH
        if (mod(ctH-1,3)==0)
            if mod(ctL, 2) == 1  % 对于偶数 ctL
                atomtype_(ctW, ctL, ctH) = 1;
            end
        end
    end
end
end
for ctW = 1:natomW
    for ctH = 1:natomH
        if (mod(ctH-1,3)==1)
            if mod(ctW, 2) == 1  % 对于奇数 ctW
                startL = 3;
                start2L = 4;
            else  % 对于偶数 ctW
                startL = 1;
                start2L = 2;
            end

            for ctL = startL:4:natomL
                atomtype_(ctW, ctL, ctH) = 1;
            end
            for ctL = start2L:4:natomL
                atomtype_(ctW, ctL, ctH) = 1;
            end
     
        end
    end
end
for ctW = 1:natomW
    for ctH = 1:natomH
        if (mod(ctH-1,3)==2)
            if mod(ctW, 2) == 1  % 对于奇数 ctW
                startL = 1;
                start2L = 4;
            else  % 对于偶数 ctW
                startL = 2;
                start2L = 3;
            end

            for ctL = startL:4:natomL
                atomtype_(ctW, ctL, ctH) = 1;
            end
            for ctL = start2L:4:natomL
                atomtype_(ctW, ctL, ctH) = 1;
            end
        end
    end
end

new_natomW = natomW + 4;
new_natomL = natomL + 4;

new_atomtype_ = 2 * ones(new_natomW, new_natomL, natomH);
for ctH = 1:natomH
    new_atomtype_(3:end-2, 3:end-2, ctH) = atomtype_(:, :, ctH);
end
%% filter for each layer
atomtype_layer1gr = zeros(natomW,natomL,'gpuArray');
atomtype_layer1or = zeros(natomW,natomL,'gpuArray');
atomtype_layer2p = zeros(natomW,natomL,'gpuArray');

%% layer 1gre=layer 3r; layer 1or=layer b; layer 2p=layer 3black;
for ctW = 1:natomW
            if mod(ctW, 2) == 1  % 对于奇数 ctW
                startL = 1;
                start2L = 3;
                start3L = 4;
            else  % 对于偶数 ctW
                startL = 3;   
                start2L = 1;
                start3L = 2;
            end
            for  ctL = startL:4:natomL
                atomtype_layer1gr(ctW, ctL) = 1;
            end   
             for ctL = start2L:4:natomL
                atomtype_layer1or(ctW, ctL) = 1;  
             end 
            for ctL = start3L:4:natomL 
                atomtype_layer2p(ctW, ctL) = 1;
            end 
end
new_atomtype_laygr = 2 * ones(new_natomW,new_natomL,'gpuArray');
new_atomtype_layor = 2 * ones(new_natomW,new_natomL,'gpuArray');
new_atomtype_layp = 2 * ones(new_natomW,new_natomL,'gpuArray');

new_atomtype_laygr (3:end-2, 3:end-2) = atomtype_layer1gr(:, :);
new_atomtype_layor (3:end-2, 3:end-2) = atomtype_layer1or(:, :);
new_atomtype_layp (3:end-2, 3:end-2) = atomtype_layer2p(:, :);


atomtype_=new_atomtype_;
atomtype_layer2p =new_atomtype_layp;
atomtype_layer1or =new_atomtype_layor;
atomtype_layer1gr=new_atomtype_laygr;
atomtype_s= zeros(new_natomW,new_natomL,natomH,'gpuArray');
ato_s= zeros(new_natomW,new_natomL,natomH,'gpuArray');
%%atomtype_s is the filter of atoms in different layer.0= no atom, 1=have atom
%%ato_s is just the opposite, this matrix is used to avoid the NaN result
%%in the rk4 calculation. 0=have atom, 1=no atom.
for i=1:new_natomW
    for j=1:new_natomL
        if atomtype_layer1gr(i,j)==2
            atomtype_layer1gr(i,j)=0;
        else
        end
        if atomtype_layer1or(i,j)==2
            atomtype_layer1or(i,j)=0;
        else
        end
        if atomtype_layer2p(i,j)==2
           atomtype_layer2p(i,j)=0;
        else
        end
    end
end

%blue in the 2nd layer is identical to orange in 1st layer
atomtype_layer2blue=atomtype_layer1or;
%red in the 3nd layer is identical to green in 1st layer
atomtype_layer3red=atomtype_layer1gr;
%black in the 3nd layer is identical to purple in 2nd layer
atomtype_layer3black=atomtype_layer2p;

%%filtr of the totally atomtype

for i=1:new_natomW
    for j=1:new_natomL
        for h=1:natomH
        if atomtype_(i,j,h)==2
            atomtype_s(i,j,h)=0;
            ato_s(i,j,h)=1;
        else
         atomtype_s(i,j,h)=1;
         ato_s(i,j,h)=0;
        end
        end
    end
end

natomW =new_natomW;
natomL=new_natomL;
clear startL start2L start3L  new_atomtype_layp  new_atomtype_layor new_atomtype_laygr new_atomtype_
clear new_natomW new_natomL
clear tmp ctW ctL ctz
if load_fixed_atom_distrib
    clear atomtype_
    load('tk.mat');
elseif save_fixed_atom_distrib%save debug data
    save('atomtypee.mat','atomtype_');%change this to save(ddebugfilename);
    error('distribution mat file has been saved, run the program again by setting load_fixed_atom_distrib=1')
end