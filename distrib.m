%% distribute atoms
% 1: have atom, 0: no atom
%according to the atom distribution, we divide the syatem with three
%different layers.The honeycomb structure is converted into a 9*9 matrix
%       *---*         [1 0 1]
%      /     \
%     *   *   *       [1 1 1]
%      \     /
%       *---*         [1 0 1]
%for layer 1
%       g*---o*       [1 0 1]
%      /      \
%     o*       g*     [1 0 1]
%      \      /
%       g*---o*       [1 0 1]
%layer 1 is divied into green (g) and orange (o) atoms to calculate
%for layer 2
%        ---b*       [0 0 1]
%      /      \
%     b*  p*         [1 1 0]
%      \      /
%        ---b*       [0 0 1]
%layer 2 is divied into blue (b) and purple (p) atoms to calculate
%for layer 3
%       r*---        [1 0 0]
%      /      \
%         b*  r*     [0 1 1]
%      \      /
%       r*---        [1 0 0]
%layer 3 is divied into black (b) and red (r) atoms to calculate
atomtype_ = zeros(natomW, natomL, natomH);
%% determine if there is atom in the first layer
for ctW=1:natomW
    for ctL = 1:natomL
        for ctH = 1:natomH
            if (mod(ctH-1,3)==0)
                if mod(ctL, 2) == 1  % odd ctL
                    atomtype_(ctW, ctL, ctH) = 1;
                end
            end
        end
    end
end
%% determine if there is atom in the second layer
for ctW = 1:natomW
    for ctH = 1:natomH
        if (mod(ctH-1,3)==1)
            if mod(ctW, 2) == 1  % odd ctW
                startL = 3;
                start2L = 4;
            else  % even ctW
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
%% determine if there is atom in the third layer
for ctW = 1:natomW
    for ctH = 1:natomH
        if (mod(ctH-1,3)==2)
            if mod(ctW, 2) == 1  % odd ctW
                startL = 1;
                start2L = 4;
            else   % even ctW
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
%% pad two column and row for the boundary
new_natomW = natomW + 4;
new_natomL = natomL + 4;
new_atomtype_ = zeros(new_natomW, new_natomL, natomH);
for ctH = 1:natomH
    new_atomtype_(3:end-2, 3:end-2, ctH) = atomtype_(:, :, ctH);
end
%% filter for each layer
atomtype_layer1gr = zeros(natomW,natomL,'gpuArray');
atomtype_layer1or = zeros(natomW,natomL,'gpuArray');
atomtype_layer2p = zeros(natomW,natomL,'gpuArray');
for ctW = 1:natomW
    if mod(ctW, 2) == 1  % odd ctW
        startL = 1;
        start2L = 3;
        start3L = 4;
    else  % even ctW
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
new_atomtype_laygr = zeros(new_natomW,new_natomL,'gpuArray');
new_atomtype_layor = zeros(new_natomW,new_natomL,'gpuArray');
new_atomtype_layp = zeros(new_natomW,new_natomL,'gpuArray');
new_atomtype_laygr (3:end-2, 3:end-2) = atomtype_layer1gr(:, :);
new_atomtype_layor (3:end-2, 3:end-2) = atomtype_layer1or(:, :);
new_atomtype_layp (3:end-2, 3:end-2) = atomtype_layer2p(:, :);

atomtype_=new_atomtype_;
natomW =new_natomW;
natomL=new_natomL;
atomtype_layer2p =new_atomtype_layp;
atomtype_layer1or =new_atomtype_layor;
atomtype_layer1gr=new_atomtype_laygr;

%blue in the 2nd layer is identical to orange in 1st layer
atomtype_layer2blue=atomtype_layer1or;
%red in the 3nd layer is identical to green in 1st layer
atomtype_layer3red=atomtype_layer1gr;
%black in the 3nd layer is identical to purple in 2nd layer
atomtype_layer3black=atomtype_layer2p;
ato_s= abs(atomtype_-1);
% ato_s is just the opposite, this matrix is used to avoid the NaN result
% in the rk4 calculation. 0=have atom, 1=no atom.
%%filtr of the totally atomtype

clear startL start2L start3L new_atomtype_layp new_atomtype_layor new_atomtype_laygr new_atomtype_
clear new_natomW new_natomL
clear tmp ctW ctL ctz
if load_fixed_atom_distrib
    clear atomtype_
    load('atomtypee.mat','atomtype_');
elseif save_fixed_atom_distrib%save debug data
    save('atomtypee.mat','atomtype_');%change this to save(ddebugfilename);
    error('distribution mat file has been saved, run the program again by setting load_fixed_atom_distrib=1')
end