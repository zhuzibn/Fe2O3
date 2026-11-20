function [mask,colors,A] = fe2o3_expected_geometry(core)
% Independent direct-index transcription of docs/structure-exchange-guide.md
% sections 2, 4, 5. Array bookkeeping only, not verified crystallography.
shape=core+[4,4,0]; mask=zeros(shape); colors=zeros(shape);
starts=[1,3,3,4,1,4;3,1,1,2,3,2];
for k=1:core(3)
    pair=2*mod(k-1,3)+(1:2);
    for r=1:core(1)
        for c=1:core(2)
            for color=pair
                if mod(c-starts(2-mod(r,2),color),4)==0
                    mask(r+2,c+2,k)=1; colors(r+2,c+2,k)=color;
                end
            end
        end
    end
end
O=[0,0]; P=[0,1;1,-1;-1,-1]; M=[0,-1;1,1;-1,1];
P2=[0,2;1,0;-1,0]; M2=[0,-2;1,0;-1,0];
at=@(xy,k) [xy,repmat(k,size(xy,1),1)];
offsets={at(O,-1),at(P2,0),at(M,1),[at(P2,1);at(M,-1)];...
    at(O,1),at(M2,0),at(P,-1),[at(P,1);at(M2,-1)];...
    at(O,-1),at(P,0),at(M2,1),[at(P,1);at(M2,-1)];...
    at(O,1),at(M,0),at(P,-1),[at(P,1);at(M,-1)];...
    at(O,1),at(M,0),at(P2,-1),[at(M,-1);at(P2,1)];...
    at(O,-1),at(P,0),at(M,1),[at(M,-1);at(P,1)]};
sites=find(mask); ids=zeros(shape); ids(sites)=1:numel(sites);
A=repmat({zeros(numel(sites))},1,4);
for n=1:numel(sites)
    [r,c,k]=ind2sub(shape,sites(n));
    for shell=1:4
        for offset=offsets{colors(sites(n)),shell}'
            q=[r;c;k]+offset;
            if any(q<1) || any(q>shape'), continue; end
            id=ids(q(1),q(2),q(3));
            if id, A{shell}(n,id)=A{shell}(n,id)+1; end
        end
    end
end
end
