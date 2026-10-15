%clear;clf;
%load("3e12.mat");
x_A=zeros(size(t));
y_A=zeros(size(t));
z_A=zeros(size(t));
x_B=zeros(size(t));
y_B=zeros(size(t));
z_B=zeros(size(t));
loop=size(t);
for i=1:loop(2)
    x_A(1,i)=mmx(8,4,2,i);
    y_A(1,i)=mmy(8,4,2,i);
    z_A(1,i)=mmz(8,4,2,i);
    x_B(1,i)=mmx(8,4,3,i);
    y_B(1,i)=mmy(8,4,3,i);
    z_B(1,i)=mmz(8,4,3,i);
end
    %3-D figure of a pair of sublattices
    quiver3(0,0,0,x_B(end),y_B(end),z_B(end),'--k','LineWidth',2.5,'Color',[1 0 0]);hold on;
    quiver3(0,0,0,x_A(end),y_A(end),z_A(end),'--k','LineWidth',2.5,'Color',[1 1 0],'MaxHeadSize',0.2);
    quiver3(0,0,0,x_B(1),y_B(1),z_B(1),'LineWidth',2.5,'Color',[0 1 0],'MaxHeadSize',0.2);
    quiver3(0,0,0,x_A(1),y_A(1),z_A(1),'LineWidth',2.5,'Color',[0 1 1],'MaxHeadSize',0.2);
    legend('M_B final','M_A final','M_B init','M_A init');
    if BDSOTRE
        quiver3(0,0,0,psjSHE(1),psjSHE(2),psjSHE(3),'LineWidth',2.5,'Color',[0 0 0],'MaxHeadSize',0.2);
        legend('σ');
    end
    xlabel('mx');ylabel('my');zlabel('mz');
    axis([-1 1 -1 1 -1 1]);
    view(70,30);
    %M-T curve
    figure();
    plot(t*1e9,x_B,t*1e9,y_B,t*1e9,z_B,'linewidth',2);
    xlabel('time(ns)');ylabel('m')
    legend('mx','my','mz');title("m-t")
    axis([0 0.2 -1 1]);

    %2-D visulization of all layers' mmy
    m_y=mmy(:,:,:,end);
    figure();
    min_val = min(m_y(:));  
    max_val = max(m_y(:));  
    for i = 1:size(mmy, 3)
        subplot(2,3, i); 
        imagesc(m_y(:,:,i), [min_val, max_val]); 
        colorbar;  
        title(['Layer ', num2str(i)]); 
        axis on;  
    end