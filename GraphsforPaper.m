%Carley Cook
% This code was written to produce graphs for the Wnt-10b and bone volume
% project
%% Data from Bennett 2005 and Bennett 2007
clear variables
close all
format long e
xdata=[-1, 5, 50]; %Wnt10b Fold Change
ydata=[-29.7,  69.2, 339]; % Bennet Data normalized BV/TV 

kM=25; %kM for MM with k2 Try 50 or 25
N=1;%[6,6,12]; %Number of Cycles
cyclelength=100; %length of cycles
tlag=14; %DDE lag from osteoblast maturation to activation
savegraphs=1; %1 for automatically saving graphs
             %2 for manual saving of graphs
BVandCells=2;%1 to produce Bone Volume and cells for each graph
Estimationcase=1; %1 to produce estimation cases on same graph
Validationcase=2; %1 to produce validation case
OCWnt=1; %1 to produce cells vs wnt
Ratio=2;
barg=2; %1 to produce a bar graph of Wnt fold changes
%% Fitted parameters
% % from bioRxiv version 1
% kg(1)=4.96e-01; %beta1adj
% 
% 
% kg(2)=6.47e-01;%alpha3adj
% 
%              
% kg(3)=2.30e-3;%beta2adj
% 
% kg(4)=2.76e+01;%K

% fitting ydata(3) to N(3) = 12;
    kg(1)=2.15e-01;%beta1adj


    kg(2)=2.98e-02;%6.47e-01;%alpha3adj = k(1) + k(2)

             
    kg(3)=1.28e-03;%beta2adj

    kg(4)=8.42e+00;%K

k=kg;

%% Initial conditions 
S0=180; % Initial Osteocytes 
P0=0; %Initial Pre-Osteoblasts
B0=0; %Initial Osteoblasts
C0=0; %Initial Osteoclasts
z0=100; %Initial Bone Volume
y0=[S0,P0,B0,C0,z0]; % Initial conditions in vector

%% Graphing of Cells and Bone Volume vs time for all cases
if BVandCells==1
    Gwntdose= [-1 0 1.8 1.8 5 50]; %Dose of wnt that will be graphed
    for i=1:length(Gwntdose)
        l='FoldChangeCells';
        l2='FoldChangeBone';
        N=6;
        xaxis=[-5, 605];
        if i==4
            N=12;
            xaxis=[-10, 1210];
        elseif i==6
            N=12;
            xaxis=[-10, 1210];
        end


    [tcalcpt,ycalcpt]=Cyclefunction(k,Gwntdose(i),y0,N,...
        cyclelength);

    %Cells
        figure(i)
        r1=tiledlayout(2,2,'TileSpacing','Compact');
        %Osteocytes
        nexttile
        plot(tcalcpt,ycalcpt(:,1),'r-');
        %xlabel('time(days)','FontSize',12) 
        ylabel('Osteocyte Cells','FontSize',12)
        title('A')
        xlim(xaxis)
        %Pre-osteoblasts
        nexttile
        plot(tcalcpt,ycalcpt(:,2),'b-');
        axis([0 N*100 0 200])
        %xlabel('time(days)','FontSize',12) 
        ylabel('Pre-osteoblast Cells','FontSize',12)
        xlim(xaxis)
        title('B')
        %Osteoblasts
        nexttile
        plot(tcalcpt,ycalcpt(:,3),'g-');
        %xlabel('time(days)','FontSize',12) 
        ylabel('Osteoblast Cells','FontSize',12)
        xlim(xaxis)
        title('C')
        %Osteoclasts
        nexttile
        plot(tcalcpt,ycalcpt(:,4),'m-');
        axis([0 N*100 0 15])
        title('D')
        %xlabel('time(days)','FontSize',12) 
        ylabel('Osteoclast Cells','FontSize',12)
        xlim(xaxis)
        xlabel(r1,'time(days)','FontSize',12)
    %Bone Volume
        figure(i+6)
        plot(tcalcpt,ycalcpt(:,5),'g-','Linewidth',2)
        xlabel('time(days)','FontSize',12) 
        ylabel('Relative bone volume (%)','FontSize',12)
        xlim(xaxis)

    if savegraphs==1
    q=string(Gwntdose(i));
    v=strcat(l,q);
    v2=strcat(l2,q);
    if i== 3
        q='18';
        v=strcat(l,q);
        v2=strcat(l2,q);
    end
    if i== 4
        q='18';
        w='12months';
        v=strcat(l,q,w);
        v2=strcat(l2,q,w);
    end
    saveas(figure(i),v,'png')
    saveas(figure(i+6),v2,'png')
    end 

    end
end
%% Graphing of cases used for parameter estimation
if Estimationcase==1
    N=[6,6,12]; %Number of Cycles
    [tcalcpt1,ycalcpt1]=Cyclefunction(k,-1,y0,N(1),cyclelength);
    [tcalcpt5,ycalcpt5]=Cyclefunction(k,5,y0,N(2),cyclelength);
    [tcalcpt50,ycalcpt50]=Cyclefunction(k,50,y0,N(3),cyclelength);
        
        %Tile the graphs for the paper
        figure(13)
        set(gcf,'Position',[100 100 985 420])
        r=tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
        nexttile

        plot(tcalcpt1,ycalcpt1(:,5),'g-','Linewidth',2)
        hold on
        plot(tcalcpt5,ycalcpt5(:,5),'b:','Linewidth',2)
        plot(tcalcpt50,ycalcpt50(:,5),'r-.','Linewidth',2)
         xdata2007=[600,600];
         ydata2007=[-29.7,  69.2]+100;
         xdata2005=1200;
         ydata2005=339+100;
        plot(xdata2007,ydata2007,'ko','Linewidth',2)
        plot(xdata2005,ydata2005,'bs','Linewidth',2)
        %plot(600,ydata(1:2)+100,'ko','Linewidth',2)
        %plot(600,ydata(2)+100,'ko','Linewidth',2)
        %plot(1200,ydata(3)+100,'bs','Linewidth',2)
        legend("-1 Fold","5 Fold","50 Fold","Bennett 2007 Data","Bennett 2005 Data",...
            'Location','Best','FontSize',12)
        xlabel('time(days)','FontSize',12) 
        ylabel('Relative bone volume (%)','FontSize',12)
        xlim([-10 1210])
        title('A')
        nexttile
        %Reproducing figure 1 from BoneCompartmentUp
         xp = linspace(-1,50,100);
         Np6=zeros(length(xp),1);
         Np12=zeros(length(xp),1);
         Np6(:)=6;
         Np12(:)=12;
         ycalcp6 = Graham2013(k,xp,y0,Np6,cyclelength);
         ycalcp12 = Graham2013(k,xp,y0,Np12,cyclelength);
         plot(xp,ycalcp6(1,:),'r','Linewidth',2);
         hold on
         plot(xp,ycalcp12(1,:),'b:','Linewidth',2);
         xdata2007=[-1,5];
         ydata2007=[-29.7,  69.2];
         xdata2005=50;
         ydata2005=339;
         plot(xdata2007,ydata2007,'ko','Linewidth',2)
         plot(xdata2005,ydata2005,'bs','Linewidth',2)
         legend("Simulation Results 6 Cycles","Simulation Results 12 Cycles","Bennett 2007 Data","Bennett 2005 Data",'Location','Best','FontSize',12)
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('BV//TV (% Change from normal Wnt-10b)','FontSize',12)
        title('B')
    if savegraphs==1
    saveas(figure(13),'EstimationResults','png')
    end
end

%% Graphing of model validation
if Validationcase==1
    [tcalcpt12,ycalcpt12]=Cyclefunction(k,1.8,y0,12,cyclelength);
    [tcalcpt12top,ycalcpt12top]=Cyclefunction(k,2.4,y0,12,cyclelength);
    [tcalcpt12bottom,ycalcpt12bottom]=Cyclefunction(k,1.2,y0,12,cyclelength);
    figure(14)
        plot(tcalcpt12,ycalcpt12(:,5),'Linewidth',2)
        hold on
        errorbar( 600 , 126.6 , 15,'o','Linewidth',2 )
        errorbar( 1200 , 136.6 , 29,'s','Linewidth',2 )
        shade(tcalcpt12top,ycalcpt12top(:,5),tcalcpt12bottom,ycalcpt12bottom(:,5),'FillType',[1 2;2 1],'Color', [.8 .8 .8]);
        legend("Simulation results","Data corresponding to 6 cycles",...
            "Data corresponding to 12 cycles",'Location',...
            'Best','FontSize',12)
        xlabel('Time (days)','FontSize',12) 
        ylabel('Relative bone volume (%)','FontSize',12)
        xlim([-10 1210])
    if savegraphs==1
    saveas(figure(14),'ValidationResults','png')
    end
end

%% Graphing important cell population vs Wnt
if OCWnt==1
    N=1;
    Gwntdose= [-1 0 1.8 5 25 50];
    
    [tcalcpt1,ycalcpt1]=Cyclefunction(k,-1,y0,N,...
        cyclelength);
    [tcalcpt0,ycalcpt0]=Cyclefunction(k,0,y0,N,...
        cyclelength);
    [tcalcpt18,ycalcpt18]=Cyclefunction(k,1.8,y0,N,...
        cyclelength);
    [tcalcpt5,ycalcpt5]=Cyclefunction(k,5,y0,N,...
        cyclelength);
    [tcalcpt25,ycalcpt25]=Cyclefunction(k,25,y0,N,...
        cyclelength);
    [tcalcpt50,ycalcpt50]=Cyclefunction(k,50,y0,N,...
        cyclelength);
    
%     %Area under pre-ob curve 
%     APO(1)=trapz(tcalcpt1,ycalcpt1(:,2));
%     APO(2)=trapz(tcalcpt0,ycalcpt0(:,2));
%     APO(3)=trapz(tcalcpt18,ycalcpt18(:,2));
%     APO(4)=trapz(tcalcpt5,ycalcpt5(:,2));
%     APO(5)=trapz(tcalcpt25,ycalcpt25(:,2));
%     APO(6)=trapz(tcalcpt50,ycalcpt50(:,2));
%     %Area under osteoblast curve
%     AO(1)=trapz(tcalcpt1,ycalcpt1(:,3));
%     AO(2)=trapz(tcalcpt0,ycalcpt0(:,3));
%     AO(3)=trapz(tcalcpt18,ycalcpt18(:,3));
%     AO(4)=trapz(tcalcpt5,ycalcpt5(:,3));
%     AO(5)=trapz(tcalcpt25,ycalcpt25(:,3));
%     AO(6)=trapz(tcalcpt50,ycalcpt50(:,3));
%     %Area under osteoclast curve
%     AC(1)=trapz(tcalcpt1,ycalcpt1(:,4));
%     AC(2)=trapz(tcalcpt0,ycalcpt0(:,4));
%     AC(3)=trapz(tcalcpt18,ycalcpt18(:,4));
%     AC(4)=trapz(tcalcpt5,ycalcpt5(:,4));
%     AC(5)=trapz(tcalcpt25,ycalcpt25(:,4));
%     AC(6)=trapz(tcalcpt50,ycalcpt50(:,4));
    
    
    figure(15)
        oc(1)=max(ycalcpt1(:,4));
        oc(2)=max(ycalcpt0(:,4));
        oc(3)=max(ycalcpt18(:,4));
        oc(4)=max(ycalcpt5(:,4));
        oc(5)=max(ycalcpt25(:,4));
        oc(6)=max(ycalcpt50(:,4));
        plot(Gwntdose,oc,'m-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Osteoclast Cells','FontSize',12)
        xlim([-5 55])
    figure(16)
        PO(1)=max(ycalcpt1(:,2));
        PO(2)=max(ycalcpt0(:,2));
        PO(3)=max(ycalcpt18(:,2));
        PO(4)=max(ycalcpt5(:,2));
        PO(5)=max(ycalcpt25(:,2));
        PO(6)=max(ycalcpt50(:,2));
        plot(Gwntdose,PO,'b-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Pre-osteoblast Cells','FontSize',12)
        xlim([-5 55])
    figure(17)
        O(1)=max(ycalcpt1(:,3));
        O(2)=max(ycalcpt0(:,3));
        O(3)=max(ycalcpt18(:,3));
        O(4)=max(ycalcpt5(:,3));
        O(5)=max(ycalcpt25(:,3));
        O(6)=max(ycalcpt50(:,3));
        plot(Gwntdose,O,'b-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Osteoblast Cells','FontSize',12)
        xlim([-5 55])    
    
    figure(18)
    r2=tiledlayout(2,2,'TileSpacing','Compact');
    
        nexttile
        plot(tcalcpt1,ycalcpt1(:,1),'Linewidth',2);
        hold on
        plot(tcalcpt0,ycalcpt0(:,1),'Linewidth',2);
        plot(tcalcpt18,ycalcpt18(:,1),'Linewidth',2);
        plot(tcalcpt5,ycalcpt5(:,1),'Linewidth',2);
        plot(tcalcpt25,ycalcpt25(:,1),'Linewidth',2);
        plot(tcalcpt50,ycalcpt50(:,1),'Linewidth',2);
        %xlabel('time(days)','FontSize',12) 
        ylabel('Osteocytes (cells)','FontSize',12)
        legend('-1', '0', '1.8', '5', '25', '50')
        xlim([-5 105])
        title('A')
        
        nexttile
    
        plot(tcalcpt1,ycalcpt1(:,2),'Linewidth',2);
        hold on
        plot(tcalcpt0,ycalcpt0(:,2),'Linewidth',2);
        plot(tcalcpt18,ycalcpt18(:,2),'Linewidth',2);
        plot(tcalcpt5,ycalcpt5(:,2),'Linewidth',2);
        plot(tcalcpt25,ycalcpt25(:,2),'Linewidth',2);
        plot(tcalcpt50,ycalcpt50(:,2),'Linewidth',2);
        %xlabel('time(days)','FontSize',12) 
        ylabel('Pre-Osteoblasts (cells)','FontSize',12)
        legend('-1', '0', '1.8', '5', '25', '50')
        xlim([-5 105])
        title('B')
        nexttile
    
        plot(tcalcpt1,ycalcpt1(:,3),'Linewidth',2);
        hold on
        plot(tcalcpt0,ycalcpt0(:,3),'Linewidth',2);
        plot(tcalcpt18,ycalcpt18(:,3),'Linewidth',2);
        plot(tcalcpt5,ycalcpt5(:,3),'Linewidth',2);
        plot(tcalcpt25,ycalcpt25(:,3),'Linewidth',2);
        plot(tcalcpt50,ycalcpt50(:,3),'Linewidth',2);
        %xlabel('time(days)','FontSize',12) 
        ylabel('Osteoblast (cells)','FontSize',12)
        legend('-1', '0', '1.8', '5', '25', '50')
        xlim([-5 105])
        title('C')
    
        nexttile
        plot(tcalcpt1,ycalcpt1(:,4),'Linewidth',2);
        hold on
        plot(tcalcpt0,ycalcpt0(:,4),'Linewidth',2);
        plot(tcalcpt18,ycalcpt18(:,4),'Linewidth',2);
        plot(tcalcpt5,ycalcpt5(:,4),'Linewidth',2);
        plot(tcalcpt25,ycalcpt25(:,4),'Linewidth',2);
        plot(tcalcpt50,ycalcpt50(:,4),'Linewidth',2);
        %xlabel('time(days)','FontSize',12) 
        ylabel('Osteoclast (cells)','FontSize',12)
        legend('-1', '0', '1.8', '5', '25', '50')
        xlim([-5 105])
        title('D')
        xlabel(r2,'time(days)','FontSize',12)
    
    if savegraphs==1
    saveas(figure(15),'OCWnt','png')
    saveas(figure(16),'POWnt','png')
    saveas(figure(17),'OWnt','png')
    saveas(figure(18),'Celldynam','png')
    end
end

%% Getting ratios for area under the curve
if Ratio==1
    N=1;
    Gwntdose=linspace(-1,50,100);
    APO = zeros(size(Gwntdose));
    AO = zeros(size(Gwntdose));
    AC = zeros(size(Gwntdose));
    PO = zeros(size(Gwntdose));
    O = zeros(size(Gwntdose));
    % fill in other relevant quantities as needed
    for i = 1:length(Gwntdose)
        [tcalcpt,ycalcpt]=Cyclefunction(k,Gwntdose(i),y0,N,...
        cyclelength);
        %Area under pre-ob curve 
        APO(i)=trapz(tcalcpt,ycalcpt(:,2));
        %Area under osteoblast curve
        AO(i)=trapz(tcalcpt,ycalcpt(:,3));
        %Area under osteoclast curve
        AC(i)=trapz(tcalcpt,ycalcpt(:,4));
        %Preosteoblast max
        PO(i)=max(ycalcpt(:,2));
        %Osteoblast max
        O(i)=max(ycalcpt(:,3));
    end
    
    figure(19)
    plot(Gwntdose,PO./O,'Linewidth',2);
    xlabel('Wnt-10b (Fold Change)','FontSize',12) 
    ylabel('Pre-Osteoblasts/Osteoblast  Max Cells','FontSize',12)
    
    figure(20)
    r3=tiledlayout(1,2,'TileSpacing','Compact');
    %figure(23)
        nexttile
        plot(Gwntdose,APO./AO,'Linewidth',2);
        %xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Pre-Osteoblasts/Osteoblast AUC','FontSize',12)
        xlim([-5 55])
        title('A')
    %figure(24)
        nexttile
        plot(Gwntdose,AC./AO,'Linewidth',2);
        %xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Osteoclast/Osteoblast AUC','FontSize',12)
        xlim([-5 55])
        title('B')
        xlabel(r3,'Wnt-10b (Fold Change)','FontSize',12)
        
    if savegraphs==1
    saveas(figure(19),'POOBmax','png')
    saveas(figure(20),'AUC','png')
    end
        
end

%% 
if barg==1
    figure(21)
        bary=[-1,0,1.8,5,50];
        barx=categorical({'No Wnt-10b','Normal Wnt-10b',...
            'Wnt-10b increase 1','Wnt-10b increase 2',...
            'Wnt-10b increase 3'});
        bar(barx,bary)
        ylabel('Wnt-10b (Fold Change)','FontSize',12)
    if savegraphs==1
        saveas(figure(21),'BarFold','png')
    end
end



%% Functions needed to produce graphs
%% Define ODE equations with variable parameters
function dydt= ODEeq(t,y,k,x)%x is a scalar wnt10b dose

Bone=1;
alpha_1=0.5;
alpha_2=0.1;
alpha_3=0.1;
beta_1=0.1;
delta=0.1;
beta_2=0.1;
alpha_4=0.1;
K_S=200;
k1=.698331;%Graham2013 paper has .7 to get the ss to 100 .6983 works for ode
         %When dde .69825
         %When dde with round(c).71575
         %When dde with round(c) in dcdt .72249
k2=0.015445;
g_31=1;
g_21=2;
g_22=1;
g_32=1;
g_41=1;
g_42=1;
g_43=-1;
g_44=1;
f_12=1;
f_14=1;
f_23=1;
f_34=1;
epsilon=1;
beta_3=0.1;
rho=20;

    S =y(1);
    P =y(2);
    B =y(3);
    C =y(4);
    z=y(5);
    %ylag=Z(:,1);
    %Blag=ylag(3);


    %Set parameter definitions
    beta1adj = k(1);
    alpha3adj=k(2)+k(1);
    beta2adj = k(3);
    K=k(4);
    %K=25;
    piwnta=(x/(x+K));
    
    %Algebraic equations needed for the ODEs
    Differentiation_of_Osteoblast_to_Osteocytes = Bone*(alpha_1)*power(B,g_31)*max((1-S/K_S),0);
    
    %Differentiation_of_MSC_cells_to_PreOsteoblast_cells = Bone*(alpha_2)*power(S,g_21)*max((1-S/K_S),0)^g_22*(1+(power(P,f_12)*beta1adj*x));
    %Differentiation_of_MSC_cells_to_PreOsteoblast_cells = Bone*(alpha_2)*power(S,g_21)*max((1-S/K_S),0)^g_22*(1+(alphanew*x));
    Differentiation_of_MSC_cells_to_PreOsteoblast_cells = Bone*(alpha_2)*power(S,g_21)*max((1-S/K_S),0)^g_22;
    
    %Proliferation_of_preosteoblasts = Bone*alpha_3*power(P,g_32)*max((1-S/K_S),0)*(1+(power(P,f_12)*beta1adj*x));
    %Proliferation_of_preosteoblasts = Bone*alpha_3*power(P,g_32)*max((1-S/K_S),0)*(1+alphanew*x);
    Proliferation_of_preosteoblasts = Bone*((alpha_3)*power(P,g_32)*max((1-S/K_S),0)+((alpha3adj*piwnta)*power(P,f_12)));
    
    %Differentiation_of_PreOsteoblast_to_mature_osteoblast = Bone*(beta_1+(beta1adj*x))*power(P,f_12)*power(C,f_14);
    Differentiation_of_PreOsteoblast_to_mature_osteoblast = Bone*((beta_1*power(P,f_12)*power(C,f_14))+((beta1adj*piwnta)*power(P,f_12)));
    
    %Apoptosis_of_preosteoblast = Bone*(delta-beta2adj*C*x)*P;
    Apoptosis_of_preosteoblast = Bone*(delta)*P;
    
    %Apoptosis_of_osteoblasts = Bone*(beta_2-((beta2adj*piwntr)))*power(B,f_23);
    Apoptosis_of_osteoblasts = Bone*((beta_2-beta2adj*piwnta)*power(B,f_23));
    
    %Differentiation_of_preosteoclast_to_osteoclasts = Bone*alpha_4*power...
    %   (S,g_41)*power(P,g_42)*min((epsilon+B),(epsilon+494.3))^g_43*...
    %   max((1-S/K_S),0)^g_44;
    Differentiation_of_preosteoclast_to_osteoclasts = Bone*alpha_4*power...
        (S,g_41)*power(P,g_42)*(epsilon+B)^g_43*...
        max((1-S/K_S),0)^g_44;
    
    Apoptosis_of_osteoclasts = Bone*beta_3*power(C,f_34);
    %Apoptosis_of_osteoclasts = Bone*beta_3*power(round(C),f_34);
    
    Resorption_of_bone = Bone*k1*C;
    
    %Formation_of_bone = Bone*((knew*Blag));
    %Formation_of_bone = Bone*((knew*B));
    %Formation_of_bone = Bone*((k2*Blag));
    Formation_of_bone = Bone*((k2*B));
    %ODEs

    %d([Osteocytes (S)])/dt = 1/Bone*Differentiation_of_Osteoblast_to_Osteocytes;
    dydt(1)=1/Bone*Differentiation_of_Osteoblast_to_Osteocytes;

    %d([Pre-Osteoblasts (P)])/dt = 1/Bone*(Differentiation_of_MSC_cells_to_PreOsteoblast_cells + Proliferation_of_preosteoblasts - Differentiation_of_PreOsteoblast_to_mature_osteoblast - Apoptosis_of_preosteoblast)
    dydt(2)=1/Bone*...
        (Differentiation_of_MSC_cells_to_PreOsteoblast_cells ...
        + Proliferation_of_preosteoblasts - ...
        Differentiation_of_PreOsteoblast_to_mature_osteoblast...
        - Apoptosis_of_preosteoblast);

    %d([Osteoblasts (B)])/dt = 1/Bone*
    %(-Differentiation_of_Osteoblast_to_Osteocytes + 
    %Differentiation_of_PreOsteoblast_to_mature_osteoblast
    %- Apoptosis_of_osteoblasts)
    dydt(3)=1/Bone*(-Differentiation_of_Osteoblast_to_Osteocytes...
        + Differentiation_of_PreOsteoblast_to_mature_osteoblast...
        - Apoptosis_of_osteoblasts);

    %d([Osteoclasts (C)])/dt = 1/Bone*(Differentiation_of_preosteoclast_to_osteoclasts - Apoptosis_of_osteoclasts)
    dydt(4)=1/Bone*(Differentiation_of_preosteoclast_to_osteoclasts...
        - Apoptosis_of_osteoclasts);

    %d([Bone volume (z)])/dt = 1/Bone*(-Resorption_of_bone + Formation_of_bone)
    dydt(5)=1/Bone*(-Resorption_of_bone + Formation_of_bone);    
    

    dydt=[dydt(1)
          dydt(2)
          dydt(3)
          dydt(4)
          dydt(5)];
    
end

%% Solve ODE using variable parameters
function yout = Graham2013(k,xdata,y0,N,cyclelength)
 
    for i = 1:length(xdata)
        
        [~,ycalc] = Cyclefunction(k,xdata(i),y0,N,cyclelength);
        
        yBV(i,1)=ycalc(end,5);
        yout(i,1)=yBV(i,:)-100;
        
            
    end
    yout = yout';
end
%% Cycle Function
     function [combined_tcalc_N_cycles, combined_ycalc_N_cycles] = Cyclefunction(k,xdata,y0,N,cyclelength)
     startindex=1; %Indices are used to combine the loops into a single column
     finalindex=1;
     %oldsol=[];
          for j= 1:N

            %tspan = (j-1)*cyclelength:0.01:j*cyclelength;
            %tspan =linspace((j-1)*cyclelength,j*cyclelength,101);
            tspan = linspace((j-1)*cyclelength,j*cyclelength,551);
            %tspan = [(j-1)*cyclelength,j*cyclelength];
             %[tcalc,ycalc] = ode23s(@(t,y) ODEeq(t,y,k,xdata(i)),tspan,y0);
            opts=odeset('NonNegative',(1:5),'reltol',1e-7);
            %opts=odeset('NonNegative',(1:5));
            [tcalc,ycalc] = ode15s(@(t,y) ODEeq(t,y,k,xdata),tspan,y0,opts);

             ycalc(:,1);
             %Reset for next time interval
             y0=ycalc(end,:);
             idx=(y0<1);
             y0(idx)=0; %Sets fractions of cells to zero
             y0(1,1)=y0(1,1)-20; %Reduces osteocytes to initate next cycle
             %oldsol.y(:,end)=y0';
             finalindex=finalindex+length(tcalc)-1;
             combined_tcalc_N_cycles(startindex:finalindex,1)=tcalc;
             combined_ycalc_N_cycles(startindex:finalindex,:)=ycalc;
             startindex=startindex+length(tcalc)-1; %equals previous final index
          end

     end