%Carley Cook
% This code was written to estimate the change in parameters when a change
% in Wnt10b occurs
%% Data from Bennett 2005 and Bennett 2007
clear variables
close all
format long e
xdata=[-1,-1, 5, 50]; %Wnt10b Fold Change
ydata=[-29.7, -41.9, 69.2, 339]; % Bennet Data normalized BV/TV

ParamY=1;% ParamY=1 for running parameter estimation code
         % ParamY=2 for just graphing
N=[4,6,6,12]; %Number of Cycles
cyclelength=100; %length of cycles
Gwntdose=50; %Dose of wnt that will be graphed
OCWnt=2; %1 to look at extra graphs
graphing=2;
ParamG=7; %1 for latin hypercube sampling
         %2 for random number sampling
         %3 for modifided latin hypercube sampling
         %4 for narrowed down parameter space
         %5 for narrowed down parameter space even more
         %6 for narrowed down parameter space for k3
         %6 for narrowed down parameter space using sd and average
%% Initial conditions 
% S0=200 for SS and K_S-rho=(180) for activation
S0=180; % Initial Osteocytes 
P0=0; %Initial Pre-Osteoblasts
B0=0; %Initial Osteoblasts
C0=0; %Initial Osteoclasts
z0=100; %Initial Bone Volume
y0=[S0,P0,B0,C0,z0]; % Initial conditions in vector

%% Parameter Estimation with k parameters, BV, and resnorm outputs
%Turning on parameter estimation
if ParamY==1
% Setting up parameter estimation with multistart

repeats = 100;
Obj_best = zeros(repeats,1);
p_best = zeros(repeats,4);
rng default
LHCparams=lhsdesign(repeats,4);
tic
parpool(8)

parfor p = 1:repeats
    A = []; B = []; A_eq = []; B_eq = []; % no inequality and equality constraints 
    lb=[0,0,0,1];  %lower bound of parameters 
    ub=[1,1,1, 100]; %upper bound of parameters
    if ParamG==1 %Trying different parameter spaces
        kguess=[LHCparams(p,1),LHCparams(p,2),LHCparams(p,3),LHCparams(p,4)*ub(4);]
    
    elseif ParamG==2
        kguess=[rand(1), rand(1), rand(1), rand(1)*100]
    
    elseif ParamG==3
        if p<=100
            kguess=[LHCparams(p,1)*10^-2,LHCparams(p,2),LHCparams(p,3),LHCparams(p,4)*ub(4);]
        elseif p>100 && p<=200
            kguess=[LHCparams(p,1),LHCparams(p,2)*10^-2,LHCparams(p,3),LHCparams(p,4)*ub(4);]
        elseif p>200 && p<=300
            ub=[100,100,100,100]; %upper bound of parameters
            kguess=[LHCparams(p,1)*10,LHCparams(p,2)*10,LHCparams(p,3)*10,LHCparams(p,4)*ub(4);]
        end
    
    elseif ParamG==4
        ub=[1,1,1,10]; %upper bound of parameters
        kguess=[LHCparams(p,1),LHCparams(p,2),LHCparams(p,3)*10^-3,LHCparams(p,4)*ub(4)];
    
    elseif ParamG==5
        lb=[.0001,.0001,0,1];
        ub=[1,1,1,7]; %upper bound of parameters
        kguess=[LHCparams(p,1),LHCparams(p,2),LHCparams(p,3)*10^-2,LHCparams(p,4)*ub(4)];
    elseif ParamG==6
        lb=[.0001,.00001,1*10^-7,1];
        ub=[1,.1,1,10]; %upper bound of parameters
        kguess=[LHCparams(p,1),LHCparams(p,2)*ub(2),LHCparams(p,3),LHCparams(p,4)*ub(4)];
    elseif ParamG==7
        lb=[.0001,.00001,1*10^-7,1];
        %k=[0.189306506149579	0.0635392435460320	0.000776198394479350	6.29750459645756];
        %SD=[0.0152408945053190	0.0154632185194991	0.00100572290174402	0.0995400947720392];
        ub=[1,.1,1,10]; %upper bound of parameters
        kguess=[normrnd(0.189306506149579,0.0152408945053190), normrnd(0.0635392435460320,0.0154632185194991), normrnd(0.000776198394479350,0.00100572290174402),normrnd(6.29750459645756,0.0995400947720392)]
    end
    kguessStore(p,:)=kguess
    %implementing lsqcurvefit
    OPTIONS = optimoptions('lsqcurvefit','StepTolerance',1e-14,...
    'FunctionTolerance',1e-14,'optimalitytolerance', 1e-6,'MaxFunctionEvaluations',10000,'Algorithm','levenberg-marquardt');
    [p_best(p,:),Obj_best(p),residual] = lsqcurvefit(@(k,xdata) Graham2013(k,xdata,y0,...
    N,cyclelength),kguess,xdata,ydata, lb, ub, OPTIONS);
    p
end
toc
delete(gcp)
[global_Obj_best,index_best] = min(Obj_best);
k = p_best(index_best,:);

end      
%% Graphing of BV vs Wnt
if graphing==1
figure(1)
 %xp = linspace(xdata(1),xdata(end),1001);
 %xp = linspace(xdata(1),xdata(end),100);
 xp = linspace(-1,50,100);
 Np4=zeros(length(xp),1);
 Np6=zeros(length(xp),1);
 Np12=zeros(length(xp),1);
 Np4(:)=4;
 Np6(:)=6;
 Np12(:)=12;
 ycalcp4 = Graham2013(k,xp,y0,Np4,cyclelength);
 ycalcp6 = Graham2013(k,xp,y0,Np6,cyclelength);
 ycalcp12 = Graham2013(k,xp,y0,Np12,cyclelength);
 plot(xp,ycalcp4(1,:),'r','Linewidth',2);
 hold on
 plot(xp,ycalcp6(1,:),'g','Linewidth',2);
 plot(xp,ycalcp12(1,:),'b','Linewidth',2);
 xdata2007=[-1,5];
 ydata2007=[-41.9,  69.2];
 xdata2005=[-1,50];
 ydata2005=[-29.7,339];
 plot(xdata2007,ydata2007,'ko','Linewidth',2)
 plot(xdata2005,ydata2005,'bs','Linewidth',2)
% Checking that the rounded paramters produce a simular simulation result
%  if ParamY==2
%  ycalcr6 = Graham2013(kr,xp,y0,Np6,cyclelength);
%  plot(xp,ycalcr6(1,:),'b','Linewidth',2)
%  ycalcr12 = Graham2013(kr,xp,y0,Np12,cyclelength);
%  plot(xp,ycalcr12(1,:),'b','Linewidth',2)
%  end
 legend("Simulation Results 4 Cycles","Simulation Results 6 Cycles","Simulation Results 12 Cycles","Bennett 2007 Data","Bennett 2005 Data",'Location','Best','FontSize',12)
xlabel('Wnt-10b (Fold Change)','FontSize',15) 
ylabel('BV//TV (% Change from normal Wnt-10b)','FontSize',15)


%Plot all cases on the same graph

%% Graphing of Cells and Bone Volume vs time
figure(2)
[tcalcpt,ycalcpt]=Cyclefunction(k,Gwntdose,y0,N,cyclelength);

tiledlayout(2,2)
%Osteocytes
nexttile
plot(tcalcpt,ycalcpt(:,1),'r-');
xlabel('time(days)') 
ylabel('Osteocyte Cells')

%Pre-osteoblasts
nexttile
plot(tcalcpt,ycalcpt(:,2),'b-');
xlabel('time(days)') 
ylabel('Pre-osteoblast Cells')

%Osteoblasts
nexttile
plot(tcalcpt,ycalcpt(:,3),'g-');
xlabel('time(days)') 
ylabel('Osteoblast Cells')

%Osteoclasts
nexttile
plot(tcalcpt,ycalcpt(:,4),'m-');
xlabel('time(days)') 
ylabel('Osteoclast Cells')

figure(3)
plot(tcalcpt,ycalcpt(:,5),'g-','Linewidth',2)
xlabel('time(days)','FontSize',15) 
ylabel('Relative bone volume (%)','FontSize',15)


    [tcalcpt1,ycalcpt1]=Cyclefunction(k,-1,y0,N(2),cyclelength);
%     [tcalcpt18,ycalcpt18]=Cyclefunction(k,1.8,y0,12,cyclelength);
    [tcalcpt5,ycalcpt5]=Cyclefunction(k,5,y0,N(3),cyclelength);
    [tcalcpt50,ycalcpt50]=Cyclefunction(k,50,y0,N(4),cyclelength);
    figure(14)
    plot(tcalcpt1,ycalcpt1(:,5),'g-','Linewidth',2)
    hold on
    %         plot(tcalcpt18,ycalcpt18(:,5),'b:','Linewidth',2)
    plot(tcalcpt5,ycalcpt5(:,5),'b:','Linewidth',2)
    plot(tcalcpt50,ycalcpt50(:,5),'r-.','Linewidth',2)
    plot(N(1)*cyclelength,ydata(1)+100,'ko','Linewidth',2)
    plot(N(2)*cyclelength,ydata(2)+100,'ko','Linewidth',2)
    plot(N(3)*cyclelength,ydata(3)+100,'ko','Linewidth',2)
    plot(N(4)*cyclelength,ydata(4)+100,'ko','Linewidth',2)
    legend("-1 Fold","5 Fold","50 Fold","Literature Data",...
        'Location','Best','FontSize',12)
    xlabel('time(days)','FontSize',15)
    ylabel('Relative bone volume (%)','FontSize',15)
if OCWnt==1
    N=1;
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
    Gwntdose= [-1 0 1.8 5 25 50];
    %Area under pre-ob curve 
    APO(1)=trapz(tcalcpt1,ycalcpt1(:,2));
    APO(2)=trapz(tcalcpt0,ycalcpt0(:,2));
    APO(3)=trapz(tcalcpt18,ycalcpt18(:,2));
    APO(4)=trapz(tcalcpt5,ycalcpt5(:,2));
    APO(5)=trapz(tcalcpt25,ycalcpt25(:,2));
    APO(6)=trapz(tcalcpt50,ycalcpt50(:,2));
    %Area under osteoblast curve
    AO(1)=trapz(tcalcpt1,ycalcpt1(:,3));
    AO(2)=trapz(tcalcpt0,ycalcpt0(:,3));
    AO(3)=trapz(tcalcpt18,ycalcpt18(:,3));
    AO(4)=trapz(tcalcpt5,ycalcpt5(:,3));
    AO(5)=trapz(tcalcpt25,ycalcpt25(:,3));
    AO(6)=trapz(tcalcpt50,ycalcpt50(:,3));
    %Area under osteoclast curve
    AC(1)=trapz(tcalcpt1,ycalcpt1(:,4));
    AC(2)=trapz(tcalcpt0,ycalcpt0(:,4));
    AC(3)=trapz(tcalcpt18,ycalcpt18(:,4));
    AC(4)=trapz(tcalcpt5,ycalcpt5(:,4));
    AC(5)=trapz(tcalcpt25,ycalcpt25(:,4));
    AC(6)=trapz(tcalcpt50,ycalcpt50(:,4));
    
    
    figure(4)
        oc(1)=max(ycalcpt1(:,4));
        oc(2)=max(ycalcpt0(:,4));
        oc(3)=max(ycalcpt18(:,4));
        oc(4)=max(ycalcpt5(:,4));
        oc(5)=max(ycalcpt25(:,4));
        oc(6)=max(ycalcpt50(:,4));
        plot(Gwntdose,oc,'m-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Osteoclast Cells','FontSize',12)

    figure(5)
        PO(1)=max(ycalcpt1(:,2));
        PO(2)=max(ycalcpt0(:,2));
        PO(3)=max(ycalcpt18(:,2));
        PO(4)=max(ycalcpt5(:,2));
        PO(5)=max(ycalcpt25(:,2));
        PO(6)=max(ycalcpt50(:,2));
        plot(Gwntdose,PO,'b-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Pre-osteoblast Cells','FontSize',12)

    figure(6)
        O(1)=max(ycalcpt1(:,3));
        O(2)=max(ycalcpt0(:,3));
        O(3)=max(ycalcpt18(:,3));
        O(4)=max(ycalcpt5(:,3));
        O(5)=max(ycalcpt25(:,3));
        O(6)=max(ycalcpt50(:,3));
        plot(Gwntdose,O,'b-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Osteoblast Cells','FontSize',12)
            
    figure(7)
        plot(tcalcpt1,ycalcpt1(:,1),'Linewidth',2);
        hold on
        plot(tcalcpt0,ycalcpt0(:,1),'Linewidth',2);
        plot(tcalcpt18,ycalcpt18(:,1),'Linewidth',2);
        plot(tcalcpt5,ycalcpt5(:,1),'Linewidth',2);
        plot(tcalcpt25,ycalcpt25(:,1),'Linewidth',2);
        plot(tcalcpt50,ycalcpt50(:,1),'Linewidth',2);
        xlabel('time(days)','FontSize',12) 
        ylabel('Osteocytes','FontSize',12)
        legend('-1', '0', '1.8', '5', '25', '50')
    figure(8)
        plot(tcalcpt1,ycalcpt1(:,2),'Linewidth',2);
        hold on
        plot(tcalcpt0,ycalcpt0(:,2),'Linewidth',2);
        plot(tcalcpt18,ycalcpt18(:,2),'Linewidth',2);
        plot(tcalcpt5,ycalcpt5(:,2),'Linewidth',2);
        plot(tcalcpt25,ycalcpt25(:,2),'Linewidth',2);
        plot(tcalcpt50,ycalcpt50(:,2),'Linewidth',2);
        xlabel('time(days)','FontSize',12) 
        ylabel('Pre-Osteoblasts','FontSize',12)
        legend('-1', '0', '1.8', '5', '25', '50')
    figure(9)
        plot(tcalcpt1,ycalcpt1(:,3),'Linewidth',2);
        hold on
        plot(tcalcpt0,ycalcpt0(:,3),'Linewidth',2);
        plot(tcalcpt18,ycalcpt18(:,3),'Linewidth',2);
        plot(tcalcpt5,ycalcpt5(:,3),'Linewidth',2);
        plot(tcalcpt25,ycalcpt25(:,3),'Linewidth',2);
        plot(tcalcpt50,ycalcpt50(:,3),'Linewidth',2);
        xlabel('time(days)','FontSize',12) 
        ylabel('Osteoblast Cells','FontSize',12)
        legend('-1', '0', '1.8', '5', '25', '50')
    figure(10)
        plot(tcalcpt1,ycalcpt1(:,4),'Linewidth',2);
        hold on
        plot(tcalcpt0,ycalcpt0(:,4),'Linewidth',2);
        plot(tcalcpt18,ycalcpt18(:,4),'Linewidth',2);
        plot(tcalcpt5,ycalcpt5(:,4),'Linewidth',2);
        plot(tcalcpt25,ycalcpt25(:,4),'Linewidth',2);
        plot(tcalcpt50,ycalcpt50(:,4),'Linewidth',2);
        xlabel('time(days)','FontSize',12) 
        ylabel('Osteoclast Cells','FontSize',12)
        legend('-1', '0', '1.8', '5', '25', '50')
    
    figure(11)
        plot(Gwntdose,PO./O,'b-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Pre-Osteoblasts/Osteoblast  Max Cells','FontSize',12)
    figure(12)
        plot(Gwntdose,APO./AO,'b-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Pre-Osteoblasts/Osteoblast area under the curve','FontSize',12)
    figure(13)
        plot(Gwntdose,AC./AO,'b-o','Linewidth',2);
        xlabel('Wnt-10b (Fold Change)','FontSize',12) 
        ylabel('Osteoclast/Osteoblast area under the curve','FontSize',12)



end

end
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
    alpha3adj=k(1)+k(2);
    %alpha3adj=k(2);
    beta2adj = k(3);
    K=k(4);
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
        
        [~,ycalc] = Cyclefunction(k,xdata(i),y0,N(i),cyclelength);
        
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
            %opts=odeset('NonNegative',(1:5),'reltol',1e-5);
            %opts=odeset('NonNegative',(1:5));
            opts=odeset('reltol',1e-7);
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