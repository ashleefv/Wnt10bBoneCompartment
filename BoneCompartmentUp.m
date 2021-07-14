%Carley Cook
% This code was written to estimate the change in parameters when a change
% in Wnt10b occurs
%% Data from Bennett 2005 and Bennett 2007
clear variables
close all
format long e
xdata=[-1, 5, 50]; %Wnt10b Fold Change
%xdata=[-1,1.8, 5]; %Wnt10b Fold Change
ydata=[-29.7,  69.2, 339]; % Bennet Data normalized BV/TV %.339 I ommited one set of OC data is because it is a repeated data point
%ydata=[-29.7, 36.6, 69.2]; % Bennet Data normalized BV/TV %.339 I ommited one set of OC data is because it is a repeated data point

ParamY=2;% ParamY=1 for running parameter estimation code
         % ParamY=2 for just graphing
N=[6,6,12]; %Number of Cylces
cyclelength=100; %length of cylces
tlag=14; %DDE lag from osteoblast maturation to activation
Gwntdose=0; %Dose of wnt that will be graphed
OCWnt=1; %1 to look at extra graphs
%% Guesses for the parameters
% from bioRxiv version 1
% Resnorm 1.675805403615064e+02
% kg(1)=4.961167679351240e-01;%4.961167679351240e-01;%4.961167679351240e-01; % 4.143062584626079e-01;
% 
% 
% kg(2)=6.472402955969400e-01-4.961167679351240e-01;%6.472402955969400e-01;%6.472402955969400e-01;%6.074262643621255e-01;%alpha3adj
% 
%              
% kg(3)=2.304268016306020e-03;%2.304268016306020e-03;%2.304268016306020e-03;%7.825811613971986e-03;%beta2adj
% 
% kg(4)=2.763795146691095e+01;%2.763795146691095e+01;%2.765003527667373e+01;%K

% fitting ydata(3) to N(3) = 12;
%First fit with tolerances all set to e-6
% resnorm =      2.736585398239347e+01;
% kg(1)=2.451829887973575e-01;%beta1adj
% 
% 
% kg(2)=2.044164503258383e-06;%alpha3adj = k(1) + k(2)
% 
%              
% kg(3)=6.237699914466416e-08;%beta2adj
% 
% kg(4)=8.626109236588162e+00;%K

%Final fit with function tolerance e-14 and step tolerance e-14
% resnorm =      2.681368156070732e+01
kg(1)=2.154735789777285e-01;%2.175764810196212e-01;%2.451829887973575e-01;%beta1adj


kg(2)=2.976859670464066e-02;%1.261121673336933e-04;%2.044164503258383e-06;%alpha3adj = k(1) + k(2)

             
kg(3)=1.282354008936270e-03;%3.175894315844930e-03;%6.237699914466416e-08;%beta2adj

kg(4)=8.424005104382097e+00;%8.533257716235623e+00;%8.626109236588162e+00;%K


kguess=kg;

%Setting guesses to k values if only graphing
if ParamY==2
    k=kg;
    %% rounding
% % from bioRxiv version 1
%     kr(1)=4.96e-01;%beta1adj
% 
% 
%     kr(2)=6.47e-01-4.96e-01;%6.47e-01;%alpha3adj = k(1) + k(2)
% 
%              
%     kr(3)=2.30e-03;%beta2adj
% 
%     kr(4)=2.76e+01;%K

% fitting ydata(3) to N(3) = 12;
    kr(1)=2.15e-01;%beta1adj


    kr(2)=2.98e-02;%6.47e-01;%alpha3adj = k(1) + k(2)

             
    kr(3)=1.28e-03;%beta2adj

    kr(4)=8.42e+00;%K
    
end


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
OPTIONS = optimoptions('lsqcurvefit','StepTolerance',1e-14,...
    'FunctionTolerance',1e-14,'optimalitytolerance', 1e-6,'MaxFunctionEvaluations',10000,'Algorithm','levenberg-marquardt');

%lb=[-0.1/50,-0.1,0];
lb=[0,0,0,1];
%ub=[12,1.34e-03,0.3];
ub=[1,1,1, 100];
%ub=[Inf,0.00015,Inf,100];
%[k,resnorm] = lsqcurvefit(@(k,xdata) Graham2013(k,xdata,y0,N,cyclelength,tlag),kguess,xdata,ydata, lb, ub)%, OPTIONS)
[k,resnorm,residual] = lsqcurvefit(@(k,xdata) Graham2013(k,xdata,y0,...
    N,cyclelength),kguess,xdata,ydata, lb, ub, OPTIONS)
end      
%% Graphing of BV vs Wnt
figure(1)
 %xp = linspace(xdata(1),xdata(end),1001);
 %xp = linspace(xdata(1),xdata(end),100);
 xp = linspace(-1,50,100);
 Np6=zeros(length(xp),1);
 Np12=zeros(length(xp),1);
 Np6(:)=6;
 Np12(:)=12;
 ycalcp6 = Graham2013(k,xp,y0,Np6,cyclelength);
 ycalcp12 = Graham2013(k,xp,y0,Np12,cyclelength);
 plot(xp,ycalcp6(1,:),'r','Linewidth',2);
 hold on
 plot(xp,ycalcp12(1,:),'r','Linewidth',2);
 xdata2007=[-1,5];
 ydata2007=[-29.7,  69.2];
 xdata2005=50;
 ydata2005=339;
 plot(xdata2007,ydata2007,'ko','Linewidth',2)
 plot(xdata2005,ydata2005,'bs','Linewidth',2)
% Checking that the rounded paramters produce a simular simulation result
%  if ParamY==2
%  ycalcr6 = Graham2013(kr,xp,y0,Np6,cyclelength);
%  plot(xp,ycalcr6(1,:),'b','Linewidth',2)
%  ycalcr12 = Graham2013(kr,xp,y0,Np12,cyclelength);
%  plot(xp,ycalcr12(1,:),'b','Linewidth',2)
%  end
 legend("Simulation Results 6 Cycles","Simulation Results 12 Cylces","Bennett 2007 Data","Bennett 2005 Data",'Location','Best','FontSize',12)
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


    [tcalcpt1,ycalcpt1]=Cyclefunction(k,-1,y0,N(1),cyclelength);
%     [tcalcpt18,ycalcpt18]=Cyclefunction(k,1.8,y0,12,cyclelength);
    [tcalcpt5,ycalcpt5]=Cyclefunction(k,5,y0,N(2),cyclelength);
    [tcalcpt50,ycalcpt50]=Cyclefunction(k,50,y0,N(3),cyclelength);
    figure(14)
    plot(tcalcpt1,ycalcpt1(:,5),'g-','Linewidth',2)
    hold on
    %         plot(tcalcpt18,ycalcpt18(:,5),'b:','Linewidth',2)
    plot(tcalcpt5,ycalcpt5(:,5),'b:','Linewidth',2)
    plot(tcalcpt50,ycalcpt50(:,5),'r-.','Linewidth',2)
    plot(N(1)*cyclelength,ydata(1)+100,'ko','Linewidth',2)
    plot(N(2)*cyclelength,ydata(2)+100,'ko','Linewidth',2)
    plot(N(3)*cyclelength,ydata(3)+100,'ko','Linewidth',2)
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
%     alpha3adj=k(2);
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