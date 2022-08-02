%Carley Cook
% This code was written see the change in model results for multiple
% parameter sets
%% Data from Bennett 2005 and Bennett 2007
clear variables
close all
format long e
xdata=[-1, -1, 5, 50]; %Wnt10b Fold Change
ydata=[-29.7, -41.9,  69.2, 339]; % Bennet Data normalized BV/TV %.339 
N=[4,6,6,12]; %Number of Cycles
cyclelength=100; %length of cycles
Gwntdose=50; %Dose of wnt that will be graphed
%% Loading in the parameter results from mutlistart lsqcurve fit
Rparams1=load('BoneCompartmentUp4ParaRelTOL7.mat','Obj_best','p_best', 'kguessStore');
Rparams2=load('BoneCompartmentUp4ParaRand.mat','Obj_best','p_best', 'kguessStore');
Rparams3=load('BoneCompartmentUp4ParaReduced.mat','Obj_best','p_best', 'kguessStore');
Rparams4=load('BoneCompartmentUp4ParaNarrow.mat','Obj_best','p_best', 'kguessStore');
Rparams5=load('BoneCompartmentUp4ParaNarrow2.mat','Obj_best','p_best', 'kguessStore');
Rparams6=load('BoneCompartmentUp4ParaK2Narrow.mat','Obj_best','p_best', 'kguessStore');
Rparams7=load('BoneCompartmentUp4ParaMeanNarrow.mat','Obj_best','p_best', 'kguessStore');
Obj_best=[Rparams1.Obj_best; Rparams2.Obj_best; Rparams3.Obj_best; Rparams4.Obj_best; Rparams5.Obj_best;  Rparams6.Obj_best; Rparams7.Obj_best]; %Resnorm of each run
p_best=[Rparams1.p_best; Rparams2.p_best; Rparams3.p_best; Rparams4.p_best; Rparams5.p_best; Rparams6.p_best; Rparams7.p_best]; %Lsqcurve fit paramerers
kguessStore=[Rparams1.kguessStore; Rparams2.kguessStore; Rparams3.kguessStore; Rparams4.kguessStore; Rparams5.kguessStore; Rparams6.kguessStore; Rparams7.kguessStore]; %initial guess parameters

%Sorted matrix of guesses
Rguess=[Obj_best, kguessStore, p_best];
Rguess=sortrows(Rguess,1);

%Sorted matrix of results
Rparams=[Obj_best, p_best];
RparamsS = sortrows(Rparams,1);
cutoff=RparamsS(1,1)+RparamsS(1,1)*.1; 
RparamsS(RparamsS(:,1)>cutoff, :)=[]; %Reducing the results to only ones within 1% of minimum resnorm
RparamsT=RparamsS; %Stores all the results
RparamsS(RparamsS(:,2)<.0001, :)=[];%Cutting results that cause params to be less than 10-5
RparamsS(RparamsS(:,3)<.0001, :)=[];
RparamsS(RparamsS(:,4)<.0001, :)=[];
RparamsS(RparamsS(:,5)<.0001, :)=[];
RparamsS(RparamsS(:,3)>=.1, :)=[];


k=RparamsS(:,2:end); %Storing the params without resnorm
 %% Calculating the average k and the standard deviation
  kselect=k(1:20,:)
  k=kselect
 paramavg=mean(k); %Results [0.122021014367919	0.174321892073088	2.16426217564921e-05	5.64691894915150]
 parasd=std(k);% Results [0.000220459131050966	0.000177680002087373	1.19185218536249e-05	0.0535222520019603]

 %k=[normrnd(paramavg(1),parasd(1),[100,1]), normrnd(paramavg(2),parasd(2),[100,1]), normrnd(paramavg(3),parasd(3),[100,1]), normrnd(paramavg(4),parasd(4),[100,1]),];
 %k(k(:,3)<0,3)=0;
%% 
 
 
 %  Checking SD
 SD3CHECKu=paramavg
 SD3CHECKu(1,3)=SD3CHECKu(1,3)+parasd(1,3)
 SD3CHECKl=paramavg
 SD3CHECKl(1,3)=SD3CHECKl(1,3)-parasd(1,3)
%  for l=1:height(k) %creating matrixes that include changes in sd for each param seperately
%  paramup(l,:)=paramavg;
%  paramlow(l,:)=paramavg;
%  paramup(l,l)=paramup(l,l)+parasd(l);
%  paramlow(l,l)=paramlow(l,l)-parasd(l);
%  end
%% Initial conditions 
% S0=200 for SS and K_S-rho=(180) for activation
S0=180; % Initial Osteocytes 
P0=0; %Initial Pre-Osteoblasts
B0=0; %Initial Osteoblasts
C0=0; %Initial Osteoclasts
z0=100; %Initial Bone Volume
y0=[S0,P0,B0,C0,z0]; % Initial conditions in vector    
%% Visually inspecting the parameter results
    % Simulation at the average k values
    [tcalcpt1avg, ycalcpt1avg]=Cyclefunction(paramavg,-1,y0,N(2),cyclelength);
    [tcalcpt5avg, ycalcpt5avg]=Cyclefunction(paramavg,5,y0,N(3),cyclelength);
    [tcalcpt50avg, ycalcpt50avg]=Cyclefunction(paramavg,50,y0,N(4),cyclelength);

    %Checking upper and lower sd of k3
%     [tcalcpt1avg2, ycalcpt1avg2]=Cyclefunction(SD3CHECKu,-1,y0,N(2),cyclelength);
%     [tcalcpt5avg2, ycalcpt5avg2]=Cyclefunction(SD3CHECKu,5,y0,N(3),cyclelength);
%     [tcalcpt50avg2, ycalcpt50avg2]=Cyclefunction(SD3CHECKu,50,y0,N(4),cyclelength);
% 
%     [tcalcpt1avg3, ycalcpt1avg3]=Cyclefunction(SD3CHECKl,-1,y0,N(2),cyclelength);
%     [tcalcpt5avg3, ycalcpt5avg3]=Cyclefunction(SD3CHECKl,5,y0,N(3),cyclelength);
%     [tcalcpt50avg3, ycalcpt50avg3]=Cyclefunction(SD3CHECKl,50,y0,N(4),cyclelength);
    
    for i=1:height(k) %Different simulation senarios
    %Simulation for the upper and lower sd of each k value
    [tcalcpt1avgup(:,i), ycalcpt1avgup(:,5*i-4:5*i)]=Cyclefunction(paramup(i,:),-1,y0,N(2),cyclelength);
    [tcalcpt5avgup(:,i), ycalcpt5avgup(:,5*i-4:5*i)]=Cyclefunction(paramup(i,:),5,y0,N(3),cyclelength);
    [tcalcpt50avgup(:,i), ycalcpt50avgup(:,5*i-4:5*i)]=Cyclefunction(paramup(i,:),50,y0,N(4),cyclelength);
    [tcalcpt1avglow(:,i), ycalcpt1avglow(:,5*i-4:5*i)]=Cyclefunction(paramlow(i,:),-1,y0,N(2),cyclelength);
    [tcalcpt5avglow(:,i), ycalcpt5avglow(:,5*i-4:5*i)]=Cyclefunction(paramlow(i,:),5,y0,N(3),cyclelength);
    [tcalcpt50avglow(:,i), ycalcpt50avglow(:,5*i-4:5*i)]=Cyclefunction(paramlow(i,:),50,y0,N(4),cyclelength);

    %Simulation for all of the selected ks
    [tcalcpt1(:,i), ycalcpt1(:,5*i-4:5*i)]=Cyclefunction(k(i,:),-1,y0,N(2),cyclelength);
    [tcalcpt5(:,i), ycalcpt5(:,5*i-4:5*i)]=Cyclefunction(k(i,:),5,y0,N(3),cyclelength);
    [tcalcpt25(:,i), ycalcpt25(:,5*i-4:5*i)]=Cyclefunction(k(i,:),25,y0,N(4),cyclelength);
    [tcalcpt50(:,i), ycalcpt50(:,5*i-4:5*i)]=Cyclefunction(k(i,:),50,y0,N(4),cyclelength);

    %Simulation for the k values for validation values
    [tcalcpt12(:,i),ycalcpt12(:,5*i-4:5*i)]=Cyclefunction(k(i,:),1.8,y0,12,cyclelength);
    [tcalcpt12top(:,i),ycalcpt12top(:,5*i-4:5*i)]=Cyclefunction(k(i,:),2.4,y0,12,cyclelength);
    [tcalcpt12bottom(:,i),ycalcpt12bottom(:,5*i-4:5*i)]=Cyclefunction(k(i,:),1.2,y0,12,cyclelength);


    %resnorm(i,1)=(ycalcpt1(2201,5*i)-(ydata(1)+100))^2+(ycalcpt1(3301,5*i)-(ydata(2)+100))^2+(ycalcpt5(3301,5*i)-(ydata(3)+100))^2+(ycalcpt50(6601,5*i)-(ydata(4)+100))^2;
    end



figure(1) %Plotting the average k and the selected k results
    plot(tcalcpt1avg,ycalcpt1avg(:,5),'r','Linewidth',2)
    hold on
    plot(tcalcpt5avg,ycalcpt5avg(:,5),'r','Linewidth',2)
    plot(tcalcpt50avg,ycalcpt50avg(:,5),'r','Linewidth',2)
    for i=1:height(k) %Plotting for each k
    v=rand(3,1); %Randomize the color for each run the color will be the same for each ks -1, 5, and 50 fold change
    txt=['run=',num2str(i)];
    q(i)=plot(tcalcpt1(:,1),ycalcpt1(:,5*i),'Linewidth',2,'DisplayName',txt,'color',v);
    plot(tcalcpt5(:,1),ycalcpt5(:,5*i),'Linewidth',2,'color',v)
    plot(tcalcpt50(:,1),ycalcpt50(:,5*i),'Linewidth',2,'color',v)
    end
    
    errorbar(N(1)*cyclelength,ydata(1)+100,19.7,'ko','Linewidth',2)
    errorbar(N(2)*cyclelength,ydata(2)+100,13.5,'ko','Linewidth',2)
    errorbar(N(3)*cyclelength,ydata(3)+100,40.1,'ko','Linewidth',2)
    errorbar(N(4)*cyclelength,ydata(4)+100,114,'ko','Linewidth',2)
    xlabel('time(days)','FontSize',15)
    ylabel('Relative bone volume (%)','FontSize',15)
    legend(q)


figure(2) %Plotting the average k and the standard deviation of each average k
    plot(tcalcpt1avg,ycalcpt1avg(:,5),'r','Linewidth',2)
    hold on
    plot(tcalcpt5avg,ycalcpt5avg(:,5),'r','Linewidth',2)
    plot(tcalcpt50avg,ycalcpt50avg(:,5),'r','Linewidth',2)
    for i=1:height(k)
    v=rand(3,1);
    txt=['Paramchange ',num2str(i)];
    o(i)=plot(tcalcpt1avgup(:,1),ycalcpt1avgup(:,5*i),'Linewidth',2,'DisplayName',txt,'color',v);
    plot(tcalcpt5avgup(:,1),ycalcpt5avgup(:,5*i),'Linewidth',2,'color',v)
    plot(tcalcpt50avgup(:,1),ycalcpt50avgup(:,5*i),'Linewidth',2,'color',v)
    plot(tcalcpt1avglow(:,1),ycalcpt1avglow(:,5*i),'Linewidth',2,'color',v)
    plot(tcalcpt5avglow(:,1),ycalcpt5avglow(:,5*i),'Linewidth',2,'color',v)
    plot(tcalcpt50avglow(:,1),ycalcpt50avglow(:,5*i),'Linewidth',2,'color',v)
    end

    errorbar(N(1)*cyclelength,ydata(1)+100,19.7,'ko','Linewidth',2)
    errorbar(N(2)*cyclelength,ydata(2)+100,13.5,'ko','Linewidth',2)
    errorbar(N(3)*cyclelength,ydata(3)+100,40.1,'ko','Linewidth',2)
    errorbar(N(4)*cyclelength,ydata(4)+100,114,'ko','Linewidth',2)
    legend(o)


 figure(3) %Plotting the average k for the validation case
        for i=1:height(k)
        v=rand(3,1);
        txt=['Paramgroup ',num2str(i)];
        g(i)=plot(tcalcpt12(:,1),ycalcpt12(:,5*i),'Linewidth',2,'DisplayName',txt,'color',v);
        hold on
        shade(tcalcpt12top(:,1),ycalcpt12top(:,5*i),tcalcpt12bottom(:,1),ycalcpt12bottom(:,5*i),'FillType',[1 2;2 1],'Color', v);
        end
        errorbar( 600 , 126.6 , 15,'o','Linewidth',2 )
        errorbar( 1200 , 136.6 , 29,'s','Linewidth',2 )
        legend(g)
        xlabel('Time (days)','FontSize',12) 
        ylabel('Relative bone volume (%)','FontSize',12)
        xlim([-10 1210])


% figure(4) %Checking k3 sd
%     plot(tcalcpt1avg,ycalcpt1avg(:,5),'r','Linewidth',2)
%     hold on
%     plot(tcalcpt5avg,ycalcpt5avg(:,5),'r','Linewidth',2)
%     plot(tcalcpt50avg,ycalcpt50avg(:,5),'r','Linewidth',2)
% 
%     plot(tcalcpt1avg2,ycalcpt1avg2(:,5),'g','Linewidth',2)
%     
%     plot(tcalcpt5avg2,ycalcpt5avg2(:,5),'g','Linewidth',2)
%     plot(tcalcpt50avg2,ycalcpt50avg2(:,5),'g','Linewidth',2)
%     plot(tcalcpt1avg3,ycalcpt1avg3(:,5),'b','Linewidth',2)
%     plot(tcalcpt5avg3,ycalcpt5avg3(:,5),'b','Linewidth',2)
%     plot(tcalcpt50avg3,ycalcpt50avg3(:,5),'b','Linewidth',2)
%     
%     errorbar(N(1)*cyclelength,ydata(1)+100,19.7,'ko','Linewidth',2)
%     errorbar(N(2)*cyclelength,ydata(2)+100,13.5,'ko','Linewidth',2)
%     errorbar(N(3)*cyclelength,ydata(3)+100,40.1,'ko','Linewidth',2)
%     errorbar(N(4)*cyclelength,ydata(4)+100,114,'ko','Linewidth',2)
%     xlabel('time(days)','FontSize',15)
%     ylabel('Relative bone volume (%)','FontSize',15)
%     legend("Average","","","Upper SD","","","Lower SD")
        


%% Plotting parameter guesses and parameter results

% figure(4)
% plot(Rguess(:,2)./Rguess(:,6),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP1/P1','FontSize',12)
% 
% figure(5)
% plot(Rguess(:,3)./Rguess(:,7),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP2/P2','FontSize',12)
% 
% figure(6)
% plot(Rguess(:,4)./Rguess(:,8),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP3/P3','FontSize',12)
% 
% figure(7)
% plot(Rguess(:,5)./Rguess(:,9),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP4/P4','FontSize',12)
% 
% figure(8)
% plot(Rguess(:,2)./Rguess(:,3),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP1/GP2','FontSize',12)
% 
% figure(9)
% plot(Rguess(:,2)./Rguess(:,4),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP1/GP3','FontSize',12)
% 
% figure(10)
% plot(Rguess(:,2)./Rguess(:,5),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP1/GP4','FontSize',12)
% 
% figure(11)
% plot(Rguess(:,3)./Rguess(:,4),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP2/GP3','FontSize',12)
% 
% figure(12)
% plot(Rguess(:,3)./Rguess(:,5),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP2/GP4','FontSize',12)
% 
% figure(13)
% plot(Rguess(:,4)./Rguess(:,5),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('GP3/GP4','FontSize',12)
% 
% figure(14)
% plot(Rguess(:,6)./Rguess(:,7),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('P1/P2','FontSize',12)
% 
% figure(15)
% plot(Rguess(:,6)./Rguess(:,8),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('P1/P3','FontSize',12)
% 
% figure(16)
% plot(Rguess(:,6)./Rguess(:,9),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('P1/P4','FontSize',12)
% 
% figure(17)
% plot(Rguess(:,7)./Rguess(:,8),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('P2/P3','FontSize',12)
% 
% figure(18)
% plot(Rguess(:,7)./Rguess(:,9),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('P2/P4','FontSize',12)
% 
% figure(19)
% plot(Rguess(:,8)./Rguess(:,9),'o')
% xlabel('Increasing Resnorm','FontSize',12)
% ylabel('P3/P4','FontSize',12)

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
            %opts=odeset('NonNegative',(1:5),'reltol',1e-7);
            opts=odeset('reltol',1e-7);
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