function [error] = costfunction(param_in,ready)

error = 0;

global BartData
global BartData2
global BartData4
global pNamesOpt  
global icOrig 
global modelName
global maxsim
global minsim
global maxsim2
global minsim2
global maxsim3
global minsim3
global maxsim4
global minsim4
global maxsim5
global minsim5

%% simulation options
simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e4;

%% input parameters, saved as log(values)
cost=param_in(1);
kscale=param_in(2);
param=param_in(3:end);

param_test = exp(param(1:end-7));
param_BMP4_ba1 = exp(param(end-6)); %BMP4 effects on ucp1
param_BMP4_ba2 = exp(param(end-5)); %BMP4 effects on ETC/FA
param_BMP4_ba3 = exp(param(end-4)); %BMP4 effects on Htot/nonmit
param_scr = exp(param(end-3)); %depletion experiments
param_ba1 = exp(param(end-2)); %ucp1
param_ba2 = exp(param(end-1)); %ETC/FA
param_ba3 = exp(param(end)); %Htot/nonmit

p = zeros(7,1);

if strcmp(modelName,'model_H1') || strcmp(modelName,'model_H3')
    ucp1_bmp4=0; % for model_H1 and model_H3
elseif strcmp(modelName,'model_H2')
    ucp1_bmp4=param_ba1*param_BMP4_ba1; % for model_H2
end

ucp1_rosi=param_ba1;
depl=1;

%% define simulation times (corresponding to experimental time points)
time = BartData.Set1.time(1:13);
    injection = [18 47];
    sim_time{1} = [time(4):time(9)]-injection(1); %[time(4:9)-injection(1)]; %
    sim_time{2} = [time(9):time(13)+2]-injection(2); % [time(9:13)-injection(2)]; %
time2 = [1;BartData4.time];
    injection2  = [13 33];
sim_time2{1} = [time2(4):time2(8)+3]-injection2(1); %[time(4:9)-injection(1)]; %
sim_time2{2} = [time2(8)+3:time2(14)]-injection2(2); % [time(9:13)-injection(2)]; %

%% all simulations should be done six times
% 0) WAT witout iso, 1)BMP4 without iso, 2) Rosi without iso,
% 3) WAT with iso, 4) BMP4 with iso, 5) Rosi with iso
A = [0 2 4 1 3 5]; %order
for i = A
iso = 0;
oligo = 0;
atra = 0;
ba1=0;
ba2=0;
ba3=0;
if i == 0
    bmp4 = 0;
    ba1 = 0;
elseif i == 3
    bmp4 = 0;    
    ba1 = 0;
elseif i == 2
    bmp4 = 0;
    ba3 = param_ba3; %Htot/nonmit
elseif i == 5
    bmp4 = 0;
    ba1 = ucp1_rosi; %ucp1
    ba2 = param_ba2; %ETC/FA
    ba3 = param_ba3; %Htot/nonmit 
elseif  i == 1
    bmp4 = 1;
    ba3 = param_ba3*param_BMP4_ba3; %Htot/nonmit
elseif  i == 4
    bmp4 = 1;
    ba1 = ucp1_bmp4; %ucp1
    ba2 = param_ba2*param_BMP4_ba2; %ETC/FA
    ba3 = param_ba3*param_BMP4_ba3; %Htot/nonmit
end
paramSS = [param_test depl iso oligo ba1 ba2 ba3 bmp4 atra]; %parameters for steady state simulation
oligo = 1;
paramOLIGO2 = [param_test depl iso oligo ba1 ba2 ba3 bmp4 atra]; %parameters for oligomycin simulation
if i <= 2
    iso = 0;
else
    iso = 1;  
end
paramISO2 = [param_test depl iso oligo ba1 ba2 ba3 bmp4 atra]; %parameters for isoproterenol simulation
iso=0;
if i <= 2
    atra = 0;
else
    atra = 1;  
end
paramATRA = [param_test depl iso oligo ba1 ba2 ba3 bmp4 atra]; %parameters for all-trans retinoic acid simulation
atra=0;
oligo = 0;
if i <= 2
    iso = 0;
else
    iso = 1;  
end
paramISO = [param_test depl iso oligo ba1 ba2 ba3 bmp4 atra]; %parameters for isoproterenol simulation (when after oligo)
oligo = 1;
paramOLIGO = [param_test depl iso oligo ba1 ba2 ba3 bmp4 atra]; %parameters for oligomycin simulation (when before iso)

%% simulate the model with the given parameters
try
  simSS = IQMPsimulate(modelName,0:1:10000,icOrig,pNamesOpt,paramSS,simOptions);
catch
    disp('Simulation 0 crashed...');
    error = inf;
    return 
end

initcondOLIGO2=simSS.statevalues(end,:); %update initial conditions
try
    simOLIGO2 =IQMPsimulate(modelName,sim_time2{1},initcondOLIGO2,pNamesOpt,paramOLIGO2,simOptions);
catch
    disp('Simulation 2 crashed...');
    error = inf;
    return 
end

initcondISO2=simOLIGO2.statevalues(end,:); %update initial conditions
try
    simISO2 =IQMPsimulate(modelName,sim_time2{2},initcondISO2,pNamesOpt,paramISO2,simOptions);
catch
    disp('Simulation 1 crashed...');
    error = inf;
    return 
end

initcondATRA=simOLIGO2.statevalues(end,:); %update initial conditions
try
    simATRA =IQMPsimulate(modelName,sim_time2{2},initcondATRA,pNamesOpt,paramATRA,simOptions);
catch
    disp('Simulation 1 crashed...');
    error = inf;
    return 
end

initcondISO=simSS.statevalues(end,:); %update initial conditions
try
    simISO =IQMPsimulate(modelName,sim_time{1},initcondISO,pNamesOpt,paramISO,simOptions);
catch
    disp('Simulation 1 crashed...');
    error = inf;
    return 
end

initcondOLIGO=simISO.statevalues(end,:); %update initial conditions
try
    simOLIGO =IQMPsimulate(modelName,sim_time{2},initcondOLIGO,pNamesOpt,paramOLIGO,simOptions);
catch
    disp('Simulation 2 crashed...');
    error = inf;
    return 
end
    SS = simSS.variablevalues(end-3*6:end,end)';
    ISO = simISO.variablevalues(:,end)';
    OLIGO = simOLIGO.variablevalues(:,end)';
    
    OCRsim = [SS ISO(2:end) OLIGO(2:end)]';
    OCRsim2 = [simSS.variablevalues(end-13:end,end)' simOLIGO2.variablevalues(2:end,end)' simISO2.variablevalues(2:end,end)']';
    OCRsim3 = [simSS.variablevalues(end-13:end,end)' simOLIGO2.variablevalues(2:end,end)' simATRA.variablevalues(2:end,end)']';

%% collection and scaling of simulation results
    if i==0
        SSnorm = SS(1);
        Data = [BartData.Set1.time BartData.Set1.BasDMSO BartData.Set1.BasDMSOstd];
        Data = Data(2:end-1,:);
    elseif  i==1
        SSnorm = SS(1);
        Data = [BartData.Set1.time BartData.Set1.BMP4DMSO BartData.Set1.BMP4DMSOstd];
        Data = Data(2:end-1,:);
        Data2 = [BartData4.time BartData4.Dmso' BartData4.Dmsosem'];
        Data3 = [BartData4.time BartData4.DmsoAtra' BartData4.DmsoAtrasem'];
    elseif i==2
        SSnorm = SS(1);
        Data = [BartData.Set1.time BartData.Set1.RosiDMSO BartData.Set1.RosiDMSOstd];
        Data = Data(2:end-1,:);
        Data2 = [BartData2.time(2:end) BartData2.Dmso(2:end)./BartData2.Dmso(4).*100 BartData2.Dmsosem(2:end)./BartData2.Dmso(4).*100];
        Data3 = [BartData2.time(2:end) BartData2.Dmso(2:end)./BartData2.Dmso(4).*100 BartData2.Dmsosem(2:end)./BartData2.Dmso(4).*100];
    elseif i==3
        Data = [BartData.Set1.time BartData.Set1.BasIso BartData.Set1.BasIsostd];
        Data = Data(2:end-1,:);
    elseif i==4
        Data = [BartData.Set1.time BartData.Set1.BMP4Iso BartData.Set1.BMP4Isostd];
        Data = Data(2:end-1,:);
        Data2 = [BartData4.time BartData4.BMP4' BartData4.BMP4sem'];
        Data3 = [BartData4.time BartData4.BMP4Atra' BartData4.BMP4Atrasem'];
    elseif i==5
        Data = [BartData.Set1.time BartData.Set1.RosiIso BartData.Set1.RosiIsostd];
        Data = Data(2:end-1,:);
        Data2 = [BartData2.time(2:end) BartData2.Rosi(2:end)./BartData2.Rosi(4).*100 BartData2.Rosisem(2:end)./BartData2.Rosi(4).*100];
        Data3 = [BartData2.time(2:end) BartData2.RosiAtra(2:end)./BartData2.RosiAtra(4).*100 BartData2.RosiAtrasem(2:end)./BartData2.RosiAtra(4).*100];
    end

if i==1 || i==4
   for t = 1:73
        if OCRsim(t)/OCRsim(1)*100>maxsim(i+1,t)
            maxsim(i+1,t)=OCRsim(t)/OCRsim(1)*100;
        elseif OCRsim(t)/OCRsim(1)*100<minsim(i+1,t)
            minsim(i+1,t)=OCRsim(t)/OCRsim(1)*100;
        end
   end
else
   for t = 1:73
        if kscale*OCRsim(t)>maxsim(i+1,t)
            maxsim(i+1,t)=kscale*OCRsim(t);
        elseif kscale*OCRsim(t)<minsim(i+1,t)
            minsim(i+1,t)=kscale*OCRsim(t);
        end
   end
end
   for t = 1:55
        if OCRsim2(t)/OCRsim2(1)*100>maxsim2(i+1,t)
            maxsim2(i+1,t)=OCRsim2(t)/OCRsim2(1)*100;
        elseif OCRsim2(t)/OCRsim2(1)*100<minsim2(i+1,t)
            minsim2(i+1,t)=OCRsim2(t)/OCRsim2(1)*100;
        end
        if OCRsim3(t)/OCRsim3(1)*100>maxsim3(i+1,t)
            maxsim3(i+1,t)=OCRsim3(t)/OCRsim3(1)*100;
        elseif OCRsim3(t)/OCRsim3(1)*100<minsim3(i+1,t)
            minsim3(i+1,t)=OCRsim3(t)/OCRsim3(1)*100;
        end
    end
    
%% simulate the scramlbe experiment
        if i >= 2 && i<3 || i>4    
        ucp1_scr=1;
        bmp4 = 0;
        paramSS_scr = paramSS;
        paramSS_scr(end-7) = param_scr;
        paramSS_scr(end-4) = ucp1_scr;
        paramISO_scr = paramISO;
        paramISO_scr(end-7) = param_scr;
        paramISO_scr(end-4) = ucp1_scr;
        paramOLIGO_scr = paramOLIGO;
        paramOLIGO_scr(end-7) = param_scr;
        paramOLIGO_scr(end-4) = ucp1_scr;
        
        try
            simSS_scr =IQMPsimulate(modelName,0:1:10000,icOrig,pNamesOpt,paramSS_scr,simOptions);
        catch
            disp('Simulation 1 crashed...');
            error = inf;
            return
        end
        SS_scr = simSS_scr.variablevalues(end-3*6:end,end)';
        initicondISO_scr =  simSS_scr.statevalues(end,:);     
        try
            simISO_scr =IQMPsimulate(modelName,sim_time{1},initcondISO,pNamesOpt,paramISO_scr,simOptions);
        catch
            disp('Simulation 1 crashed...');
            error = inf;
            return
        end
        ISO_scr = simISO_scr.variablevalues(:,end)';
        initcondOLIGO_scr = simISO_scr.statevalues(end,:);
        try
            simOLIGO_scr = IQMPsimulate(modelName,sim_time{2},initcondOLIGO_scr,pNamesOpt,paramOLIGO_scr,simOptions);
        catch
            disp('Simulation 1 crashed...');
            error = inf;
            return
        end
        OLIGO_scr = simOLIGO_scr.variablevalues(:,end)';
        end

    
%% simulate the depletion experiment
       if i >= 2 && i<3 || i>4 
        ucp1_depletion=0.1;
        bmp4 = 0;
        paramSS_dep = paramSS;
        paramSS_dep(end-7) = param_scr;
        paramSS_dep(end-4) = ucp1_depletion;
        paramISO_dep = paramISO;
        paramISO_dep(end-7) = param_scr;
        paramISO_dep(end-4) = ucp1_depletion;
        paramOLIGO_dep = paramOLIGO;
        paramOLIGO_dep(end-7) = param_scr;
        paramOLIGO_dep(end-4) = ucp1_depletion;
        
        try
            simSS_dep =IQMPsimulate(modelName,0:1:10000,icOrig,pNamesOpt,paramSS_dep,simOptions);
        catch
            disp('Simulation 1 crashed...');
            error = inf;
            return
        end
        SS_dep = simSS_dep.variablevalues(end-3*6:end,end)';
        initicondISO_dep =  simSS_dep.statevalues(end,:);     
        try
            simISO_dep =IQMPsimulate(modelName,sim_time{1},initcondISO,pNamesOpt,paramISO_dep,simOptions);
        catch
            disp('Simulation 1 crashed...');
            error = inf;
            return
        end
        ISO_dep = simISO_dep.variablevalues(:,end)';
        initcondOLIGO_dep = simISO_dep.statevalues(end,:);
        try
            simOLIGO_dep = IQMPsimulate(modelName,sim_time{2},initcondOLIGO_dep,pNamesOpt,paramOLIGO_dep,simOptions);
        catch
            disp('Simulation 1 crashed...');
            error = inf;
            return
        end
        OLIGO_dep = simOLIGO_dep.variablevalues(:,end)';

       
    	OCRsim4=100*([SS_scr ISO_scr(2:end) OLIGO_scr(2:end)])./SS_scr(1);        
        OCRsim5=100*([SS_dep ISO_dep(2:end) OLIGO_dep(2:end)])./SS_dep(1);
        for t = 1:73
        if OCRsim4(t)/OCRsim4(1)*100>maxsim4(i+1,t)
            maxsim4(i+1,t)=OCRsim4(t)/OCRsim4(1)*100;
        elseif OCRsim4(t)/OCRsim4(1)*100<minsim4(i+1,t)
            minsim4(i+1,t)=OCRsim4(t)/OCRsim4(1)*100;
        end
        if OCRsim5(t)/OCRsim5(1)*100>maxsim5(i+1,t)
            maxsim5(i+1,t)=OCRsim5(t)/OCRsim5(1)*100;
        elseif OCRsim5(t)/OCRsim5(1)*100<minsim5(i+1,t)
            minsim5(i+1,t)=OCRsim5(t)/OCRsim5(1)*100;
        end
        end
        end
        
%% plot all results

        greenDMSO = [0 0.6 0];
        greenISO = [0 0.2 0];
        redDMSO = [1 0 0];
        redISO = [0.6 0 0];
        brownDMSO = [0.8 0.6 0.4];
        brownISO = [0.4 0.2 0];
        color_uncert=[0.8 0.8 0.8];
        
        if i<1
            color=greenDMSO;
            marker='o';
        elseif i<=1 && i<2
            color=brownDMSO;
            marker='o';
        elseif i >= 2 && i<3
            color=redDMSO;
            marker='^';
        elseif i >=3 && i<4
            color=greenISO;
            marker='s';
        elseif i>=4 && i<5
            color=brownISO;
            marker='s';
        else
            color = redISO;
             marker='v';
        end
        
        plottime = [1:72];
        plottime2 = [1:55];
        
    if (nargin>1) && (ready == 1)

        if i==0 || i==3 || i==2 || i==5 
        figure(6)
        errorbar(time(2:13), Data(1:12,2), Data(1:12,3), 'Color', color, 'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',8)
        hold on
        p(i+1) = plot(1:72, kscale*OCRsim(1:72), 'Color', color,'Linewidth',2);
        xlabel('Time (min)','FontSize',20)
        ylabel('OCR (pmol/min) / Cells','FontSize',20)
        set(gca,'TickDir','out','YTick',[0 0.01 0.02 0.03 0.04 0.05],'YTickLabel',[0 10 20 30 40 50],'XTick',[1 12 24 36 47 58 70 81 95 107 119],'FontSize',20);
        axis([0 71 0 0.04])
        box off
        end
        
        if i==2 || i==5
        figure(1)
        subplot(1,2,1)
        plot(1:73, 100*([SS_scr ISO_scr(2:end) OLIGO_scr(2:end)])./SS_scr(1),'Color', color,'LineWidth',2)
        hold on
        xlabel('Time (min)','FontSize',17)
        ylabel('OCR (pmol/min) / Cells','FontSize',17)
        set(gca,'TickDir','out','YTick',[40 80 120 160],'XTick',[1 12 24 36 47 58 70 81 95 107 119],'FontSize',17);
        axis([0 71 40 160])
        box off
        
        subplot(1,2,2)
        plot(1:73, 100*([SS_dep ISO_dep(2:end) OLIGO_dep(2:end)])./SS_dep(1),'Color', color,'LineWidth',2)
        hold on
        xlabel('Time (min)','FontSize',17)
        ylabel('OCR (pmol/min) / Cells','FontSize',17)
        set(gca,'TickDir','out','YTick',[40 80 120 160],'XTick',[1 12 24 36 47 58 70 81 95 107 119],'FontSize',17);
        axis([0 71 40 160])
        box off
        
        figure(3)
        errorbar(Data2(:,1),Data2(:,2),Data2(:,3), 'Color', color, 'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',8)
        hold on
        p(i+1) = plot(plottime2, OCRsim2./OCRsim2(1).*100, 'Color', color,'Linewidth',2);
        xlabel('Time (min)','FontSize',20)
        ylabel('OCR (% Baseline)','FontSize',20)
        set(gca,'TickDir','out','YTick',[0 50 100 150 200],'XTick',[9 18 26 35 43 51],'FontSize',20);
        axis([0 54 0 130])
        box off
        
        figure(5)
        errorbar(Data3(:,1),Data3(:,2),Data3(:,3), 'Color', color, 'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',8)
        hold on
        p(i+1) = plot(1:55, OCRsim3./OCRsim3(1).*100, 'Color', color,'Linewidth',2);
        xlabel('Time (min)','FontSize',20)
        ylabel('OCR (% Baseline)','FontSize',20)
        set(gca,'TickDir','out','YTick',[0 50 100 150 200],'XTick',[9 18 26 35 43 51],'FontSize',20);
        axis([0 55 0 130])
        box off
        end
        
        if i==1 || i==4
        figure(4)
        errorbar(Data2(:,1),Data2(:,2),Data2(:,3), 'Color', color, 'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',8)
        hold on
        p(i+1) = plot(plottime2, OCRsim2./OCRsim2(1).*100, 'Color', color,'Linewidth',2);
        xlabel('Time (min)','FontSize',20)
        ylabel('OCR (% Baseline)','FontSize',20)
        set(gca,'TickDir','out','YTick',[0 50 100 150 200],'XTick',[9 18 26 35 43 51],'FontSize',20);
        axis([0 54 0 120])
        box off
        
        figure(7)
        errorbar(time(2:13), Data(1:12,2)./Data(1,2).*100, Data(1:12,3)./Data(1,2).*100, 'Color', color, 'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',8)
        hold on
        p(i+1) = plot(1:72, OCRsim(1:72)./OCRsim(1).*100, 'Color', color,'Linewidth',2);
        xlabel('Time (min)','FontSize',20)
        ylabel('OCR (pmol/min) / Cells','FontSize',20)
        set(gca,'TickDir','out','YTick',[0 50 100 150 200],'XTick',[1 12 24 36 47 58 70 81 95 107 119],'FontSize',20);
        axis([0 71 0 210])
        box off
        
        figure(8)
        errorbar(Data3(:,1),Data3(:,2),Data3(:,3), 'Color', color, 'LineStyle', 'n', 'Marker',marker,'MarkerFaceColor',color, 'MarkerSize',8)
        hold on
        p(i+1) = plot(1:55, OCRsim3./OCRsim3(1).*100, 'Color', color,'Linewidth',2);
        xlabel('Time (min)','FontSize',20)
        ylabel('OCR (% Baseline)','FontSize',20)
        set(gca,'TickDir','out','YTick',[0 50 100 150 200],'XTick',[9 18 26 35 43 51],'FontSize',20);
        axis([0 55 0 120])
        box off
        end        
    end

end

end
