function [] = plotscript(minsim,minsim2,minsim3,minsim4,minsim5,maxsim,maxsim2,maxsim3,maxsim4,maxsim5)

global BartData
global BartData2
global BartData4
global pNamesOpt  
global icOrig 
global modelName
global stateNames
global variableNames
global FID

time = BartData.Set1.time(1:13);
    injection = [18 47];
    sim_time{1} = [time(4):time(9)]-injection(1); 
    sim_time{2} = [time(9):time(13)+2]-injection(2);
time2 = [1;BartData4.time];
    injection2  = [13 33];
sim_time2{1} = [time2(4):time2(8)+3]-injection2(1); 
sim_time2{2} = [time2(8)+3:time2(14)]-injection2(2); 
    
% 0) WAT witout iso, 1)BMP4 without iso, 2) Rosi without iso,
% 3) WAT with iso, 4) BMP4 with iso, 5) Rosi with iso
A = [0 2 4 1 3 5]; %order
for i = A    
    normOCRsim_max = maxsim(i+1,:);
    normOCRsim2_max = maxsim2(i+1,:);
    normOCRsim3_max = maxsim3(i+1,:); 
    normOCRsim4_max = maxsim4(i+1,:); 
    normOCRsim5_max = maxsim5(i+1,:); 
    normOCRsim_min = minsim(i+1,:);
    normOCRsim2_min = minsim2(i+1,:);
    normOCRsim3_min = minsim3(i+1,:); 
    normOCRsim4_min = minsim4(i+1,:); 
    normOCRsim5_min = minsim5(i+1,:); 
    
%%%%%%%%%%%%%%%%%%%%%%
%%%% Data and Scaling
%%%%%%%%%%%%%%%%%%%%%%

    if i==0
        Data = [BartData.Set1.time BartData.Set1.BasDMSO BartData.Set1.BasDMSOstd];
        Data = Data(2:end-1,:);
    elseif  i==1
        Data = [BartData.Set1.time BartData.Set1.BMP4DMSO BartData.Set1.BMP4DMSOstd];
        Data = Data(2:end-1,:);
        Data2 = [BartData4.time BartData4.Dmso' BartData4.Dmsosem'];
        Data3 = [BartData4.time BartData4.DmsoAtra' BartData4.DmsoAtrasem'];
    elseif i==2
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
        
%%%%%%%%%%%%%
%%%% Plotting
%%%%%%%%%%%%%

        
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
            marker='^';
        elseif i >= 2 && i<3
            color=redDMSO;
            marker='^';
        elseif i >=3 && i<4
            color=greenISO;
            marker='s';
        elseif i>=4 && i<5
            color=brownISO;
            marker='v';
        else
            color = redISO;
             marker='v';
        end
        
        plottime = [0:time(end)];
        plottime2 = [0:time2(end)];
    
        %%%%%%%%%%%%%%% plot OCR
        if i==2 || i==5          
        figure(1)
        subplot(1,2,1)
        area1=[1:72,fliplr(1:72)];
        area2=[normOCRsim4_max(1:72),fliplr(normOCRsim4_min(1:72))];
        f=fill(area1,area2,[0.9 0.9 0.9]);
        set(f,'EdgeColor','none')
        hold on
        subplot(1,2,2)
        area1=[1:72,fliplr(1:72)];
        area2=[normOCRsim5_max(1:72),fliplr(normOCRsim5_min(1:72))];
        f=fill(area1,area2,[0.9 0.9 0.9]);
        set(f,'EdgeColor','none')
        hold on
        end
        
        if i==0 || i==3 || i==2 || i==5 
        figure(6)
        area1=[1:72,fliplr(1:72)];
        area2=[normOCRsim_max(1:72),fliplr(normOCRsim_min(1:72))];
        f=fill(area1,area2,[0.9 0.9 0.9]);
        set(f,'EdgeColor','none')
        hold on
        end
        
        if i==2 || i==5 
        figure(3)
        area1=[1:55,fliplr(1:55)];
        area2=[normOCRsim2_max(1:55),fliplr(normOCRsim2_min(1:55))];
        f=fill(area1,area2,[0.9 0.9 0.9]);
        set(f,'EdgeColor','none')
        hold on
        
        figure(5)
        area1=[1:55,fliplr(1:55)];
        area2=[normOCRsim3_max(1:55),fliplr(normOCRsim3_min(1:55))];
        f=fill(area1,area2,[0.9 0.9 0.9]);
        set(f,'EdgeColor','none')
        hold on
        end
        
        if i==1 || i==4
        figure(4)
        area1=[1:55,fliplr(1:55)];
        area2=[normOCRsim2_max(1:55),fliplr(normOCRsim2_min(1:55))];
        f=fill(area1,area2,[0.9 0.9 0.9]);
        set(f,'EdgeColor','none')
        hold on
        
        figure(7)
        area1=[1:72,fliplr(1:72)];
        area2=[normOCRsim_max(1:72),fliplr(normOCRsim_min(1:72))];
        f=fill(area1,area2,[0.9 0.9 0.9]);
        set(f,'EdgeColor','none')
        hold on
        
        figure(8)
        area1=[1:55,fliplr(1:55)];
        area2=[normOCRsim3_max(1:55),fliplr(normOCRsim3_min(1:55))];
        f=fill(area1,area2,[0.9 0.9 0.9]);
        set(f,'EdgeColor','none')
        hold on
        end        
end
end

