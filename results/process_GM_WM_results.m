%% define path
%% 3_1 results
clearvars
fig=false;
fig2=true;
for iii=1:5
    switch iii
        case 1
            misspool='no amine';
        case 2
            misspool='no NOE 1p7';
        case 3
            misspool='no NOE 3p5';
        case 4
            misspool='no amide';
        case 5
            misspool='all pools';
    end
            
    cd(['FTEST_3_' num2str(iii)])
    
    ALLresults=nan(12,23,100);
    %% over all runs
    for r=1:100
        try        ALLresults(1,:,r)=csvread(['results1kp100R_',num2str(r),'Sub01_GM']);catch;end
        try        ALLresults(2,:,r)=csvread(['results1kp100R_',num2str(r),'Sub02_GM']);catch;end
        try        ALLresults(3,:,r)=csvread(['results1kp100R_',num2str(r),'Sub03_GM']);catch;end
        try        ALLresults(4,:,r)=csvread(['results1kp100R_',num2str(r),'Sub04_GM']);catch;end
        try        ALLresults(5,:,r)=csvread(['results1kp100R_',num2str(r),'Sub05_GM']);catch;end
        try        ALLresults(6,:,r)=csvread(['results1kp100R_',num2str(r),'Sub06_GM']);catch;end
        
        try        ALLresults(7,:,r)=csvread(['results1kp100R_',num2str(r),'Sub01_WM']);catch;end
        try        ALLresults(8,:,r)=csvread(['results1kp100R_',num2str(r),'Sub02_WM']);catch;end
        try        ALLresults(9,:,r)=csvread(['results1kp100R_',num2str(r),'Sub03_WM']);catch;end
        try        ALLresults(10,:,r)=csvread(['results1kp100R_',num2str(r),'Sub04_WM']);catch;end
        try        ALLresults(11,:,r)=csvread(['results1kp100R_',num2str(r),'Sub05_WM']);catch;end
        try        ALLresults(12,:,r)=csvread(['results1kp100R_',num2str(r),'Sub06_WM']);catch;end
        
    end
    
    
    %% find best run
    
    SSq = squeeze(ALLresults(:,23,:)); %12rows x 100cols
    for s=1:12
        [~,bestspec(s)] = nanmin(SSq(s,:));
        RR(s,:) = ALLresults(s,:,bestspec(s));
    end
    
    GMSSq = RR(1:6,23);
    WMSSq = RR(7:12,23);
    
    
    %% sort result - 12rows x 23 cols
    
    results=RR;
    results(:,1) = RR(:,1)*100;
    results(:,5) = RR(:,5)*100;
    results(:,9) = RR(:,9)*100;
    results(:,13) = RR(:,13)*100;
    results(:,17) = RR(:,17)*100;
    
    results(:,3) = RR(:,3)*1000;
    results(:,7) = RR(:,7)*1000;
    results(:,11) = RR(:,11)*1000;
    results(:,15) = RR(:,15)*1000;
    results(:,19) = RR(:,19)*1000;
    
    %close all
    
    
    %% graphs
    if fig
        cmap=[jet(99); 0,0,0];
        
        
        %% hists
        
        xplot = [1:6,9:14];
        
        for pool=1:5
            for variable=1:3
                
                figure,hold on
                plot(xplot(1:6),results(1:6,4*(pool-1)+variable),'bx','LineWidth',2,'MarkerSize',18)
                plot(xplot(7:12),results(7:12,4*(pool-1)+variable),'mx','LineWidth',2,'MarkerSize',18)
                xlim([0 15]),xlabel('GM            Subject            WM')
                set(gca,'XTick',xplot,'XTickLabel',[1:6,1:6])
                
                if pool==1
                    prefix='MT ';
                elseif pool==2
                    prefix='Amide ';
                elseif pool==3
                    prefix='NOE ';
                elseif pool==4
                    prefix='NOE (-1.7ppm) ';
                elseif pool==5
                    prefix='Creatine ';
                end
                if variable==1
                    suffix='pool size (%)';
                elseif variable==2
                    suffix='exchange rate (Hz)';
                elseif variable==3
                    suffix='T_2 (ms)';
                end
                ylabel([prefix suffix])
            end
        end
        
    end
    %% look
    
    
    load('freqs.mat')
    freqs=sort(freqs(2:64));
    
    data(1,:)=csvread('spec_Sub01_GM');
    data(2,:)=csvread('spec_Sub02_GM');
    data(3,:)=csvread('spec_Sub03_GM');
    data(4,:)=csvread('spec_Sub04_GM');
    data(5,:)=csvread('spec_Sub05_GM');
    data(6,:)=csvread('spec_Sub06_GM');
    data(7,:)=csvread('spec_Sub01_WM');
    data(8,:)=csvread('spec_Sub02_WM');
    data(9,:)=csvread('spec_Sub03_WM');
    data(10,:)=csvread('spec_Sub04_WM');
    data(11,:)=csvread('spec_Sub05_WM');
    data(12,:)=csvread('spec_Sub06_WM');
    
    B1inh=[0.616480,0.616963,0.605006,0.674842,0.646471,0.560476,...
        0.629720,0.625226,0.610191,0.696420,0.658614,0.574476];
    
    if fig2
        %% paper figure
        for d=1:12
            
            if d==1
                figure('units','normalized','outerposition',[0 0 1 1]),hold on
                col=winter(5);
                titstr=['GM spectra for 5 pools, ' misspool];
                d2=0;
            elseif d==7
                figure('units','normalized','outerposition',[0 0 1 1]),hold on
                d2=6;
                titstr=['WM spectra for 5 pools, ' misspool];
                col=[1 0 1;1 0.125 0.75;1 0.25 0.5;1 0.375 0.25;1 0.5 0];
            end
            
            subplot(2,3,d-d2);hold on
            plot(freqs/300,data(d,1:63),'ko')
            plot(freqs/300,data(d,64:126),'ko')
            plot(freqs/300,data(d,127:189),'ko')
            plot(freqs/300,data(d,190:252),'ko')
            plot(freqs/300,data(d,253:315),'ko')
            
            sim=BMsim(RR(d,:),[0.33 0.67 1 1.33 1.67]*B1inh(d),freqs,3);
            for p=1:5
                plot(freqs/300,sim(:,p),'Color',col(p,:),'LineWidth',1)
            end
            xlabel('ppm'),ylabel('S/S_0')
            xlim([-6 6])
            set(gca,'XDir','reverse')
            if (d==6 || d==12)
                suptitle(titstr)
                saveas(gcf,['D:\Code\PSO\FTESTS\' titstr '.png'])
            end
        end
        %% end paper figure
    end
    
    %% calc averages $\pm$ variance
    
    for p=1:5
        GM_M0_mean(p) = mean(results(1:6,1+(p-1)*4));
        GM_M0_std(p) = std(results(1:6,1+(p-1)*4));
        GM_k_mean(p) = mean(results(1:6,2+(p-1)*4));
        GM_k_std(p) = std(results(1:6,2+(p-1)*4));
        GM_T2_mean(p) = mean(results(1:6,3+(p-1)*4));
        GM_T2_std(p) = std(results(1:6,3+(p-1)*4));
        
        WM_M0_mean(p) = mean(results(7:12,1+(p-1)*4));
        WM_M0_std(p) = std(results(7:12,1+(p-1)*4));
        WM_k_mean(p) = mean(results(7:12,2+(p-1)*4));
        WM_k_std(p) = std(results(7:12,2+(p-1)*4));
        WM_T2_mean(p) = mean(results(7:12,3+(p-1)*4));
        WM_T2_std(p) = std(results(7:12,3+(p-1)*4));
    end
    
    %% print table
    
    TABLE = {[num2str(GM_M0_mean(1),'%.2f'),' ',char(177),' ',num2str(GM_M0_std(1),'%.2f')],...
        [num2str(GM_M0_mean(2),'%.2f'),' ',char(177),' ',num2str(GM_M0_std(2),'%.2f')],...
        [num2str(GM_M0_mean(3),'%.2f'),' ',char(177),' ',num2str(GM_M0_std(3),'%.2f')],...
        [num2str(GM_M0_mean(4),'%.2f'),' ',char(177),' ',num2str(GM_M0_std(4),'%.2f')],...
        [num2str(GM_M0_mean(5),'%.2f'),' ',char(177),' ',num2str(GM_M0_std(5),'%.2f')];...
        ...
        [num2str(WM_M0_mean(1),'%.2f'),' ',char(177),' ',num2str(WM_M0_std(1),'%.2f')],...
        [num2str(WM_M0_mean(2),'%.2f'),' ',char(177),' ',num2str(WM_M0_std(2),'%.2f')],...
        [num2str(WM_M0_mean(3),'%.2f'),' ',char(177),' ',num2str(WM_M0_std(3),'%.2f')],...
        [num2str(WM_M0_mean(4),'%.2f'),' ',char(177),' ',num2str(WM_M0_std(4),'%.2f')],...
        [num2str(WM_M0_mean(5),'%.2f'),' ',char(177),' ',num2str(WM_M0_std(5),'%.2f')];...
        ...
        [num2str(GM_k_mean(1),'%.2f'),' ',char(177),' ',num2str(GM_k_std(1),'%.2f')],...
        [num2str(GM_k_mean(2),'%.2f'),' ',char(177),' ',num2str(GM_k_std(2),'%.2f')],...
        [num2str(GM_k_mean(3),'%.2f'),' ',char(177),' ',num2str(GM_k_std(3),'%.2f')],...
        [num2str(GM_k_mean(4),'%.2f'),' ',char(177),' ',num2str(GM_k_std(4),'%.2f')],...
        [num2str(GM_k_mean(5),'%.2f'),' ',char(177),' ',num2str(GM_k_std(5),'%.2f')];...
        ...
        [num2str(WM_k_mean(1),'%.2f'),' ',char(177),' ',num2str(WM_k_std(1),'%.2f')],...
        [num2str(WM_k_mean(2),'%.2f'),' ',char(177),' ',num2str(WM_k_std(2),'%.2f')],...
        [num2str(WM_k_mean(3),'%.2f'),' ',char(177),' ',num2str(WM_k_std(3),'%.2f')],...
        [num2str(WM_k_mean(4),'%.2f'),' ',char(177),' ',num2str(WM_k_std(4),'%.2f')],...
        [num2str(WM_k_mean(5),'%.2f'),' ',char(177),' ',num2str(WM_k_std(5),'%.2f')];...
        ...
        [num2str(GM_T2_mean(1),'%.3f'),' ',char(177),' ',num2str(GM_T2_std(1),'%.3f')],...
        [num2str(GM_T2_mean(2),'%.3f'),' ',char(177),' ',num2str(GM_T2_std(2),'%.3f')],...
        [num2str(GM_T2_mean(3),'%.3f'),' ',char(177),' ',num2str(GM_T2_std(3),'%.3f')],...
        [num2str(GM_T2_mean(4),'%.3f'),' ',char(177),' ',num2str(GM_T2_std(4),'%.3f')],...
        [num2str(GM_T2_mean(5),'%.3f'),' ',char(177),' ',num2str(GM_T2_std(5),'%.3f')];...
        ...
        [num2str(WM_T2_mean(1),'%.3f'),' ',char(177),' ',num2str(WM_T2_std(1),'%.3f')],...
        [num2str(WM_T2_mean(2),'%.3f'),' ',char(177),' ',num2str(WM_T2_std(2),'%.3f')],...
        [num2str(WM_T2_mean(3),'%.3f'),' ',char(177),' ',num2str(WM_T2_std(3),'%.3f')],...
        [num2str(WM_T2_mean(4),'%.3f'),' ',char(177),' ',num2str(WM_T2_std(4),'%.3f')],...
        [num2str(WM_T2_mean(5),'%.3f'),' ',char(177),' ',num2str(WM_T2_std(5),'%.3f')]...
        };
    
    %%
    pwd
    %disp(TABLE)
    disp([GMSSq' ])
    disp([WMSSq'])
    
    %% t-test (2tailed)
    
    GMres = results(1:6,:);
    WMres = results(7:12,:);
    
    for vv=1:23
        [~,p(vv)]=ttest2(GMres(:,vv),WMres(:,vv));
    end
    
    cd ..
end



