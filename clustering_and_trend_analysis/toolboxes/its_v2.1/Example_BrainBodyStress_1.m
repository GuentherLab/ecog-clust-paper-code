%% LINEAR ANALYSIS OF BRAIN-BODY INTERACTIONS IN REST/STRESS STATES
%% MUTUAL INFORMATION AND TRANSFER ENTROPY between BODY and BRAIN
clear; close all; clc;

P=3; %RR, resp, PAT
Q=4; %delta, theta, alpha, beta

num_el=14;
num_subj = 18;
pfilter=0.94; % AR detrending

p=3; %model order for information transfer

%% Load data
load('BrainBodyStress.mat');

%% ANALYSIS
for i_cond=1:3
    data_cond= series(:,i_cond);
            
    for i_subj=1:num_subj
        for i_el=1:num_el
            for i_target=1:7 % internal cycle: data preparation    
                % arange data for the given brain electrode
                data_body = data_cond{i_subj}(:,1:3);
                elIn = 4+(4*(i_el-1));
                elFin = elIn+3;
                data_brain = data_cond{i_subj}(:,elIn:elFin);
                data = [data_body,data_brain];
                [N,M]=size(data); assert(M==P+Q);
                % preprocessing: filter, normalize
                dataf=zeros(N,M);
                for m=1:M
                    dataf(:,m)=AR_filter(data,m,pfilter);
                    Z(:,m)=dataf(:,m)-mean(dataf(:,m));
                    Z(:,m)=(dataf(:,m)-mean(dataf(:,m)))/std(dataf(:,m));
                end
            end
            
            %% mutual information between districts
            ii=1:P;
            jj=P+1:M;
            [MI_XY,pXY]=its_MIlin_V(Z,ii,jj);
            
            results_MI(i_subj,i_el,i_cond)=MI_XY;
            results_pvalMI(i_subj,i_el,i_cond)=pXY;
            
            %% information transfer between districts
            out=its_BTElin_V(Z,ii,jj,p); % heart->brain
            results_TExy(i_subj,i_el,i_cond)=out.TE;
            results_pvalTExy(i_subj,i_el,i_cond)=out.pval_TE;
            
            out=its_BTElin_V(Z,jj,ii,p); % brain->heart
            results_TEyx(i_subj,i_el,i_cond)=out.TE;
            results_pvalTEyx(i_subj,i_el,i_cond)=out.pval_TE;
            
            
            end
        end
end

MIavg=squeeze(mean(results_MI,1));         
for i_el=1:num_el
    for i_cond=1:3
        num_pval(i_el,i_cond)=sum(results_pvalMI(:,i_el,i_cond)<0.05);
    end
end

TExyavg=squeeze(mean(results_TExy,1));
TEyxavg=squeeze(mean(results_TEyx,1));
for i_el=1:num_el
    for i_cond=1:3
        num_pvalxy(i_el,i_cond)=sum(results_pvalTExy(:,i_el,i_cond)<0.05);
        num_pvalyx(i_el,i_cond)=sum(results_pvalTEyx(:,i_el,i_cond)<0.05);
    end
end


%% Electrode locations
%j=1 ==> AF3;
label_Elettrodo(1)={'01 - AF3'}; label_ElettrodoN(1)={'AF3'};
%j=2 ==> F7;
label_Elettrodo(2)={'02 - F7'}; label_ElettrodoN(2)={'F7'};
%j=3 ==> F3;
label_Elettrodo(3)={'03 - F3'}; label_ElettrodoN(3)={'F3'};
%j=4 ==> FC5;
label_Elettrodo(4)={'04 - FC5'}; label_ElettrodoN(4)={'FC5'};
%j=5 ==> T7;
label_Elettrodo(5)={'05 - T7'}; label_ElettrodoN(5)={'T7'};
%j=6 ==> P7;
label_Elettrodo(6)={'06 - P7'}; label_ElettrodoN(6)={'P7'};
%j=7 ==> O1;
label_Elettrodo(7)={'07 - O1'}; label_ElettrodoN(7)={'O1'};
%j=8 ==> O2;
label_Elettrodo(8)={'08 - O2'}; label_ElettrodoN(8)={'O2'};
%j=9 ==> P8;
label_Elettrodo(9)={'09 - P8'}; label_ElettrodoN(9)={'P8'};
%j=10 ==> T8;
label_Elettrodo(10)={'10 - T8'}; label_ElettrodoN(10)={'T8'};
%j=11 ==> FC6;
label_Elettrodo(11)={'11 - FC6'}; label_ElettrodoN(11)={'FC6'}; 
%j=12 ==> F4;
label_Elettrodo(12)={'12 - F4'}; label_ElettrodoN(12)={'F4'};
%j=13 ==> F8;
label_Elettrodo(13)={'13 - F8'}; label_ElettrodoN(13)={'F8'};
%j=14 ==> AF4
label_Elettrodo(14)={'14 - AF4'}; label_ElettrodoN(14)={'AF4'};

%posElettrodi
posX(1) = 3.0;posY(1) = 9.1;
posX(2) = 1.9;posY(2) = 7.6;
posX(3) = 3.8;posY(3) = 7.3;
posX(4) = 2.7;posY(4) = 5.5;
posX(5) = 0.7;posY(5) = 5.6;
posX(6) = 1.5;posY(6) = 3.7;
posX(7) = 3.5;posY(7) = 0.7;
posX(8) = 6.5;posY(8) = 0.7;
posX(9) = 8.5;posY(9) = 3.7;
posX(10) = 9.3;posY(10) = 5.6;
posX(11) = 7.5;posY(11) = 5.5;
posX(12) = 6.5;posY(12) = 7.3;
posX(13) = 8.1;posY(13) = 7.6;
posX(14) = 7.0;posY(14) = 9.1;

posXn =posX./10;
posYn =posY./10;

K=100;
[Xi,Yi]=meshgrid(1:1:K,1:K);%crea la griglia
th = 0:pi/50:2*pi; x= K/2; y= K/2;
xunit = K/2 * cos(th) + x;
yunit = K/2 * sin(th) + y;


%% FIGURES MI (average over subjects + n. of subjects with significant MI)
condlabel{1}='REST'; condlabel{2}='MENTAL'; condlabel{3}='GAME';
colormap_minimum = 1.1*min(MIavg(:));
colormap_maximum = 0.9*max(MIavg(:));
for i_cond=1:3
    figure(1);subplot(2,3,i_cond); hold on; colormap(gca,'jet')
    mappa=griddata(posXn*K,posYn*K,MIavg(:,i_cond),Xi,Yi,'v4');%interpolazione
    imagesc([0 K],[0 K],mappa); colorbar
    caxis([colormap_minimum colormap_maximum]);
    pbaspect([1 1 1]);
    set(gca,'FontSize',12,'YTickLabel',[],'XTickLabel',[]);
    set(gca,'YDir','normal');
    for j=1:num_el
        plot(posXn(j)*K,posYn(j)*K,'ok','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',6);
    end
    xlim([0 K]); ylim([0 K]);
    plot(xunit, yunit,'--k','LineWidth',3);
	title([condlabel{i_cond} ',I_{XY}'])
    
    subplot(2,3,i_cond+3); hold on; colormap(gca,'parula')
    mappa=griddata(posXn*K,posYn*K,num_pval(:,i_cond),Xi,Yi,'v4');%interpolazione
    imagesc([0 K],[0 K],mappa); colorbar
    caxis([0 num_subj]);
    pbaspect([1 1 1]);
    set(gca,'FontSize',12,'YTickLabel',[],'XTickLabel',[]);
    set(gca,'YDir','normal');
    for j=1:num_el
        plot(posXn(j)*K,posYn(j)*K,'ok','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',6);
    end
    xlim([0 K]); ylim([0 K]);
    plot(xunit, yunit,'--k','LineWidth',3);
	title([condlabel{i_cond} ',n.of subjects'])
    
end


%% FIGURES TE (average over subjects)
colormap_minimum = 1*min(min(TExyavg(:)),min(TEyxavg(:)));
colormap_maximum = 1*max(max(TExyavg(:)),max(TEyxavg(:)));
for i_cond=1:3
    figure(2);subplot(2,3,i_cond); hold on; colormap(jet)
    mappa=griddata(posXn*K,posYn*K,TExyavg(:,i_cond),Xi,Yi,'v4');%interpolazione
    imagesc([0 K],[0 K],mappa); colorbar
    caxis([colormap_minimum colormap_maximum]);
    pbaspect([1 1 1]);
    set(gca,'FontSize',12,'YTickLabel',[],'XTickLabel',[]);
    set(gca,'YDir','normal');
    for j=1:num_el
        plot(posXn(j)*K,posYn(j)*K,'ok','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',6);
    end
    xlim([0 K]); ylim([0 K]);
    plot(xunit, yunit,'--k','LineWidth',3);
	title([condlabel{i_cond} ',T_{X \rightarrow Y}'])    
    
    figure(2);subplot(2,3,i_cond+3); hold on; colormap(jet)
    mappa=griddata(posXn*K,posYn*K,TEyxavg(:,i_cond),Xi,Yi,'v4');%interpolazione
    imagesc([0 K],[0 K],mappa); colorbar
    caxis([colormap_minimum colormap_maximum]);
    pbaspect([1 1 1]);
    set(gca,'FontSize',12,'YTickLabel',[],'XTickLabel',[]);
    set(gca,'YDir','normal');
    for j=1:num_el
        plot(posXn(j)*K,posYn(j)*K,'ok','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',6);
    end
    xlim([0 K]); ylim([0 K]);
    plot(xunit, yunit,'--k','LineWidth',3);
	title([condlabel{i_cond} ',T_{Y \rightarrow X}']) 
end




   