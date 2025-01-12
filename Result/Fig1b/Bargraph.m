clc;clear;close all;
x = 1:3;
load("ReDate_OFDM_500.mat"); ReDate_500 = MEA_nmse; load("ReDate_OFDM_560.mat"); ReDate_560 = MEA_nmse; load("ReDate_OFDM_600.mat"); ReDate_600 = MEA_nmse;
load("MeDate_OFDM_500.mat"); MeDate_500 = MEA_nmse; load("MeDate_OFDM_560.mat"); MeDate_560 = MEA_nmse; load("MeDate_OFDM_600.mat"); MeDate_600 = MEA_nmse;

load("ReDate_FBMC_500.mat"); ReDate_FBMC_500 = Men; load("ReDate_FBMC_560.mat"); ReDate_FBMC_560 = Men; load("ReDate_FBMC_600.mat"); ReDate_FBMC_600 = Men;
load("MeDate_FBMC_500.mat"); MeDate_FBMC_500 = Men; load("MeDate_FBMC_560.mat"); MeDate_FBMC_560 = Men; load("MeDate_FBMC_600.mat"); MeDate_FBMC_600 = Men;

err_500 = abs(ReDate_500-MeDate_500); err_560 = abs(ReDate_560-MeDate_560); err_600 = abs(ReDate_600-MeDate_600);
err_FBMC_500 = abs(ReDate_FBMC_500-MeDate_FBMC_500); err_FBMC_560 = abs(ReDate_FBMC_560-MeDate_FBMC_560); err_FBMC_600 = abs(ReDate_FBMC_600-MeDate_FBMC_600);

Date_OFDM = [abs(ReDate_500),abs(MeDate_500),err_500;
        abs(ReDate_560),abs(MeDate_560),err_560;
        abs(ReDate_600),abs(MeDate_600),err_600];
Date_FBMC = [abs(ReDate_FBMC_500),abs(MeDate_FBMC_500),err_FBMC_500;
        abs(ReDate_FBMC_560),abs(MeDate_FBMC_560),err_FBMC_560;
        abs(ReDate_FBMC_600),abs(MeDate_FBMC_600),err_FBMC_600];
errdown_OFDM = Date_OFDM/50;
errup_OFDM   = Date_OFDM/50;
errdown_FBMC = Date_FBMC/50;
errup_FBMC   = Date_FBMC/50;
C1 = [180 180 180]./255;
C2 = [255 255 255]./255;
C3 = [255 110 220]./255;
figure('Color',[1 1 1]);
subplot(1,2,1)
GO_OFDM = bar(Date_OFDM,1,'EdgeColor','k');
GO_OFDM(1).FaceColor = C1;
GO_OFDM(2).FaceColor = C2;
GO_OFDM(2).FaceColor = C2;
hold on;
errorbar([1 2 3],Date_OFDM(:,2),errdown_OFDM(:,2),errup_OFDM(:,2),'k','LineStyle','None','LineWidth',0.5);
%errorbar([1.225 2.225 3.225],Date_OFDM(:,3),errdown_OFDM(:,3),errup_OFDM(:,3),'k','LineStyle','None','LineWidth',0.5);
errorbar([0.775 1.775 2.775],Date_OFDM(:,1),errdown_OFDM(:,1),errup_OFDM(:,1),'k','LineStyle','None','LineWidth',0.5);
set(gca,'box','off');
set(gca,'XGrid','off','YGrid','on');
set(gca,'TickDir','out','TickLength',[0.01,0.01],'XMinorTick','off','YMinorTick','off')
set(gca,'Xticklabel',{'500km/h','560km/h','600km/h'});
ylabel('NMSE ($-$dB) ','Interpreter','latex');
legend([GO_OFDM(1),GO_OFDM(2),GO_OFDM(3)],'Real-${\bf R}_{\bf H}$','Measured-${\bf R}_{\bf H}$','Err','FontSize',12,'Fontname','Times new roman ','Interpreter','latex');
title('OFDM');
set(gca,'FontSize',15,'Fontname','Times new roman ');

subplot(1,2,2)
GO_FBMC = bar(Date_FBMC,1,'EdgeColor','k');
GO_FBMC(1).FaceColor = C1;
GO_FBMC(2).FaceColor = C2;
GO_FBMC(2).FaceColor = C2;
hold on;
errorbar([1 2 3],Date_FBMC(:,2),errdown_FBMC(:,2),errup_FBMC(:,2),'k','LineStyle','None','LineWidth',0.5);
%errorbar([1.225 2.225 3.225],Date_FBMC(:,3),errdown_FBMC(:,3),errup_FBMC(:,3),'k','LineStyle','None','LineWidth',0.5);
errorbar([0.775 1.775 2.775],Date_FBMC(:,1),errdown_FBMC(:,1),errup_FBMC(:,1),'k','LineStyle','None','LineWidth',0.5);
set(gca,'box','off');
set(gca,'XGrid','off','YGrid','on');
set(gca,'TickDir','out','TickLength',[0.01,0.01],'XMinorTick','off','YMinorTick','off')
set(gca,'Xticklabel',{'500km/h','560km/h','600km/h'});
ylabel('NMSE ($-$dB)','Interpreter','latex');
legend([GO_FBMC(1),GO_FBMC(2),GO_FBMC(3)],'Real-${\bf R}_{\bf H}$','Measured-${\bf R}_{\bf H}$','Err','FontSize',12,'Fontname','Times new roman ','Interpreter','latex');
title('FBMC');
set(gca,'FontSize',15,'Fontname','Times new roman ');



