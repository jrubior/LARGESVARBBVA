% this file plots the figures
clear variables;
currdir = pwd;
restoredefaultpath;
% replicate Figure 1 in KM (2012)

fixed_rf=0;
label_R = 'kilianmurphy';
prior_type = 'minnesota';
ndraws  =1e4;
prior_only=0;
spec= ['rgibbs_',label_R,'prior_only_',num2str(prior_only),'prior_',prior_type,'fixed_rf_',num2str(fixed_rf),'_ndraws',num2str(ndraws)];




L =load(['matfiles/',filesep,spec,'.mat']);




nhorizons = size(L.L,1);
nvar      = size(L.L,2);
nshocks   = size(L.L,3);



L_q50=nan(nhorizons,nvar,nshocks); % store conjF quantile 50th
L_q16=nan(nhorizons,nvar,nshocks); % store conjF quantile 16th
L_q84=nan(nhorizons,nvar,nshocks); % store conjF quantile 84th

Lcum_q50=nan(nhorizons,nvar,nshocks); % store conjF quantile 50th
Lcum_q16=nan(nhorizons,nvar,nshocks); % store conjF quantile 16th
Lcum_q84=nan(nhorizons,nvar,nshocks); % store conjF quantile 84th


for ii=1:nhorizons
    for jj=1:nvar
        for kk=1:nshocks
            

        L_q50(ii,jj,kk) = quantile(L.L(ii,jj,kk,:),0.5);
        L_q16(ii,jj,kk) = quantile(L.L(ii,jj,kk,:),0.16);
        L_q84(ii,jj,kk) = quantile(L.L(ii,jj,kk,:),0.84);

        Lcum_q50(ii,jj,kk) = quantile(L.cumL(ii,jj,kk,:),0.5);
        Lcum_q16(ii,jj,kk) = quantile(L.cumL(ii,jj,kk,:),0.16);
        Lcum_q84(ii,jj,kk) = quantile(L.cumL(ii,jj,kk,:),0.84);

        end
    end
end




% figure settings
ftsizeaxis   = 11;
ftsizexlabel = 11;
ftsizetitle  = 11;
ftsizelegend = 11;
ftlinewidth  = 1.0;
medianwidth  = 1.0;


plot_uv_scale = 1;
horizon       = L.horizon;


%% figure 1

close all

hFig = figure('name','output');
set(hFig, 'Position', [20 20 900 900])

subplot(3,4,1)
variable   = 1;
shock      = 1;
auv16=Lcum_q16(:,variable,shock)'*plot_uv_scale;
auv50=Lcum_q50(:,variable,shock)'*plot_uv_scale;
buv84=Lcum_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Oil Production (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Flow Supply Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-2 1])



subplot(3,4,2)
variable   = 2;
shock      = 1;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Real Activity (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Flow Supply Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-5 10])

subplot(3,4,3)
variable   = 3;
shock      = 1;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Real Price of Oil (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Flow Supply Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-5 10])


subplot(3,4,4)
variable   = 4;
shock      = 1;
auv16=Lcum_q16(:,variable,shock)'*plot_uv_scale;
auv50=Lcum_q50(:,variable,shock)'*plot_uv_scale;
buv84=Lcum_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Inventories (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Flow Supply Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-20 20])






subplot(3,4,5)
variable   = 1;
shock      = 2;
auv16=Lcum_q16(:,variable,shock)'*plot_uv_scale;
auv50=Lcum_q50(:,variable,shock)'*plot_uv_scale;
buv84=Lcum_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Oil Production (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Flow Demand Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-1 2])



subplot(3,4,6)
variable   = 2;
shock      = 2;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Real Activity (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Flow Demand Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-5 10])

subplot(3,4,7)
variable   = 3;
shock      = 2;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Real Price of Oil (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Flow Demand Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-5 10])


subplot(3,4,8)
variable   = 4;
shock      = 2;
auv16=Lcum_q16(:,variable,shock)'*plot_uv_scale;
auv50=Lcum_q50(:,variable,shock)'*plot_uv_scale;
buv84=Lcum_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Inventories (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Flow Demand Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-20 20])



subplot(3,4,9)
variable   = 1;
shock      = 3;
auv16=Lcum_q16(:,variable,shock)'*plot_uv_scale;
auv50=Lcum_q50(:,variable,shock)'*plot_uv_scale;
buv84=Lcum_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Oil Production (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Speculative Demand Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-1 2])



subplot(3,4,10)
variable   = 2;
shock      = 3;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Real Activity (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Speculative Demand Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-5 10])

subplot(3,4,11)
variable   = 3;
shock      = 3;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Real Price of Oil (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Speculative Demand Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-5 10])


subplot(3,4,12)
variable   = 4;
shock      = 3;
auv16=Lcum_q16(:,variable,shock)'*plot_uv_scale;
auv50=Lcum_q50(:,variable,shock)'*plot_uv_scale;
buv84=Lcum_q84(:,variable,shock)'*plot_uv_scale;
hold on
hq16=plot([0:horizon],auv16,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hq84=plot([0:horizon],buv84,'linestyle','--','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','--','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Inventories (%)','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title('Speculative Demand Shock','FontSize',ftsizetitle)
box off
legend off
grid on
ylim([-20 20])


if ~exist('pngfiles', 'dir')
    mkdir('pngfiles');
    disp('Created directory: pngfiles');
else
    disp('Directory pngfiles already exists');
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'.png'],'-dpng');
