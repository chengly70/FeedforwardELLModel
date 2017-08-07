%script to plot effective corr vs std rates (Fig 4E)

fr_min=0;
fr_max=125; %in Hz
d_frq=2.5; %mesh incr of frq discr (Hz)
frq = (fr_min+d_frq/2: d_frq: fr_max-d_frq/2)'; %mesh for frq vals in connMat J_di

% Parms (conn, q,th, presy_Frq_tuned, Ne/ei, etc.), from
% gen_hetConn_Parms.m, [ftuning_aff.m, type files]
load Parms_clust

[frq_v,eff_corr_Net1]=calc_effCorrDrive(qPreF,ThresE,presy_Frq_tuned,J_di,frq);


figure
hold on
load d_F_1clus %load model results, network 1
plot(eff_corr_Net1,std(fRt_vc),'.-','color',[0 .5 0],'MarkerSize',24,'LineWidth',2)

load Parms_wedge
[frq_v,eff_corr_Net2]=calc_effCorrDrive(qPreF,ThresE,presy_Frq_tuned,J_di,frq);

load d_F_2wedge %load model results, network 1
plot(eff_corr_Net2,std(fRt_vc),'b.-','MarkerSize',24,'LineWidth',2)
set(gca,'XLim',[-1 1]) %manually reverse axis
set(gca,'FontSize',18)