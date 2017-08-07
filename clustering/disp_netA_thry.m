%script to display/calc application of network architecture theory, Fig 4D

freq=[5;120]; 

open ../Data_std.fig

data_stdRt=[std(mean(fr_5Hz_normal,2)); std(mean(fr_120Hz_normal,2))];
figure
hold on
plot(freq,data_stdRt,'k.-','LineWidth',2,'MarkerSize',24)

load d_F_1_int_clus %load model results, network 1
plot(frq,std(fRt_vc),'.-','color',[0 .5 0],'MarkerSize',24,'LineWidth',2)


load d_F_2_int_wedge %load model results, network 2
plot(frq,std(fRt_vc),'b.-','MarkerSize',24,'LineWidth',2)
set(gca,'XLim',[0 125])
set(gca,'FontSize',18)
