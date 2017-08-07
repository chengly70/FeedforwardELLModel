%script to plot cycle histogram; match to fast freq stim

load dCyclHist_Bg

MnCyc_hist120=mean(fCyc_rt)';

open high_frq.fig
hold on
d_xc120=x_cyc(2)-x_cyc(1);

%b/c of delay, peak is at pi for Bg:
mx_pkMd=round((pi-x_cyc(1))/d_xc120)+1;
%index for peak of data
mx_data=round((4.9440-x_cyc(1))/d_xc120)+1;
%assuming shft_amt < 0
shft_amt=mx_pkMd-mx_data; %shifting all indices by this amount

nw_Indx=mod((1:length(x_cyc))'+shft_amt,length(x_cyc))+1;

plot(x_cyc,MnCyc_hist120(nw_Indx),'b','LineWidth',4)

% for the 5Hz freq
load dCyclHist_Sm

MnCyc_hist5=mean(fCyc)';

open low_frq.fig
hold on
d_xc5=x_cyc(2)-x_cyc(1);

%b/c of delay, peak is at 2.25 for Sm:
mx_pkMd=round((2.25-x_cyc(1))/d_xc5)+1;
%index for peak of data
mx_data=round((4.351-x_cyc(1))/d_xc5)+1;
%assuming shft_amt < 0
shft_amt=mx_pkMd-mx_data; %shifting all indices by this amount
nw_Indx=mod((1:length(x_cyc))'+shft_amt,length(x_cyc))+1;

plot(x_cyc,MnCyc_hist5(nw_Indx),'b','LineWidth',4)