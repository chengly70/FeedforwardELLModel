%script to display/calc application of theory using specific Corr

freq=[5;120]; 

open Data_std.fig %contains std of Data

load d_Sm_int_July17 %load model results
[m,ind5]=min(abs(std(mean(fr_5Hz_normal,2))-std(fRt_vc)));
corr_preTh(ind5)

mod_stdRt(1,1)=std(fRt_vc(:,ind5));

load d_Bg_int_July17 %load model results
[m,ind120]=min(abs(std(mean(fr_120Hz_normal,2))-std(fRt_vc)));
corr_preTh(ind120)

mod_stdRt(2,1)=std(fRt_vc(:,ind120));
plot(freq,mod_stdRt,'c.-','MarkerSize',24)
set(gca,'FontSize',18)