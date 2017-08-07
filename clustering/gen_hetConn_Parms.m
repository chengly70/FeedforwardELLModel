%Generate the uncorr qPreF,ThresE, coupling matrix, presy_Frq_tuned, Ne/ei
% comment out the appropriate lines on lines 24-28 
% (ftuning_aff<->Parms_clust_int) and (ftuning_wedge<->Parms_wedge_int)

Nei=100;  %size of each E (granule) and I (interneuron) presyn

Ne=1000; %# E-cells of interest

fr_min=0;
fr_max=125; %in Hz
d_frq=2.5; %mesh incr of frq discr (Hz)
frq = (fr_min+d_frq: d_frq: fr_max-d_frq)'; %mesh for frq vals in connMat J_di

%altering intrins & net. heterog. relationship
perc_th=(0.05: 0.9/(Ne-1): 0.95)';
sig_nm=.2;
qPreF=(0.5: 1/(Ne-1): 1.5)';
ThresE0=exp(norminv(perc_th,-sig_nm^2/2,sig_nm)); %log-normal with mean exp(0)=1
load ../storedPerm
Thres1=ThresE0(strdPerm); %crtW_corr function works better with random samples of same vector
Thres1=Thres1(randperm(Ne)); %Randomize order to get diff (qPreF,ThresE) combo.
ThresE=crtW_corr(qPreF,Thres1,0); %no correl between (q,th); can chance

%[presy_Frq_tuned,J_di,W_di]=ftuning_aff(qPreF,ThresE,Ne,Nei);
[presy_Frq_tuned,J_di,W_di]=ftuning_wedge(qPreF,ThresE,Ne,Nei); %line54; use 100(1), 70(2)

%save Parms_clust_int qPreF ThresE sig_nm perc_th Nei Ne W_di J_di presy_Frq_tuned
save Parms_wedge_int qPreF ThresE sig_nm perc_th Nei Ne W_di J_di presy_Frq_tuned