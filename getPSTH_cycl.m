%script to get PSTH/Cycle histogram; specify corr(q,thres) in 'corrPrc'
% and freq in 'frq'; we only do this separately for large freq b/c the time step in 
%mx file is 1ms, too small for 120hz cycle histogram, so we use 0.1ms here; 
%for small freq, the PSTH/Cycle histogram is calc in scpt_SMvC_Aprl17.m via
%mx_tt_FFpyr.c

Nei=200;  %size of BOTH E and I presyn

Ne=1000; %# E-cells of interest

frq=120; %small or big chirp
corrPrc=0; %specify (q,Thres) correlation

dc_val=0.4;
Amp=0.4;
scl_Iinp=1.2;

% [W_ie,W_ee,J_ie,J_ee]=genConn(Ne,Nei,fr_ItoE,fr_EtoE,1);
% --- loading so have same connect matrix across sims ---
load conRnd_A
fr_ItoE=0.2; %make sure it matches what was loaded
fr_EtoE=0.2;
g_ie=1/(fr_ItoE*Nei);
g_ee=2.3/(fr_EtoE*Nei);
gieRaw=g_ie*fr_ItoE*Nei;
geeRaw=g_ee*fr_EtoE*Nei;

[I_aff,I_ie,I_ee,tV,nuI,nuE]=getSinInp(frq,dc_val,Amp,scl_Iinp,Nei);
%incorp. delay (20ms) in I_ie, I_ee 'feedback'
t_delay=.02; %in sec
dt=tV(2)-tV(1);
num_shift=round(t_delay/dt)+1; %tV(1)=dt, not 0
I_ieS=zeros(size(I_ie));
I_ieS(:,num_shift:end)=I_ie(:,1:end-num_shift+1);
I_eeS=zeros(size(I_ee));
I_eeS(:,num_shift:end)=I_ee(:,1:end-num_shift+1);

%altering intrins & net. heterog. relationship
perc_th=(0.05: 0.9/(Ne-1): 0.95)';
sig_nm=.1;
qPreF=(0.5: 1/(Ne-1): 1.5)';
ThresE0=exp(norminv(perc_th,-sig_nm^2/2,sig_nm)); %log-normal with mean exp(0)=1
load storedPerm
Thres1=ThresE0(strdPerm); %crtW_corr function works better with random samples of same vector

ThresE=crtW_corr(qPreF,Thres1,corrPrc);

tic
    %[nu_e,f_rate,f_isi,t_vc]=mx_tt_FFpyr(I_aff,I_ieS,I_eeS,tV,W_ee,g_ie,W_ee,g_ee,ThresE,qPreF);
    [ThresE,qPreF,nuI,nuE,f_rate,nu_e]=ff_LTsim(J_ie,J_ee,g_ie,g_ee,frq,dc_val,Amp,scl_Iinp,Nei,Ne,corrPrc);
toc

dts=0.0001; %in sec, same as dts in mx_tt_FFpyr.c
t_PSTH=(dts:dts:1)'; %time in each realiz of mx_tt_FFpyr.c
f_PSTH=(sum(nu_e)./Ne)';

nm_cyc=round(1/frq/dts);
fCyc_rt=nu_e(:,1:nm_cyc);
for j=2:frq %b/c 1 sec in time in mx_tt_FFpyr.c
	fCyc_rt = (fCyc_rt*(j-1)+nu_e(:,(j-1)*nm_cyc+1:j*nm_cyc))./j; %running sum
end
x_cyc=(0 : 2*pi/(nm_cyc-1) : 2*pi)'; %scale so between 0 & 2*pi


save dCyclHist_Bg t_PSTH f_PSTH x_cyc fCyc_rt