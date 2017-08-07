function [ThresE,qPreF,nuI,nuE,f_rate,nu_t]=ff_LTsim_BGcyc(J_ie,J_ee,g_ie,g_ee,frq,dc_val,Amp,Nei,Ne,corrPrc)
%Input: J_ie,J_ee are saved connect matrices from genConn.m, g_ie/g_ee are conduct vals, 
%cutoff_frq (Hz); dc_val (limit of direct input); amp_rng (amplitude
%range); 
%function to determine effective current input to E cell; 
%I_aff band-limited Random Amplitude Modulation (RAM)
%[units of norm. volt [0,1], assuming feedforward architechure;
%Nei is the number of of each presyn E and presyn I neurons. 
% addressing Reviewer 1's interneuron concern; no longer receives direct
% sinusoidal input, but we have new parameter s_EE
% function specifically for 120Hz stim Cycle Hist (need smaller time-step
% for visual, etc)


N_rlz = 42;  %number of realizations

perc_th=(0.05: 0.9/(Ne-1): 0.95)';

sig_nm=0.1;
qPreF=(0.5: 1/(Ne-1): 1.5)';
ThresE0=exp(norminv(perc_th,-sig_nm^2/2,sig_nm)); %log-normal with mean exp(0)=1
load storedPerm
Thres1=ThresE0(strdPerm); %crtW_corr function works better with random samples of same vector
ThresE=crtW_corr(qPreF,Thres1,corrPrc);

%Choice of t_end insures a whole number of cycles; only a problem when freq has decimal
t_end=1; %in sec; 
dt=.0001; %in sec
dts=.0001; %sub-sampling in sims 
tV=(dt:dt:t_end)';
L_t=length(tV);
spT=round(dts/dt);

%delay time from FF input
t_dely=0.02; %in sec
nm_gBck=floor(t_dely/dt); %number of time-steps go back for time-delay

%biophys params
%Nei=200; %numb of Inhib neurons
nuI=zeros(Nei,1); %firing rate of presy I-cells
nuE=zeros(Nei,1); %firing rate of presy E-cells
f_rate=zeros(Ne,1); %firing rate of Target
nu_t=zeros(Ne,round(t_end/dts)); %PSTH of Target
tau_m=0.01; %mem. time const. in sec
t_ref=0.0005; %ref period in sec
t_r=0.002; %Both synapses; rise time
t_d=0.01;  %Both synapses; decay time
ampS=2; %jump in asyn upon spiking
sig_a=1;
tau_eta=0.005;
sqtn=1/sqrt(dt*tau_eta);
%for target superf Pyram
sig_BP=0.75;
ref_s=.001;
Isyn=-.5;
Esyn=6.5;
vmin=-.5;


v_I=rand(Nei,1);%zeros(Nei,1);
i_syn=zeros(Nei,nm_gBck);  %track previous ones up to t_dely
i_asyn=zeros(Nei,1);
eta_i=zeros(Nei,1);
v_Ei=rand(Nei,1);
e_syn=zeros(Nei,nm_gBck);  %track previous ones up to t_dely
e_asyn=zeros(Nei,1);
eta_e=zeros(Nei,1);

v_Pyr=rand(Ne,1);
eta_Pyr=zeros(Ne,1);
Gee=zeros(Ne,1);
Gie=zeros(Ne,1);

%1 new param; interneurons
s_EE=2.7/Nei;

% Afferent input I_aff to E (& I) cells; frozen RAM
I_aff=dc_val+Amp*sin(2*pi*tV*frq);
I_aff(I_aff<0)=0; %threshold so f-rate incr

%pass in connect (random Erdos-Reyni),conducts passed in
%put in weights g_ in J_ matrices
J_ie=g_ie*J_ie;
J_ee=g_ee*J_ee;

nrnSpace=floor(t_ref/dt);
nrnSpaceP=floor(ref_s/dt);

for k=1:N_rlz
    if(k==1)
        TmIspk=-nrnSpace*ones(Nei,1);
        TmEspk=-nrnSpace*ones(Nei,1);
        TmPspk=-nrnSpaceP*ones(Ne,1);
    else
        TmIspk=TmIspk-L_t;
        TmEspk=TmEspk-L_t;
        TmPspk=TmPspk-L_t;
    end
for j=2:L_t
    %I-cells (Interneuron model); NO direct sin iput; E-input from granule cells only
    gEE=s_EE*sum(e_syn(:,1)); %allE-to-allI (FF)
   inRefracI = (j-TmIspk >= nrnSpace);
   v_I = v_I + dt/tau_m*(-v_I -gEE*(v_I-Esyn)+sig_a*eta_i).*inRefracI;
   
   eta_i = eta_i + dt*(-eta_i/tau_eta+sqtn*randn(Nei,1));
   %Two-step update:
   i_syn(:,1:end-1)=i_syn(:,2:end); %shift all cols to left
   i_syn(:,end)=i_syn(:,end-1)+dt/t_d*(-i_syn(:,end-1)+i_asyn); %usual Euler-Step
   i_asyn=i_asyn+dt/t_r*(-i_asyn);
   
   spk_Inr=(v_I>1);
   nuI(spk_Inr)=nuI(spk_Inr)+1;%store firing rates
   i_asyn(spk_Inr)=i_asyn(spk_Inr)+ampS; %jump upon spiking
   TmIspk(spk_Inr)=j;
   v_I(spk_Inr)=0; %reset
   
   %presy E-cells
   inRefracE = (j-TmEspk >= nrnSpace);
   v_Ei = v_Ei + dt/tau_m*(-v_Ei+I_aff(j)+sig_a*eta_e).*inRefracE;
   
   eta_e = eta_e + dt*(-eta_e/tau_eta+sqtn*randn(Nei,1));
   %Two-step update:
   e_syn(:,1:end-1)=e_syn(:,2:end); %shift all cols to left
   e_syn(:,end)=e_syn(:,end-1)+dt/t_d*(-e_syn(:,end-1)+e_asyn); %usual Euler-Step
   e_asyn=e_asyn+dt/t_r*(-e_asyn);
   
   spk_Inr=(v_Ei>1); %use same var spk_Inr
   nuE(spk_Inr)=nuE(spk_Inr)+1;%store firing rates
   e_asyn(spk_Inr)=e_asyn(spk_Inr)+ampS; %jump upon spiking
   TmEspk(spk_Inr)=j;
   v_Ei(spk_Inr)=0; %reset
   
   %target superfic Pyram
   Gee=qPreF.*(J_ee*e_syn(:,1)); %only use 1st col b/c of t_dely
   Gie=qPreF.*(J_ie*i_syn(:,1)); %only use 1st col b/c of t_dely
   inRefracP = (j-TmPspk >= nrnSpaceP);
   v_Pyr = v_Pyr + dt/tau_m*(-v_Pyr+I_aff(j)+sig_BP*eta_Pyr-Gee.*(v_Pyr-Esyn)-Gie.*(v_Pyr-Isyn)).*inRefracP;
   
   eta_Pyr = eta_Pyr + dt*(-eta_Pyr/tau_eta+sqtn*randn(Ne,1));
   
   spk_Inr=(v_Pyr>ThresE); %use same var spk_Inr
   f_rate(spk_Inr)=f_rate(spk_Inr)+1;%store firing rates
   nu_t(spk_Inr,ceil(j/spT))=nu_t(spk_Inr,ceil(j/spT))+1;
   TmPspk(spk_Inr)=j;
   v_Pyr(spk_Inr)=0; %reset

end
end
nuI=nuI./(N_rlz*t_end);
nuE=nuE./(N_rlz*t_end);
f_rate=f_rate./(N_rlz*t_end);
nu_t=nu_t./(dts*N_rlz);


end