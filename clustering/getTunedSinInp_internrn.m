function [I_aff,I_ie,I_ee,tV,nuI,nuE]=getTunedSinInp_internrn(frq,dc_val,Amp,Nei,aff_strgn)
%Input: cutoff_frq (Hz); dc_val (limit of direct input); amp_rng (amplitude
%range); Nei is
%size of each E & I Presyn Pop, aff_strgn is an Nei x 1 vector for the freq-tuning of inputs
% [aff_strgn is created in scpt_vryFrq.m,etc.
%-- Function to determine effective current input to E cell ---
%I_aff band-limited Random Amplitude Modulation (RAM)
%[units of norm. volt [0,1], assuming feedforward architechure;
%Nei is the number of each presyn E and presyn I neurons. 
% addressing Reviewer 1's interneuron concern; no longer receives direct
% sinusoidal input, but we have new parameter s_EE
%!!! ALL Thresholds are 1 !!!

if(size(aff_strgn,1)~=Nei || size(aff_strgn,2)~=1)
    disp('aff_strgn must be an Nei x 1 vector');
    return;
end

%Choice of t_end insures a whole number of cycles; only a problem when freq has decimal
t_end=1; %in sec; 
dt=.0001; %in sec
tV=(dt:dt:t_end)';

%biophys params
%Nei=200; %numb of Inhib neurons
nuI=zeros(Nei,1); %firing rate of I-cells
nuE=zeros(Nei,1); %firing rate of E-cells
tau_m=0.01; %mem. time const. in sec
t_ref=0.0005; %ref period in sec
t_r=0.002; %Both synapses; rise time
t_d=0.01;  %Both synapses; decay time
ampS=2; %jump in asyn upon spiking
sig_a=1;
tau_eta=0.005;
Esyn=6.5;
sqtn=1/sqrt(dt*tau_eta);

v_I=rand(Nei,1);%zeros(Nei,1);
i_syn=zeros(Nei,1);
i_asyn=zeros(Nei,1);
eta_i=zeros(Nei,1);
I_ie=zeros(Nei,length(tV)); %output of i pop, without scalar (g_ie)
v_Ei=rand(Nei,1);
e_syn=zeros(Nei,1);
e_asyn=zeros(Nei,1);
eta_e=zeros(Nei,1);
I_ee=zeros(Nei,length(tV)); %output of E pop, without scalar (g_ee)

%1 new param; interneurons
s_EE=2.7/Nei;

% Afferent input I_aff to E (& I) cells; frozen RAM
I_aff=dc_val+Amp*sin(2*pi*tV*frq);
I_aff(I_aff<0)=0; %threshold so f-rate incr

nrnSpace=floor(t_ref/dt);
TmIspk=-nrnSpace*ones(Nei,1);
TmEspk=-nrnSpace*ones(Nei,1);

for j=2:length(tV)
    %I-cells (Interneuron model); NO direct sin iput; E-input from granule cells only
    gEE=s_EE*sum(I_ee(:,j-1)); %allE-to-allI (FF)
   inRefracI = (j-TmIspk >= nrnSpace);
   v_I = v_I + dt/tau_m*(-v_I -gEE*(v_I-Esyn) +sig_a*eta_i).*inRefracI;
   
   eta_i = eta_i + dt*(-eta_i/tau_eta+sqtn*randn(Nei,1));
   i_syn=i_syn+dt/t_d*(-i_syn+i_asyn);
   i_asyn=i_asyn+dt/t_r*(-i_asyn);
   
   spk_Inr=(v_I>1);
   nuI(spk_Inr)=nuI(spk_Inr)+1;%store firing rates
   i_asyn(spk_Inr)=i_asyn(spk_Inr)+ampS; %jump upon spiking
   TmIspk(spk_Inr)=j;
   v_I(spk_Inr)=0; %reset
   
   I_ie(:,j)=i_syn;
   
   %E-cells
   inRefracE = (j-TmEspk >= nrnSpace);
   v_Ei = v_Ei + dt/tau_m*(-v_Ei+I_aff(j)*aff_strgn+sig_a*eta_e).*inRefracE;
   
   eta_e = eta_e + dt*(-eta_e/tau_eta+sqtn*randn(Nei,1));
   e_syn=e_syn+dt/t_d*(-e_syn+e_asyn);
   e_asyn=e_asyn+dt/t_r*(-e_asyn);
   
   spk_Inr=(v_Ei>1); %use same var spk_Inr
   nuE(spk_Inr)=nuE(spk_Inr)+1;%store firing rates
   e_asyn(spk_Inr)=e_asyn(spk_Inr)+ampS; %jump upon spiking
   TmEspk(spk_Inr)=j;
   v_Ei(spk_Inr)=0; %reset
   
   I_ee(:,j)=e_syn;
end
nuI=nuI./t_end;
nuE=nuE./t_end;