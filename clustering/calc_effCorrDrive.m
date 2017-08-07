function [frq_v,eff_corr]=calc_effCorrDrive(qPreF,ThresE,presy_Frq_tuned,J_di,frq)
%function calculates the effective (q,th) correl that is activated with each freq
%Corr is calc by weighing # presy. inputs to j (/Ne) as # copies of (q,th)
%similar to calc_effCorr.m but each input is WEIGHTED by a certain amount
%INPUT: presy_Frq_tuned is Nei x 1, J_di is Ne x Nei, with 0's and 1's (not MEX format)
%OUTPUT: (frq_v,eff_corr)

Nei=length(presy_Frq_tuned); %should be # cols of J_di
Ne=length(qPreF); %should be same as length of ThresE

%--- !MAKE sure this weighting scheme of aff_strgn is the same as scpt_Frq[].m ---
% see these parameters and lines in for-loop
d_frq=frq(2)-frq(1);
fr_min=frq(1)-d_frq/2;
fr_max=frq(end)+d_frq/2; %in Hz
frqGrid_presy=presy_Frq_tuned(presy_Frq_tuned>=0); %grid of freq values (gen ConnMatr)
x_rv=(frqGrid_presy-fr_min)./(fr_max-fr_min); %transform to beta (0,1)
sc=20; %beta distrib parameter
scl=115/frq(end); 
affStr_max=1.5;
%don't consider any Nei that have aff_strgn below min_affThres
min_affThres=0.01; %MUST be a whole number!!
mult_neiC=1/min_affThres;

%Outputs
frq_v=unique(presy_Frq_tuned); %sorted freq (-1,1.25,..,123.75)
len_frq=length(frq_v);
eff_corr=zeros(len_frq,1);

for j=1:len_frq
    %set to empty to start for each frq
    qEff=[];  %Effect. # of q's activated
    thEff=[]; %Effect. # of theta's activated
    
    frq_b=(frq(j)*scl-fr_min)/(fr_max-fr_min); 
    y_rv=betapdf(x_rv,sc,sc*(1-frq_b)/frq_b); %eval betapdf at frqGrid(transformed)=x_rv
    aff_strgn=zeros(Nei,1); %strength of afferent (sin) onto EGp (E&I cells)
    aff_strgn(presy_Frq_tuned>=0)=y_rv./(max(y_rv))*affStr_max; %relative strength, max=affStr_max
    aff_strgn(presy_Frq_tuned==-1)=0.2*affStr_max; %untuned ones set to xx*affStr_max
    aff_strgn(aff_strgn<min_affThres)=0; %cut-off below min_affThres
    numPresAct=round(mult_neiC*aff_strgn); %number of presyn activated
    
    f_ind=(numPresAct>0); %index of all presy with frq tuning frq-v(j)
    numPresAct=numPresAct(f_ind); %trim it down to a sum(f_ind) x 1 vector
    subMat_act=J_di(:,f_ind); %sub-matrix activated
    for k=1:sum(f_ind) %total # Nei cells activated
        qEff=[qEff; repmat(qPreF(subMat_act(:,k)==1),numPresAct(k),1)];
        thEff=[thEff; repmat(ThresE(subMat_act(:,k)==1),numPresAct(k),1)];
    end
    
    Mat_corr=corrcoef(qEff,thEff);
    eff_corr(j)=Mat_corr(1,2); %any of the 2 off-diagonals
end