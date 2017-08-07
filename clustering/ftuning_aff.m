function [presy_Frq_tuned,J_di,W_di]=ftuning_aff(qPreF,ThresE,Ne,Nei)
%determine the strength of I_aff and connect matrices (_ee,_ie) based on (qPreF,ThresE) relationship
%Based on pure distance to line slope: m_max->corr=1, -m_max->corr=-1, missing many (q,th) in upper(Quad1/2) and lower(Quad3/4)
%Have clustering (target cells with more presyn inputs or way less), cells
%closer to sing_pt have more, cells further away have less
%Assuming fr_min < frq < fr_max!

%Main Outputs
presy_Frq_tuned=zeros(Nei,1); %freq tuning of presyn pop (EGp)
J_di=sparse(Ne,Nei); %connec matrix for both E/I presyn (disynaptic)

%range of frequencies (fr_max is low)
fr_min=0;
fr_max=125; %in Hz
d_frq=2.5; %mesh incr of frq discr (Hz)
frq = (fr_min+d_frq/2: d_frq: fr_max-d_frq/2)'; %mesh for frq vals in connMat J_di
len_frqv=length(frq);
num_inOneF=floor(Nei/len_frqv); %number presyNrns tuned to frq(j)
presy_Frq_tuned(1:num_inOneF*len_frqv)= reshape(repmat(frq',num_inOneF,1),num_inOneF*len_frqv,1);
if(num_inOneF*len_frqv<Nei) %len_frqv doesn't divide Nei, so some left over at the end
    presy_Frq_tuned(num_inOneF*len_frqv+1:end)=-1; %set ones at end (-1) to respond to all freq
end

sing_pt=.5*[max(qPreF)+min(qPreF); max(ThresE)+min(ThresE)]; %determin sing pt (intersec) for mapping to corr value
unt_slp=(max(ThresE)-min(ThresE))/(max(qPreF)-min(qPreF)); %slope that gets mapped to corr=1
%the effective corr want to mimic, based on data (low=1, high=-1)
eff_Corr = -2*(frq-fr_min)./(fr_max-fr_min) + 1; %len_frqv x 1
eff_Slp = eff_Corr*unt_slp; %slope of line(s) (q,th) driven the hardest, len_frqv x 1

%-- Calc distances to line of neurons, for each frq --
dist_line=zeros(Ne,len_frqv); %distances to line(q,th) values that driven hardest
for j=1:len_frqv
    vc_eff = [1;eff_Slp(j)]; %vector for line (q,th) driven hardest (through origin)
    shf_or = [0; sing_pt(2)-sing_pt(1)*eff_Slp(j)]; %shift so everything through origin
    %put in matrices to calc distance to line of (q,th)
    mat_qt = [qPreF'; ThresE'] - repmat(shf_or,1,Ne); %2 x Ne
    scls_orthProj = mat_qt'*vc_eff./(vc_eff'*vc_eff); %Ne x 1; scalars for othog proj
    mat_normVc = mat_qt - repmat(vc_eff,1,Ne).*repmat(scls_orthProj',2,1); %2 x Ne
    dist_line(:,j) = sqrt(sum(mat_normVc.^2))';
end

%Generate connectivity matrices corresponding to clustering (based on (qPreF,ThresE) )
%J_di indicates it is disynaptic, so J_ee and J_ie are the SAME matrices (W_ too)
%select (qPreF,ThresE) values corresponding to effective 'correlation' (dist to line)
wdth_vr=0.1; %cells within this band driven/connect. max
for j=1:len_frqv
    dist_vr=dist_line(:,j);
    dist_vr(dist_line(:,j)<wdth_vr)=wdth_vr; %make set of cells driven same (at max)
    %aspr_effS=exp(-(dist_vr-wdth_vr)./.2); %Aspirational Effective Strength; made to be between 0 and 1
    
    prob_conn=exp(-(dist_vr-wdth_vr)./.01); %Ne x 1, prob of connection
    id_PresyNrns = (presy_Frq_tuned==frq(j)); %logical for presy neurons corresp. to this frq
    rndConn_mat=rand(Ne,sum(id_PresyNrns));
    J_di(:,id_PresyNrns) = 1*( rndConn_mat < repmat(prob_conn,1,sum(id_PresyNrns)) );  
end
id_PresyNrns = (presy_Frq_tuned==-1); %left-over neurons (no frq pref)
if(sum(id_PresyNrns)>0) %Some left over at the end, no frq tuning
    %find ones that don't receive any input and let the ones at the end
    %sparsely connect to them
    id_noInp = (sum(J_di,2)==0);
    if( sum(id_noInp)>0 ) %at least 1 SupPyr not receiving input
        rndConn_mat=rand(sum(id_noInp),sum(id_PresyNrns));
        nmGuaran_conn=min(sum(id_PresyNrns),round(.05*Nei)); %guaranteeting at least this many conn
        rndConn_mat(:,end-nmGuaran_conn+1:end)=0; %guaranteeting at least nmGuaran_conn conn
        J_di(id_noInp,id_PresyNrns) = 1*( rndConn_mat < 0.05 ); %5% conn, sparse
    else %all SupPyr receive input; make ones at end Erdos-Reyni
        rndConn_mat=rand(Ne,sum(id_PresyNrns));%sum(id_PresyNrns)=Nei-num_inOneF*len_frqv
        J_di(:,id_PresyNrns) = 1*( rndConn_mat < 0.05 ); %5% conn, sparse
    end
end

%augment J_di for mex program
l_wdi=max(sum(full(J_di)'));
W_di=zeros(Ne,l_wdi);
for j=1:Ne
   indx=find(J_di(j,:));
   W_di(j,1:length(indx))=indx;
end
%in C: 0,1,..,N-1 and -1 is where stop in C for-loop
W_di=W_di-1;

