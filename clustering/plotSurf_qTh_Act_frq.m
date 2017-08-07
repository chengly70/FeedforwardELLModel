%script to plot effective strength or # connection as freq varies

%load Parms_clust_int
load Parms_wedge_int

fr_min=0;
fr_max=125; %in Hz
d_frq=2.5; %mesh incr of frq discr (Hz)
frq = (fr_min+d_frq/2: d_frq: fr_max-d_frq/2)'; %mesh for frq vals in connMat J_di
len_frqv=length(frq);
affStr_max=1.5;

frqGrid_presy=presy_Frq_tuned(presy_Frq_tuned>=0); %grid of freq values (gen ConnMatr)
x_rv=(frqGrid_presy-fr_min)./(fr_max-fr_min); %transform to beta (0,1)
sc=20; %param for beta distr
scl=115/frq(end);
ind_ne=zeros(Ne,1);
gray_max=0.9; %max for gray-scale color, MUST be less than 1 (1 is white)

one_vec=ones(Ne,1);   

%for surf to make visible
oneS=ones(10,10);
qMat=repmat((qPreF(1): 1/9: qPreF(end))',1,10);
min_th=min(ThresE);
max_th=max(ThresE);
thMat=repmat((min_th: (max_th-min_th)/9: max_th),10,1);

figure
hold on
grid on 
for kind=1:3
    if(kind==1)
        k=2; 
    elseif(kind==2)
        k=25;
    else
        k=50;
    end
    %generate different strengths based on beta-Frq-tuning (tuningCrvs_beta.m)
    frq_b=(frq(k)*scl-fr_min)/(fr_max-fr_min); %assuming same as on lines4,5
    y_rv=betapdf(x_rv,sc,sc*(1-frq_b)/frq_b); %eval betapdf at frqGrid(transformed)=x_rv
    aff_strgn=zeros(Nei,1); %strength of afferent (sin) onto EGp (E&I cells)
    aff_strgn(presy_Frq_tuned>=0)=y_rv./(max(y_rv))*affStr_max; %relative strength, max=affStr_max
    aff_strgn(presy_Frq_tuned==-1)=0.2*affStr_max; %untuned ones set to xx*affStr_max
    
    eff_Inp_act=J_di*aff_strgn; %effective matrix activated
    eff_Inp_act=eff_Inp_act./max(eff_Inp_act); %normalize so in [0,1]
    eff_Inp_act=round(eff_Inp_act*10)./10; %truncate to 1 decimal
    uniq_effIA=unique(eff_Inp_act); %unique list of weighted values
    
    %cc=parula(length(uniq_effIA));
    cc=jet(length(uniq_effIA));
    for j=1:length(uniq_effIA)
        ind_ne = (eff_Inp_act == uniq_effIA(j)); %logical of (q,th) pairs
        %plot(qPreF(ind_ne),ThresE(ind_ne),'.','MarkerSize',18,'color', ones(1,3)*(gray_max*(1-uniq_effIA(j))) )
        plot3(frq(k)*one_vec(ind_ne),qPreF(ind_ne),ThresE(ind_ne),'.','MarkerSize',18,'color', cc(j,:) )
    end
    surf((frq(k)+1)*oneS,qMat,thMat)
    
end
set(gca,'FontSize',18)
xlabel('Frequency (Hz)')
ylabel('q')
zlabel('\theta')

