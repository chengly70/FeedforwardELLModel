%using a beta distribution for the frequency tuning of cells in EGp

load Parms_clust_int %assuming same presy_Frq_tunes for all Parms_
%same result for Parms_wedge_int.mat; just using the nubmers Nf, etc

fr_min=0;
fr_max=125; %in Hz

xrv=(0.001: .998/(fr_max-fr_min-1): 0.999)';
%xrv=(5/(fr_max-fr_min) : 100/(fr_max-fr_min)/50 : 105/(fr_max-fr_min))';

fr_v=xrv*(fr_max-fr_min)+fr_min;

os_ind=5;
cc=jet(length(fr_v));
sc=20;

figure
hold on
for j=(os_ind+1):(length(fr_v)-3)
    yrv=betapdf(xrv,sc,sc*(1-xrv(j))/xrv(j));
    tunCrv(:,j-os_ind)=yrv./(fr_max-fr_min);
    if(mod(j,5)==1)
        plot(fr_v,tunCrv(:,j-os_ind),'color',cc(j-os_ind,:),'LineWidth',4)
    end
end

% -- determ strength of aff input based on freq tuning --
fr_v=xrv*(fr_max-fr_min)+fr_min; %can change xrv, range of beta [0,1]
frq=5; %small chirp, freq in Hz
frq_b=(frq-fr_min)/(fr_max-fr_min);
yrv=betapdf(xrv,sc,sc*(1-frq_b)/frq_b);
tungCrv_low=yrv./(fr_max-fr_min);
frq=117.778; %big chirp
frq_b=(frq-fr_min)/(fr_max-fr_min);
yrv=betapdf(xrv,sc,sc*(1-frq_b)/frq_b);
tungCrv_hi=yrv./(fr_max-fr_min);

xrvF=(presy_Frq_tuned-fr_min)./(fr_max-fr_min);
cc1=jet(5);%jet(length(unique(presy_Frq_tuned)));
tunCrvFrq=zeros(length(xrvF),5);
frq_v=(2.5/2:2.5:125-2.5/2)';
%other Figure; used in paper
figure
hold on
kind=1;
for k=2:12:50
    frq_S=(frq_v(k)-fr_min)./(fr_max-fr_min);
    
    yrv=betapdf(xrvF,sc,sc*(1-frq_S)/frq_S);
    tunCrvFrq(:,kind)=yrv./(fr_max-fr_min);
    tunCrvFrq(:,kind)=1.5*tunCrvFrq(:,kind)./max(tunCrvFrq(:,kind));
    
    plot(tunCrvFrq(:,kind),'color',cc1(kind,:),'LineWidth',4)
    kind=kind+1;
end
set(gca,'FontSize',18)
set(gca,'XLim',[0.9 Nei+.1])
xlabel('Presynaptic Cell Index l')
ylabel('C(l,\phi )')



