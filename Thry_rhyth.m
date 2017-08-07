%Test the theory

% !!these lines have to match what was used in sims!!
%altering intrins & net. heterog. relationship
Ne=1000;
perc_th=(0.05: 0.9/(Ne-1): 0.95)';
sig_nm=0.1;
qPreF=(0.5: 1/(Ne-1): 1.5)';
ThresE0=exp(norminv(perc_th,-sig_nm^2/2,sig_nm)); %log-normal with mean exp(0)=1
load storedPerm
Thres1=ThresE0(strdPerm); %crtW_corr function works better with random samples of same vector

corr_preTh=(-.9:.05:.9)';
len_c=length(corr_preTh);

rng_thry=zeros(len_c,1);
std_thry=zeros(len_c,1);

for k=1:len_c
    ThresE=crtW_corr(qPreF,Thres1,corr_preTh(k));
    rng_thry(k,1)=max(qPreF./ThresE)-min(qPreF./ThresE);
    std_thry(k,1)=std(qPreF./ThresE);
end

figure
plot(corr_preTh,std_thry,'k','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Correlation (q,Thres)')
ylabel('Theory for Std. Dev. of Rates')
