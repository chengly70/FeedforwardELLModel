function [new_vc]=crtW_corr(x,y,rho_Target)
%[new_vc]=crtW_corr(x,y,rho_Target)
%given vectors: x,y of equal length and specified correl -1<rho_Target<1,
%returns vector new_vc (with mean&std y) so corr(x,new_vc)=rho_Target

if(rho_Target<-1 || rho_Target>1 || length(x)~=length(y) || length(x)==1)
    disp('x y need to be same length; rho_Target must be in range');
    return
end

ang_tg=acos(rho_Target); %(0,pi)

%center
y0=y-mean(y);
x0=x-mean(x);

%project onto x
vproj_x=x0'*y0/(norm(x0)^2)*x0;
ort_comp=y0-vproj_x;

%make unit length
x1=x0./(norm(x0));
ort1=ort_comp./(norm(ort_comp));

%new_vc=x1/tan(ang_tg)+ort1 + mean(y); %adding on mean doesn't change corr
new_vc=x1*cos(ang_tg)+ort1*sin(ang_tg) + mean(y);
std_c1=std(new_vc);

%scale new_vc to have same std as y
scal=std(y)/std_c1; %linearly related if mult. both (x1,ort1) equally
new_vc=x1*scal*cos(ang_tg)+ort1*scal*sin(ang_tg) + mean(y);