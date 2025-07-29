function make_normalizing_integrals_stage2(simdatafile,tablefile,folderout,i)

load(simdatafile,'data');

load(tablefile,'m_rel_of_xb','ntbl','vtbl');

N0=3718274;%613285; % CHANGE HERE IF TENNESSEN
x=1:N0;x=x/(2*N0);
minx=0.01;maxp=5e-8;
xrng=[10.^(-3:0.02:-2.9),10.^(-2.9:0.1:-1),0.1:0.025:0.5];xrng=unique(xrng);xrng(end)=xrng(end)+1e-10;


    N=data.N{i};
    idx=find((N>0).*(x>=minx));
    n=N(idx);
    x=x(idx);
    v=2*x.*(1-x);
    runs=data.runs(i);
    
    s=10^data.s(i);
    
    for j=1:length(xrng)-1
        idb{j}=find((x>=xrng(j)).*(x<xrng(j+1)));
    end
    
    tic;
	mc_lst=0:0.001:10;
    f=[];
    fv=[];
    for log10mc=mc_lst
        disp(log10mc);
        mc=10^log10mc;
        %np=zeros(size(v));
        %vp=zeros(size(v));
        fb=zeros(1,length(xrng)-1);
        fvb=zeros(1,length(xrng)-1);
        for j=1:length(xrng)-1
            fb(j)=(1/runs)*n(idb{j})*spline(-10:0.001:10,ntbl(j,:),log10(mc*s*v(idb{j})))';
            fvb(j)=(1/runs)*n(idb{j})*(v(idb{j}).*spline(-10:0.001:10,vtbl(j,:),log10(mc*s*v(idb{j}))))';
        end
        
        f(end+1)=sum(fb);%[f,(1/runs)*n*np'];
        fv(end+1)=sum(fvb);%fv=[fv,(1/runs)*n*(v.*vp)'];
    end
    toc;
 
 sp=spline(mc_lst,f);
 spv=spline(mc_lst,fv);

 save([folderout,'/',num2str(i),'.mat'],'sp','spv','f','fv');  %this makes the file splines.mat from tables.mat

end

