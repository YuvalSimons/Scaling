function inf=do_inference(data,p_ics,phit_ics,pout_i)


%grid of C values
cvec=3:0.01:9;

% grid of selection coefficients
logs=data.logs;
mingrid=data.mingrid;
maxgrid=data.maxgrid; 
s=10.^logs;
ls2i=@(ls) max(1,16*(ls+6)+1);



mask=zeros([data.nt,length(data.SNPs.x)]);
for trt=1:data.nt
    mask(trt,data.SNPs.id==trt)=1;
end
f=@(v) reshape(v,[length(data.SNPs.x),length(cvec)]);

y2p=@(y) y_to_p(y,logs,mingrid,maxgrid);
%LL=@(rho) get_LL(rho,data,mask,p_ics,phit_ics);
%if ~isfield(data,'eps') data.eps=0; end
%if ~isfield(data,'k') data.k=4; end

%sh=2*(data.SNPs.m_rel.*data.SNPs.x.*(1-data.SNPs.x))*s';
%zmat=repmat(abs(data.SNPs.z),1,length(s));

%data.shtot=s.*data.htot;



%xrng=[10.^(-3:0.02:-2.9),10.^(-2.9:0.1:-1),0.1:0.025:0.5];xrng=unique(xrng);xrng(end)=xrng(end)+1e-10;
%xb=arrayfun(@(xi) sum(xi>=xrng),data.SNPs.x);
%data.xb=xb;
%data.dens=(1/120000000)*data.dens(:,xb)';
%data.dens=[repmat(data.dens(:,1),1,sum(data.logs<-6)),data.dens];


%dsp='iter';%'final';
%opts=optimset('maxiter',5000,'MaxFunEvals',10000,'TolFun',1e-4,'TolX',1e-4,'display','iter');

if ~exist('yin','var')
yin=zeros(data.k,1);
end


%Parameters for the simulated annealing
Tinit=100;%10000;
Tend=0.01;
D=1000./length(data.SNPs.x);% The more SNPs you the smaller the step you have to take to find the maximum since LL~(number of hits)*p(x,b|hit) %0.01;
Ts=Tend*(Tinit/Tend).^(0:(1/11):1);
betas=1./Ts;
LL0=get_LL(y2p(yin));%,data,mask,p_ics,phit_ics);
ymax=yin;
LLmax=LL0;

yin=repmat(yin,1,length(Ts));
%k=0.01*D/2;
sz=size(yin);

LLin=LL0*ones(size(Ts));

maxevals=data.funevals;%24000;

%figure;
%plot(T,LLin-LL0,'.k'); hold on; xlim([Tend,Tinit]); ylim([-3000,3000]); set(gca,'xscale','log'); 


%LLins=zeros(ceil(maxevals/length(Ts)),length(Ts));
%LLbelow=[];
LLlst=zeros(ceil(maxevals/length(Ts)),length(Ts));
acceptedlst=zeros(ceil(maxevals/length(Ts)),length(Ts));
Dlst=zeros(ceil(maxevals/length(Ts)),1);
Tlst=zeros(ceil(maxevals/length(Ts)),length(Ts));
switchlst=zeros(ceil(maxevals/length(Ts)),2);
evals=0;

i=1;
while evals<maxevals
    yout=yin+sqrt(D*Ts).*randn(sz);
    pout=y2p(yout);
    for j=1:length(Ts)
        rho=pout(:,j);
%         p_ic=f(p_ics*rho);
%         phit_ic=f(phit_ics*rho);
%         p=p_ic./phit_ic;
%         p(~isfinite(p))=0;
%         A=mask*log(p+realmin);
%         
%         penalty=0;
%         
%         L=((data.traits.n).^2)./(mask*phit_ic); % per trait, and C - L is number of hits /mean(phit) where mean is taken over trait hits
%         penalty=-data.traits.n.*L/data.targetsize;
%         
%         vperbp=rho'*data.shtot; 
%         V=vperbp*L.*(10.^cvec)./data.traits.ReportedM_med;
%         penalty=penalty-data.traits.n.*V/1;
%         
%         
%         ll=max(A+data.eps*penalty,[],2,'omitnan');
        LLout(j)=get_LL(rho);%sum(ll);
    end
    dLL=LLout-LLin;
    accepted=(exp(betas.*dLL)>rand(size(Ts)));
    idx=find(accepted);
    yin(:,idx)=yout(:,idx);%(1-accepted).*yin+accepted.*yout;
    LLin(idx)=LLout(idx);%(1-accepted).*LLin+accepted.*LLout;

    indxs=randsample(1:length(Ts),2); 
    k=indxs(1); l=indxs(2);
    if exp(-(LLin(k)-LLin(l)).*(betas(k)-betas(l)))>rand()
        Ts(indxs)=fliplr(Ts(indxs));
        betas(indxs)=fliplr(betas(indxs));
        %LLin(indxs)=fliplr(LLin(indxs));
        switchlst(i,:)=[k l];
    else
        switchlst(i,:)=[0 0];
    end

    if prod(1-accepted)==1
        D=D/1.25;
    elseif prod(accepted)==1
        D=D*1.25;
    end

    acceptedlst(i,:)=accepted; %acceptedlst=[acceptedlst;accepted];
    LLlst(i,:)=LLin;
    Tlst(i,:)=Ts;
    Dlst(i)=D;

    [m,id]=max(LLin);
    if m>LLmax
        LLmax=m;
        ymax=yin(:,id);
    end

    evals=evals+length(Ts);
    %disp([num2str(evals),' ',num2str(LLmax-LL0)]);
    if mod(i,100)==1        
        disp([num2str(evals),' ',num2str(LLin-LL0)]);
        disp([num2str(evals),' ',num2str(Ts,2)]);
        
        %mean_accepted=mean(acceptedlst((i-99):i,:),'all');
        disp(['D=',num2str(D)]);
    end

    i=i+1;
end

[inf.LL,inf.c,inf.L,inf.V,inf.Pout]=get_LL(y2p(ymax));%,data,mask,p_ics,phit_ics);
inf.LL=inf.LL-LL0;
inf.yout=ymax;
inf.ys={yin};
inf.ps={y2p(yin)};
inf.p=y2p(ymax);

inf.logs=logs;
inf.mingrid=mingrid; inf.maxgrid=maxgrid;
inf.mu=logs'*inf.p/sum(inf.p);
inf.sig=sqrt(((logs-inf.mu).^2)'*inf.p/sum(inf.p));

inf.LL0=LL0;
inf.LLlst={LLlst-LL0}; 
inf.acceptedlst={acceptedlst};
inf.Tlst={Tlst};
inf.switchlst={switchlst};

function [LL,c,L,V,Pout]=get_LL(rho)%,data,mask,p_ics,phit_ics,pout_i)%get_LL(data,cvec,pis,nrm,CDFminz,CDFa,nref,p)

cvec=3:0.01:9;

%rho(find((logs<-6)+(logs>-2)))=0;

%rho(~isfinite(rho))=0;

LL=0;
cmax=zeros(data.nt,1);

f=@(v) reshape(v,[length(data.SNPs.x),length(cvec)]);
p_ic=f(p_ics*rho);
phit_ic=f(phit_ics*rho);

p=p_ic./phit_ic;
p=data.pout*pout_i+(1-data.pout)*p;
p(~isfinite(p))=0;

A=mask*log(p+realmin);

penalty=0; L=0; V=0;

L=((data.traits.n).^2)./(mask*phit_ic); % per trait, and C - L is number of hits /mean(phit) where mean is taken over trait hits
penalty=-data.traits.n.*L/data.targetsize;

vperbp=rho'*data.shtot; 
V=vperbp*L.*(10.^cvec)./data.traits.ReportedM_med;
penalty=penalty-data.traits.n.*V/1;


[ll,kmax]=max(A+data.eps*penalty,[],2,'omitnan');
LL=sum(ll);

if nargout>1 
c=cvec(kmax);
L=arrayfun(@(trt) L(trt,kmax(trt)),1:data.nt);
V=arrayfun(@(trt) V(trt,kmax(trt)),1:data.nt);
P=zeros(size(data.SNPs.x));
P=arrayfun(@(j) p(j,kmax(data.SNPs.id(j))),1:length(data.SNPs.x))';
Pout={(data.pout*pout_i)./P};%+(1-data.pout)*p;
end


%minz=sqrt(2) *erfcinv(5e-8);

%p1=p1(sub2ind(size(p1),(1:length(data.x))',kmax(data.id)));
%rhon=rho;
%rhon(11:end)=0;
%pn=f(mat*rhon)./(rho'*nrm);
%pn=pn(sub2ind(size(pn),(1:length(data.x))',kmax(data.id)));






end


end

function [LL,c,L,V]=old_get_LL(rho,data,mask,p_ics,phit_ics,pout_i)%get_LL(data,cvec,pis,nrm,CDFminz,CDFa,nref,p)

cvec=3:0.01:9;

%rho(find((logs<-6)+(logs>-2)))=0;

%rho(~isfinite(rho))=0;

LL=0;
cmax=zeros(data.nt,1);

f=@(v) reshape(v,[length(data.SNPs.x),length(cvec)]);
p_ic=f(p_ics*rho);
phit_ic=f(phit_ics*rho);

p=p_ic./phit_ic;
p(~isfinite(p))=0;

A=mask*log(p+realmin);

penalty=0; L=0; V=0;

L=((data.traits.n).^2)./(mask*phit_ic); % per trait, and C - L is number of hits /mean(phit) where mean is taken over trait hits
penalty=-data.traits.n.*L/data.targetsize;

vperbp=rho'*data.shtot; 
V=vperbp*L.*(10.^cvec)./data.traits.ReportedM_med;
penalty=penalty-data.traits.n.*V/1;


[ll,kmax]=max(A+data.eps*penalty,[],2,'omitnan');
c=cvec(kmax);
sum(ll);

L=arrayfun(@(trt) L(trt,kmax(trt)),1:data.nt);
V=arrayfun(@(trt) V(trt,kmax(trt)),1:data.nt);

%minz=sqrt(2) *erfcinv(5e-8);

%p1=p1(sub2ind(size(p1),(1:length(data.x))',kmax(data.id)));
%rhon=rho;
%rhon(11:end)=0;
%pn=f(mat*rhon)./(rho'*nrm);
%pn=pn(sub2ind(size(pn),(1:length(data.x))',kmax(data.id)));


LL=sum(ll);



end

function p=y_to_p(y,logs,mingrid,maxgrid)

grid=[mingrid:(maxgrid-mingrid)/length(y(:,1)):maxgrid];%[-7.5:6/length(y):-1.5];%[-7.5:1.5:-1.5];%[-7.5,-6,-4.5,-3,-1.5];%[-8:-1];
idx=round(length(y(:,1))/2+1);
y=[y(1:idx-1,:);zeros(1,length(y(1,:)));y(idx:end,:)];

logrho=-spline(grid,y',logs')';%spline    %The transpose madness here is becuase Matlab is inconsistent in it's treating a vector and matrix y in spline
logrho=700*tanh(logrho/700);

%ids=find(logs);%find((logs>-6).*(logs<=-2));
p=zeros(size(logrho));
%p(ids)=exp(logrho(ids));

p=exp(logrho);

p(~isfinite(p))=0;
p=p./sum(p);
end

function p=experimental_y_to_p(y)

logs=(-8:1/16:-1)';

grid=[-6:4/length(y):-2];%[-7.5:6/length(y):-1.5];%[-7.5:1.5:-1.5];%[-7.5,-6,-4.5,-3,-1.5];%[-8:-1];
idx=round(length(y)/2+1);
y=-[y(1:idx-1);0;y(idx:end)];

y(2)=y(1)+abs(y(2)-y(1)); %second knot must be above first knot
y(end-1)=y(end)+abs(y(end-1)-y(end)); %second before last knot must be above last knot

%yflip=y; yflip(idx)=[]; %outputing y such that 
% the second knot is above the first knot
% and the second before last knot is above the last knot

logrho=spline(grid,y,logs);%spline  
logrho=700*tanh(logrho/700); % can't allow it to get so big that we overflow

ids=find((logs>-6).*(logs<=-2));
p=zeros(size(logrho));
p(ids)=exp(logrho(ids));
p(~isfinite(p))=0;
p=p/sum(p);
end



