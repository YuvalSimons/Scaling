function [p_ics,pout_i]=make_ptable(data,options)

arguments
data struct
options.splinefile char = 'splines.mat'
options.tablefile char = 'tables.mat'
end

%p_ics is intially a 3-d table of probability densities for z_i and x_i with the first index being snp i, the second being C
%parameter, and the third the selection coefficient

%phit_ics is intially a 3-d table of probabilities of being a hit given mj of SNP i with the first index being snp i, the second being C
%parameter, and the third the selection coefficient

%probability of seing snp i with z_i and x_i conditional on it being a hit
%is then (sum over s p_ics*f_S)/(sum over s phit_ics*f_S)

%grid of C values
cvec=3:0.01:9;

minz=sqrt(2) *erfcinv(5e-8);

% grid of selection coefficients
if isfield(data,'logs') logs=data.logs; else logs=(-8:1/16:-1)'; end;
s=10.^logs;
ls2i=@(ls) min(81,max(1,16*(ls+6)+1));

% for log10(s)<-6 we use the integrals and SFS of s=10^-6 and reduce
% log10(c) by (-6)-log10(s)
%ls2i=@(ls) max(1,16*(ls+6)+1);


dens_file='dens.mat'; %Replace file if you're changing the MAF grid
load(dens_file,'dens');
norm_file=options.splinefile;%'splines.mat'; %Replace file if you're changing the normalizing integrals
load(norm_file,'sp');
tbl_file=options.tablefile;%'tables.mat'; %Replace file if you're changing the normalizing integrals
load(tbl_file,'ldunder_of_xb');

%data.zmat=repmat(abs(data.z),1,length(logs));
data.SNPs.z=abs(data.SNPs.z);

minz=sqrt(2) *erfcinv(5e-8);
p_ics=zeros(length(data.SNPs.x),length(cvec),length(logs));
pout_i=zeros(length(data.SNPs.x),1);
for i=1:length(data.SNPs.x)
    %if (mod(i,100)==0) disp(i); end

%Get the probability density table p for hit i
%if mod(i,100)==0 disp(i); end
sh=2*(data.SNPs.mj(i)*data.SNPs.m_rel(i).*data.SNPs.x(i).*(1-data.SNPs.x(i)))*(10.^logs);
sig=sqrt(1+10.^(cvec).*sh);
Pa=normpdf(data.SNPs.z(i),0,sig);
newdens=zeros(size(logs));
newdens=dens(ls2i(logs),data.SNPs.xb(i));
%idmid=find((logs>=-6).*(logs<=-1));
%newdens(logs<-6)=dens(1,data.SNPs.xb(i));
newdens=newdens.*ldunder_of_xb(data.SNPs.xb(i));
pout_i(i)=newdens(1)*normpdf(data.SNPs.z(i),0,sqrt(1+data.zout^2))/erfc(erfcinv(5e-8)*sqrt(1/(1+data.zout^2)));
p_ics(i,:,:)=bsxfun(@times,Pa,newdens)';

end
pout_i=pout_i/(286847/120000000); %normalize by p(x>0.01) for s=1e-6. Replace with code if this process works


end
