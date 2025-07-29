function pres=get_pres(data,inf,options)

if isfield(data,'logs') logs=data.logs; else logs=(-8:1/16:0)'; end;
ls2i=@(ls) min(81,max(1,16*(ls+6)+1));

s=10.^logs;

minz=sqrt(2) *erfcinv(5e-8);


dens_file='dens.mat'; %Replace fiel if you're changing the MAF grid
load(dens_file,'dens');
dens=dens(ls2i(logs),:);%[repmat(dens(1,:),sum(logs<-6),1);dens;repmat(dens(end,:),sum(logs>-1),1)]; %size of matrix- (size of logs)*(no. of MAF bins)

c=10.^inf.c;
c=c(data.SNPs.id);
c=c(:); %make column vector - size (no. of hits)

p=inf.p;

m_rel=data.SNPs.m_rel.*data.SNPs.mj;

x=data.SNPs.x;
z=data.SNPs.z;

h=m_rel.*(2*x.*(1-x));
ch=c.*h; %column vector - size (no. of hits)

sig=sqrt(1+ch.*s'); %matrix - size (no. of hits)x(size of logs)

pxs=dens(:,data.SNPs.xb)'; %p(x,s) matrix - size (no. of hits)x(size of logs)

% The math behind the next equations:
% p(x,z)=\ins_s p(z|s,x)p(x|s)p(s)
% p(x,hit)=\int_(z>minz) p(x,z)= \ins_s \int_(z>minz) p(z|s,x)p(x|s)p(s)
% =\ins_s p(hit|s,x)p(x|s)p(s)
% p(z|x,hit)=p(x,z)/p(x,hit)
% CDF(z|x,hit)=\int_(z>z'>minz) p(x,z')/p(x,hit)
% \ins_s \int_(z>z'>minz) p(z'|s,x)p(x|s)p(s) / \ins_s \int_(z>minz) p(z|s,x)p(x|s)p(s)
% Since p(z'|s,x) is normal then the \int_(z>z'>minz) p(z'|s,x) =
% normcdf(z|s,x)-normcdf(minz|s,x)

zcdf=((normcdf(repmat(z,1,length(s)),0,sig)-normcdf(minz,0,sig)).*pxs)*p;
minzcdf=((1-normcdf(minz,0,sig)).*pxs)*p;

pres=1-zcdf./minzcdf;%(zcdf-minzcdf)./(1-minzcdf); 
% For TSD this gives us a (no of snps)x(no of traits sized array). 
% Where pres is evaluated with the f(s) (here p) for each possible trait.

if length(inf.p(1,:))>1 % if TSD
pres=pres(data.SNPs.id==(1:data.nt)); % we subset the matrix to include for each SNP only the right p.
%pres=pres(:); %make column vector
end

end