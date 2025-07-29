function [data,idin]=preprocess(data,minx,maxld,minsnps,ex)


% Minimal MAF bin
if ~exist('minx','var')
    minxb=15;
else
    xrng=[10.^(-3:0.02:-2.9),10.^(-2.9:0.1:-1),0.1:0.025:0.5];xrng=unique(xrng);xrng(end)=xrng(end)+1e-10;
    minxb=find(xrng>=minx,1,'first');
end

% Minimal maximum LDscore
if ~exist('maxld','var')
    maxld=300;
end

% Traits to exclude
if ~exist('ex','var')
    ex=[];
end

% Minimal number of hits per trait
if ~exist('minsnps','var')
    minsnps=100;
end

% Remove hits that with low MAF, high LDscore, low z-score, or that belong
% to excluded traits
minz=sqrt(2) *erfcinv(5e-8);
idin=find((data.SNPs.ld<maxld).*(~ismember(data.SNPs.id,ex)).*(abs(data.SNPs.z)>=minz));
n=arrayfun(@(i) sum(data.SNPs.id(idin)==i),1:data.nt)'; %We need to check which traits don't have enough hits after culling hits

%Also, Remove traits with fewer than minsnps hits or where more than 10% of hits
%are on the same block
m=[];f=[];
for i=1:data.nt
[m(i),f(i)]=mode(data.SNPs.blk(data.SNPs.id==i)); f(i)=f(i)/sum(data.SNPs.id==i);
% f is the fraction of hits at the block with the largest number of hits
end
trts=find((f>0.1)'+(n<minsnps));

idin=find((data.SNPs.ld<maxld).*(~ismember(data.SNPs.id,[ex,trts])).*(abs(data.SNPs.z)>=minz));
data.SNPs=data.SNPs(idin,:); 
data.traits.n=arrayfun(@(i) sum(data.SNPs.id==i),1:data.nt)';
newid=cumsum(data.traits.n~=0).*(data.traits.n~=0);
data.SNPs.id=newid(data.SNPs.id);
data.nt=data.nt-sum(data.traits.n==0);
data.traits=data.traits(data.traits.n~=0,:);


% idin=find((data.SNPs.ld>=maxld).*ismember(data.SNPs.id,ex).*(abs(data.SNPs.z)<minz));
% 
% idin=(1:length(data.SNPs.x))';
% idout=find((data.SNPs.ld>=maxld)+ismember(data.SNPs.id,[ex,trts])+(abs(data.SNPs.z)<minz));
% idin(idout)=[];
% data.SNPs(idout,:)=[]; 
% data.traits.
% newid=cumsum(data.traits.n~=0).*(data.traits.n~=0);
% data.SNPs.id=newid(data.SNPs.id);
% data.nt=data.nt-sum(data.traits.n==0);
% data.traits=data.traits(trts,:);
% 
% 
% data.traits.n=arrayfun(@(i) sum(data.SNPs.id==i),1:data.nt)';
% 
% 
% 
% 
% 
% data.traits(trts,:)=[];
% idout=find();
% idin(idout)=[];
% data.SNPs(idout,:)=[];
% data.SNPs.id=data.SNPs.id-sum(trts <data.SNPs.id,2);
% data.nt=length(unique(data.SNPs.id));
% data.traits.n=arrayfun(@(i) sum(data.SNPs.id==i),1:data.nt)';


disp(['size of dataset is ',num2str(length(data.SNPs.x))]);

end

