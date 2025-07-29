function [phit_ics]=make_phittable(data,options)

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

% grid of selection coefficients
if isfield(data,'logs') logs=data.logs; else logs=(-8:1/16:-1)'; end;
s=10.^logs;

% for log10(s)<-6 we use the integrals and SFS of s=10^-6 and reduce
% log10(c) by (-6)-log10(s)
ls2i=@(ls) min(81,max(1,16*(ls+6)+1));


dens_file='dens.mat'; %Replace fiel if you're changing the MAF grid
load(dens_file,'dens');
norm_file=options.splinefile;%'splines.mat'; %Replace file if you're changing the normalizing integrals
load(norm_file,'sp');
tbl_file=options.tablefile;%'splines.mat'; %Replace file if you're changing the normalizing integrals
load(tbl_file,'ldunder_of_xb');

%data.zmat=repmat(abs(data.z),1,length(logs));
data.SNPs.z=abs(data.SNPs.z);

minz=sqrt(2) *erfcinv(5e-8);
phit_ics=zeros(length(data.SNPs.x),length(cvec),length(logs));


for sidx=1:length(logs(logs<=-1))
    %disp(sidx);
    phit_ics(:,:,sidx)=ppval(sp(1+max(0,sidx-find(logs==-6))),cvec+log10(data.SNPs.mj)+min(0,logs(sidx)+6));
end


%for i=1:length(data.SNPs.x)
%    if (mod(i,100)==0) disp(i); end
%Get the normalization table phit for hit i
%nrm=cell2mat(arrayfun(@(c) arrayfun(@(ls) ppval(sp(ls2i(ls)),c+min(0,ls+6)),logs),cvec+log10(data.SNPs.mj(i)), 'uniformoutput', false))';
%nrm(~isfinite(nrm))=0;
%phit_ics(i,:,:)=nrm; %or should it be nrm'?
%end



end
