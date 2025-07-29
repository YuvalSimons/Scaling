function make_normalizing_integrals_stage1(trait_file,LDscore_file,fileout,maxld)

if ~exist('maxld','var') maxld=300; end

tbl=readtable(trait_file,'FileType','text');

if length(tbl.Properties.VariableNames)==12 tbl(:,4)=[]; end
tbl.Properties.VariableNames={'SNP','minor_allele','x','low_con','n','AC','ytx','a','se','z','p'};
tbl.low_con=[]; tbl.AC=[]; tbl.ytx=[];
isnotXchr=cellfun(@(x) x(1)~='X', tbl.SNP);
tbl=tbl(isnotXchr,:);
splitSNP=split(tbl.SNP, ':');
tbl.chr = str2double(splitSNP(:, 1)); 
tbl.bp = str2double(splitSNP(:, 2));
tbl.loc=tbl.chr*1e9+tbl.bp;

tbl.m=1./(2*tbl.x.*(1-tbl.x).*(tbl.se.^2));

x=min(tbl.x,1-tbl.x);
tbl(~isfinite(tbl.m),:)=[];

m_med=median(tbl.m);

tbl.m_rel=tbl.m/m_med;


LDscores=readtable(LDscore_file);%'/Users/yuvalsim/PhD/Inference/Writeup/Github_files/ldscores/ukbb_new_ldscores.txt'%'/Users/yuvalsim/Downloads/TWINSUK.beagle.anno.csq.shapeit.20131101.quad.corrected.score.ld.txt','FileType','text','delimiter',' ');
LDscores.loc=LDscores.chr*1e9+LDscores.bp;
[Loca,Locb] = ismember(tbl.loc,LDscores.loc);
tbl.ld=ones(size(tbl.x));
tbl.ld(Loca)=LDscores.ldscore(Locb(Locb>0));

eps=erfcinv(5e-8);
aprop=@(mcsx) eps*sqrt(1./(1+mcsx));
%vprop=@(mcsx) mcsx.*(erfc(aprop(mcsx))+sqrt(4/pi).*(mcsx./(1+mcsx)).*aprop(mcsx).*exp(-aprop(mcsx).^2));
vprop_ap=@(ap) (erfc(ap)+sqrt(4/pi).*(1-(ap/eps).^2).*ap.*exp(-ap.^2));
vprop=@(mcsx) vprop_ap(aprop(mcsx));
nprop=@(mcsx) erfc(aprop(mcsx));

xrng=[10.^(-3:0.02:-2.9),10.^(-2.9:0.1:-1),0.1:0.025:0.5];xrng=unique(xrng);xrng(end)=xrng(end)+1e-10;
tbl.xb=zeros(size(tbl.x));
tbl.xb(tbl.x>=0.001)=arrayfun(@(xi) sum(xi>=xrng),tbl.x(tbl.x>=0.001));

m_rel_of_xb={}; ld_of_xb={}; ntbl=[];vtbl=[];
for i=1:40
    disp(i);
    m_rel_of_xb{i}=tbl.m_rel(tbl.xb==i);
    ld_of_xb{i}=tbl.ld(tbl.xb==i);
    ntbl(i,:)=arrayfun(@(log_mcsx) mean((ld_of_xb{i}<maxld).*nprop(m_rel_of_xb{i}*10^log_mcsx)),-10:0.001:10);
    disp(i);
    vtbl(i,:)=arrayfun(@(log_mcsx) mean((ld_of_xb{i}<maxld).*vprop(m_rel_of_xb{i}*10^log_mcsx)),-10:0.001:10); % We're not using this in the ccurrent inference, but I think it's nice to have predictions for the proportion of variance captured
    disp(i);
end

ldunder_of_xb=arrayfun(@(i) mean((ld_of_xb{i}<maxld)),1:40);

save(fileout,'m_rel_of_xb','ld_of_xb','ldunder_of_xb','m_med','ntbl','vtbl'); %this makes the file tables.mat
