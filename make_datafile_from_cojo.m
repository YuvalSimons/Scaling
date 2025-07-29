function data=make_datafile_from_cojo(cojo_folder,trait_list,Mmed_list,block_file,datafile,LDscore_file,name_list,h2_list)

% cojo_folder: folder with unified cojo files, one file per trait
% trait_list:list: list of such traits to include
% Mmed_list: list of median reported study sizes 
% block_file: File with genomic blocks, used for bootstrapping and jackknifing. Berisa(2016) blocks file supplied with code.
% datafile: output file, if not provided none saved (optional)
% LDscore_file: File with LDscores. If none provided, all LDscores set to 1. (optional)
% name_list: common language names for those trait (optional)
% h2_list: estimates of heritability for those traits, to compare our estimates with (optional)


opts = delimitedTextImportOptions('VariableNames',{'trait'},'delimiter',' ');
traits=readtable(trait_list,opts,'ReadVariableNames' ,false);
%opts = delimitedTextImportOptions('VariableNames',{'M_med','ReportedM_med'},'delimiter',',');
Mmeds=csvread(Mmed_list);%readtable(Mmed_list,opts,'ReadVariableNames' ,false);
M_med=Mmeds(:,1);
if exist('name_list','var')
    if ~isempty(name_list)
    opts = delimitedTextImportOptions('VariableNames',{'name'},'delimiter',' ');
    names=readtable(name_list,opts,'ReadVariableNames' ,false); 
    end
end
if exist('h2_list','var')
    if ~isempty(h2_list)
    %opts = delimitedTextImportOptions('VariableNames',{'h2'},'delimiter',' ');
    h2=csvread(h2_list);%readtable(h2_list,opts,'ReadVariableNames' ,false); 
    end
end

SNPs=[];


fct=[];

% For each trait
for i=1:length(traits.trait)
    trt=traits.trait{i};
    disp([num2str(i),' ',trt]);
    % Read the cojo file
    gwas=readtable([cojo_folder,trt,'.jma.cojo'],'FileType','text');%readtable(['/Users/yuvalsim/PhD/Inference/Writeup/cojo2022/',trt,'.jma.cojo'],'FileType','text');
    gwas.Properties.VariableNames([1])={'chr'};
    gwas.Properties.VariableNames([5])={'x'};
    gwas.x=min(gwas.x,1-gwas.x); %Convert to MAF
    gwas.id=i*ones(size(gwas.x)); %Create trait ID
    gwas.m_rel=1./(M_med(i).*2*(gwas.se.^2).*gwas.x.*(1-gwas.x));
    gwas.mj=(gwas.se./gwas.bJ_se).^2;
    %trtgwas(trtgwas.m_rel<0.8,:)=[]; DO IN PREPROCESS
    %trtgwas(trtgwas.x<0.01,:)=[];DO IN PREPROCESS
    
    %minz=sqrt(2) *erfcinv(5e-8); DO IN PREPROCESS
    %trtgwas.z=trtgwas.bJ./trtgwas.bJ_se; DO IN PREPROCESS
    %trtgwas(abs(trtgwas.z)<minz,:)=[]; DO IN PREPROCESS
    
    SNPs=[SNPs;gwas];
end

SNPs.loc=SNPs.chr*1e9+SNPs.bp; %One unique number for the genomic location of each SNP. Makes life easy. RS numbers suck.
xrng=[10.^(-3:0.02:-2.9),10.^(-2.9:0.1:-1),0.1:0.025:0.5];xrng=unique(xrng);xrng(end)=xrng(end)+1e-10;
SNPs.xb=arrayfun(@(xi) sum(xi>=xrng),SNPs.x);


%Get block number for each SNP
if ~exist(block_file,'file')
    block_file='blocks/blocks.tsv';
end   
tmp=readmatrix(block_file,'filetype','text','TrimNonNumeric',true);
blocks=array2table(tmp,'variablenames',{'chr','down','up'});
SNPs.blk=zeros(size(SNPs.loc));
blknum=length(blocks.chr);
for i=1:blknum
    idx=find((blocks.chr(i)==SNPs.chr).*(blocks.down(i)<=SNPs.bp).*(blocks.up(i)>SNPs.bp));
    SNPs.blk(idx)=i;
end


if exist('LDscore_file','var')
LDscores=readtable(LDscore_file);%'/Users/yuvalsim/PhD/Inference/Writeup/Github_files/ldscores/ukbb_new_ldscores.txt'%'/Users/yuvalsim/Downloads/TWINSUK.beagle.anno.csq.shapeit.20131101.quad.corrected.score.ld.txt','FileType','text','delimiter',' ');
LDscores.loc=LDscores.chr*1e9+LDscores.bp;
[Loca,Locb] = ismember(SNPs.loc,LDscores.loc);
ld=LDscores.ldscore(Locb(Locb>0));
SNPs=SNPs(Loca,:);
SNPs.ld=ld;
else
    SNPs.ld=ones(size(SNPs.x));
end

data=[];
SNPs.a=abs(SNPs.bJ); SNPs.a_se=abs(SNPs.bJ_se); SNPs.z=SNPs.a./SNPs.a_se;  SNPs.v=2*(SNPs.bJ.^2).*SNPs.x.*(1-SNPs.x);  SNPs.dv=2*(SNPs.bJ_se.^2).*SNPs.x.*(1-SNPs.x);

data.SNPs=SNPs;

data.nt=max(data.SNPs.id);

traits.n=arrayfun(@(trt)  sum(data.SNPs.id==trt),1:data.nt)';
traits.M_med=M_med; traits.ReportedM_med=Mmeds(:,2);
if exist('names','var') traits.name=names.name; end
if exist('h2s','var') traits.h2=h2s.h2; end
data.traits=traits;

if exist('datafile','var') 
    if ~isempty(datafile)
        save(datafile,'data'); 
    end
end

end