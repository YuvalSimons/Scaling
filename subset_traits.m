function [data]=subset_traits(data,trts)
if length(trts)>0
    idin=find(ismember(data.SNPs.id,trts));
    data.SNPs=data.SNPs(idin,:);
    data.traits.n=arrayfun(@(i) sum(data.SNPs.id==i),1:data.nt)';
    newid=cumsum(data.traits.n~=0).*(data.traits.n~=0);
    data.SNPs.id=newid(data.SNPs.id);
    data.nt=data.nt-sum(data.traits.n==0);
    data.traits=data.traits(trts,:);
end

end
