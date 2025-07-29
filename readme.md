Use:
infer_this_data_final('datafile.mat','Resultsfolder/resultfilename.mat','SSD')

or

infer_this_data_final('datafile.mat','Resultsfolder/resultfilename.mat','TSD')


Options:
boot 		bootstrap number (0 indicates no bootstrap)
jack 		jackknife number (-1 indicates no jackknife)
eps		penalty strength
k		no. of knots minus one
minx		minimal included MAF
maxld		maximal included LD score
minx		minimal included MAF
minsnps	minimal number of hits per trait
ex		excluded traits
splinefile	file with precalculated normalizing splines
tablefile	file with precalculated prob of inclusion at the LD threshold for different MAF bins
funevals	maximal number of likelihood function evaluations
maxgrid	location of knot with largest log10(s)
mingrid	location of knot with smallest log10(s)
maxlogs	maximal included log10(s)
minlogs	minimal included log10(s)
targetsize	estimated size of functional genome for the penalty
pout		Proportion of outliers SNPs (set to 0 for inference without accounting for outliers)
zout		mean z-score of outliers


 
Inference Files:
infer_this_data.m – main inference file
do_inference.m – the simulated annealing maximization
jack_this_data.m/boot_this_data.m – jackknife/bootstrap data
preprocess.m – pre process the data, that is remove variants based on MAF, z and LD score cutoff, and traits with less than 100 hits
get_pres.m – get the residual p-values
make_phittable.m  - get the probability of variant i being a hit given parameter C and study sizes                                                                                                                                                                  
make_ptable.m - get the probability of seeing effect size b and MAF x at variant i being a hit given parameter C and study sizes        
subset_traits.m – choose a subset of traits from data                                                                                                                                                          

Prep files:
make_datafile_from_cojo.m  - make data file from cojo files                                                                                                                                                         
make_normalizing_integrals_stage1.m  - makes auxiliary table, file tables.mat
make_normalizing_integrals_stage2.m  - uses tables.mat to make the normalizing splines, needs to be run for each selection coefficient separately
collate_splines.m – collates the normalizing splines for different selection coefficients into splines.mat

Data files:
data96.mat - data file
data15.mat - data file with subset of 15 nearly-independent traits
forward_sim_data.mat - data file with results of forward simulations (as described in Simons 2018)
splines.mat  - normalizing splines for LD=500
tables.mat  - prob of inclusion for different MAF bins and LD=500 (+ other tables)
dens.mat – table with probability density of MAFs, i.e. p(MAF in bin|s)/size(MAF bin) – needed data from forward_sim_data.mat exported to a smaller file to save memory
shtot.csv – s*mean heterozygosity for all s
blocks.tsv – approximately independent genomic blocks

 
Data format:
Data.SNPs – table with columns: 
chr (chromosome), SNP (chr:bp:ref:alt),bp (bsepair), refA (reference allele), x (MAF), b (effect size), se (standard error for effect size), p (p-value pre COJO), n (reported sample size), freq_geno (frequency in COJO LD reference sample), bJ (joint effect from cojo), bJ_se (standard error for joint effect), pJ (p-value for joint effect), LD_r (r for LD between this SNP and the next in COJO), id (trait id), m_rel (relative study size from x and se), mJ (reduction in study size in bJ_se), loc (1e9*chr+bp – a unique identifier number), xb (MAF bin of x), blk (genome block), ld (LD score), a (absolute value of bJ), a_se, (absolute value of bJ_se), z (a/a_se)), v (2*a^2*x*(1-x)), dv (2*a_se^2*x*(1-x))
Data.traits – table with columns: 
trait (trait code in UKbb), n (number of included hits), M_med (median effective study size), ReportedM_med (median reported study size), name (common trait name)
Data.nt – number of traits


Results format:

bootidx – indices of included SNPs in the bootstrap/jackknife
outidx – indices of excluded SNPs in the jackknife
pres – residual p-values of excluded SNPs in the jackknife
inf.LL – log likelihood (1 x no. of traits)
inf.LL0 – log-likelihood of initial guess (1 x no. of traits)

inf.c –C parameter (1 x no. of traits)
inf.L –estimated target size (1 x no. of traits)
inf.V –estimated heritability (1 x no. of traits)
inf.yout – parameters of spline – size: (no. of knots-1) x 1[SSD]; (no. of knots-1) x (no. of traits) [TSD] 
inf.p –  inferred distribution of selection coefficients, normalized to 1 - size: (size of log(s) vector) x 1[SSD]; (size of log(s) vector) x (no. of traits) [TSD] 
inf.logs –  log10(s) vector used
inf.mingrid –  location of lowest knot
inf.maxgrid –  location of highest knot
inf.mu – mean log10(selection coefficient) - size: 1 x 1[SSD]; 1 x (no. of traits) [TSD] 
inf.sig – standard deviation of log10(selection coefficient) - size: 1 x 1[SSD]; 1 x (no. of traits) [TSD] 
inf.Pout – probability that a SNP is an outlier - size: (no. of hits) x 1 [SSD]; (no. of traits) vectors of size (no. of hits in trait) x 1 each [TSD]

inf.ys, inf.ps, inf.LLlst, inf.acceptedlst, inf.switchlst, inf.Tlst – statistics from the maximization process

