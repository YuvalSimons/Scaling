# Usage

```matlab
infer_this_data_final('datafile.mat','Resultsfolder/resultfilename.mat','SSD')
```
or
```matlab
infer_this_data_final('datafile.mat','Resultsfolder/resultfilename.mat','TSD')
```

---

# Options

| Option       | Description                                                                 |
|--------------|-----------------------------------------------------------------------------|
| `boot`       | Bootstrap number (0 indicates no bootstrap)                                 |
| `jack`       | Jackknife number (-1 indicates no jackknife)                                |
| `eps`        | Penalty strength                                                            |
| `k`          | No. of knots minus one                                                      |
| `minx`       | Minimal included MAF                                                        |
| `maxld`      | Maximal included LD score                                                   |
| `minsnps`    | Minimal number of hits per trait                                            |
| `ex`         | Excluded traits                                                             |
| `splinefile` | File with precalculated normalizing splines                                 |
| `tablefile`  | File with precalculated prob of inclusion at the LD threshold for MAF bins  |
| `funevals`   | Maximal number of likelihood function evaluations                           |
| `maxgrid`    | Location of knot with largest `log10(s)`                                    |
| `mingrid`    | Location of knot with smallest `log10(s)`                                   |
| `maxlogs`    | Maximal included `log10(s)`                                                 |
| `minlogs`    | Minimal included `log10(s)`                                                 |
| `targetsize` | Estimated size of functional genome for the penalty                         |
| `pout`       | Proportion of outlier SNPs (set to 0 for inference without outliers)        |
| `zout`       | Mean z-score of outliers                                                    |

---

# Inference Files

- `infer_this_data.m` – main inference file  
- `do_inference.m` – simulated annealing maximization  
- `jack_this_data.m` / `boot_this_data.m` – jackknife/bootstrap data  
- `preprocess.m` – removes variants based on MAF, z, LD score cutoff, and traits with <100 hits  
- `get_pres.m` – get residual p-values  
- `make_phittable.m` – probability variant *i* is a hit given `C` and study sizes  
- `make_ptable.m` – probability of observing effect size *b*, MAF *x*, at variant *i*  
- `subset_traits.m` – choose subset of traits from data  

---

# Prep Files

- `make_datafile_from_cojo.m` – make data file from COJO files  
- `make_normalizing_integrals_stage1.m` – generates auxiliary table `tables.mat`  
- `make_normalizing_integrals_stage2.m` – uses `tables.mat` to make splines per selection coefficient  
- `collate_splines.m` – collates splines for different selection coefficients into `splines.mat`  

---

# Data Files

- `data96.mat` – full data file  
- `data15.mat` – subset with 15 nearly-independent traits  
- `forward_sim_data.mat` – forward simulation results (Simons 2018)  
- `splines.mat` – normalizing splines for LD=500  
- `tables.mat` – inclusion probabilities by MAF bin and LD=500  
- `dens.mat` – `p(MAF in bin | s) / size(MAF bin)`  
- `shtot.csv` – `s × mean heterozygosity` for all `s`  
- `blocks.tsv` – approximately independent genomic blocks  

---

# Data Format

### `Data.SNPs` – table with columns:
- `chr`: Chromosome  
- `SNP`: `chr:bp:ref:alt`  
- `bp`: Basepair position  
- `refA`: Reference allele  
- `x`: MAF  
- `b`: Effect size  
- `se`: Standard error for effect size  
- `p`: p-value (pre-COJO)  
- `n`: Reported sample size  
- `freq_geno`: Frequency in COJO LD reference sample  
- `bJ`: Joint effect from COJO  
- `bJ_se`: Standard error for joint effect  
- `pJ`: p-value for joint effect  
- `LD_r`: LD between this SNP and next in COJO  
- `id`: Trait ID  
- `m_rel`: Relative study size from `x` and `se`  
- `mJ`: Reduction in study size in `bJ_se`  
- `loc`: `1e9 * chr + bp` (unique identifier)  
- `xb`: MAF bin  
- `blk`: Genome block  
- `ld`: LD score  
- `a`: |`bJ`|  
- `a_se`: |`bJ_se`|  
- `z`: `a / a_se`  
- `v`: `2 * a^2 * x * (1 - x)`  
- `dv`: `2 * a_se^2 * x * (1 - x)`  

### `Data.traits` – table with columns:
- `trait`: Trait code (UKBB)  
- `n`: Number of included hits  
- `M_med`: Median effective study size  
- `ReportedM_med`: Median reported study size  
- `name`: Common trait name  

### `Data.nt` – number of traits  

---

# Results Format

- `bootidx`: Indices of included SNPs in bootstrap/jackknife  
- `outidx`: Indices of excluded SNPs in jackknife  
- `pres`: Residual p-values of excluded SNPs  
- `inf.LL`: Log likelihood (1 × number of traits)  
- `inf.LL0`: Log likelihood of initial guess (1 × number of traits)  
- `inf.c`: `C` parameter (1 × number of traits)  
- `inf.L`: Estimated target size (1 × number of traits)  
- `inf.V`: Estimated heritability (1 × number of traits)  
- `inf.yout`: Spline parameters  
  - SSD: (`knot count - 1`) × 1  
  - TSD: (`knot count - 1`) × number of traits  
- `inf.p`: Inferred distribution of selection coefficients, normalized  
  - SSD: (log(s) vector size) × 1  
  - TSD: (log(s) vector size) × number of traits  
- `inf.logs`: `log10(s)` vector used  
- `inf.mingrid`: Location of lowest knot  
- `inf.maxgrid`: Location of highest knot  
- `inf.mu`: Mean `log10(selection coefficient)`  
  - SSD: 1 × 1  
  - TSD: 1 × number of traits  
- `inf.sig`: Std dev of `log10(selection coefficient)`  
  - SSD: 1 × 1  
  - TSD: 1 × number of traits  
- `inf.Pout`: Outlier probability per SNP  
  - SSD: (number of hits) × 1  
  - TSD: One vector per trait: (number of hits in trait) × 1  
- `inf.ys`, `inf.ps`, `inf.LLlst`, `inf.acceptedlst`, `inf.switchlst`, `inf.Tlst`: Maximization stats  

---
