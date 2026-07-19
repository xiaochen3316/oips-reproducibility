# Mathematical contract for OIPS cluster-v2

This document is normative for the released analysis. Configuration values are taken from `config/manuscript.yaml`; code paths are implemented in `clustering.py`, `scoring.py`, `posthoc.py`, and `statistics.py`.

## Conventions

| Symbol | Definition and range | Missing-value rule | Layer |
|---|---|---|---|
| `i`, `j` | Pocket records or same-tool units. | Not applicable. | Static. |
| `c` | A cluster-v2 candidate. | Not applicable. | Static. |
| `t` | A target (uppercase four-character PDB identifier). | Not applicable. | Both. |
| `R(x)` | Canonical residue set for object `x`, represented by chain, residue number, and insertion code; residue name is ignored for comparison. | The empty set is retained as empty. | Both. |
| `z(x)` | Three-dimensional pocket center in Å. | Unavailable unless all three coordinates are present. | Both. |
| `NA` | An unavailable numeric value, serialized as an empty CSV field. | Never treated as zero unless a formula below explicitly defines a fallback. | Both. |
| `1[condition]` | Indicator equal to 1 when the condition is true and 0 otherwise. Range `{0,1}`. | Not applicable. | Both. |
| `clamp_[0,100](x)` | `max(0,min(100,x))`. Range `[0,100]`. | Applied only to finite `x`. | Static. |

Unless stated otherwise, `≤` and `≥` thresholds are inclusive, whereas `>` and `<` audit thresholds are strict.

## Residue and center similarity

For nonempty canonical residue sets `A=R(i)` and `B=R(j)`, residue intersection-over-union and containment are

\[
J(A,B)=\frac{|A\cap B|}{|A\cup B|},
\qquad
K(A,B)=\frac{|A\cap B|}{\min(|A|,|B|)}.
\]

`J` and `K` lie in `[0,1]`. If either set is empty, both values are `NA`; the implementation does not define empty/empty similarity as 1. The center distance is

\[
d(i,j)=\lVert z(i)-z(j)\rVert_2
=\sqrt{\sum_{a\in\{x,y,z\}}(z_a(i)-z_a(j))^2},
\]

with range `[0,∞)` Å. It is `NA` if either center is unavailable. These functions are used in static clustering and again, without changing the static result, for post-hoc mappings.

## Same-tool compatibility

Let `T(i)` be the tool name and `H(i,j)` indicate that both records are DoGSiteScorer or DoGSite3 records with the same `P_<number>` hierarchy root. Two raw records are same-tool duplicates exactly when `T(i)=T(j)` and at least one branch below is true:

\[
\begin{aligned}
D_{same}(i,j)={}&H(i,j)\land
  [d\le6\ \lor\ J\ge0.35\ \lor\ K\ge0.70]\\
&\lor J\ge0.70\\
&\lor [K\ge0.85\land(d\text{ is NA}\ \lor\ d\le8)]\\
&\lor [d\le2.5\land J\ge0.25]\\
&\lor [R(i)=R(j)=\varnothing\land d\le1.5].
\end{aligned}
\]

Here `d=d(i,j)` in Å and `J=J(R(i),R(j))`, `K=K(R(i),R(j))`. A comparison involving `NA` is false except for the explicit `d is NA` allowance in the containment branch. Every displayed bound is inclusive.

A record may join same-tool group `G` only if

\[
\forall j\in G:\ D_{same}(i,j)=1,
\qquad
\Delta(G\cup\{i\})\le8\ \text{Å}.
\]

This is a complete-link duplicate requirement plus an inclusive group-diameter constraint. `D_same` is Boolean. The group diameter `Δ` is defined below. Missing centers are omitted from `Δ`; zero or one available center gives `Δ=0`.

## Cross-tool compatibility and diameter cap

For units from different tools, the pairwise compatibility predicate is

\[
D_{cross}(i,j)=
1[d\le6]
\lor 1[d\le10\land J\ge0.20]
\lor 1[J\ge0.35].
\]

`D_cross` is Boolean. `d` is in Å and `J` is in `[0,1]`; unavailable terms make only their own branch false. The 6, 10, 0.20, and 0.35 thresholds are inclusive. Cross-tool candidate assignment requires at least one compatible link to the existing candidate and a valid medoid-constrained group.

For the subset `Z_c` of available unit centers in candidate `c`, the center diameter is

\[
\Delta_c=
\begin{cases}
0, & |Z_c|\le1,\\
\max_{u,v\in Z_c}\lVert u-v\rVert_2, & |Z_c|\ge2.
\end{cases}
\]

`Δ_c` has range `[0,∞)` Å before validation. A formal candidate must satisfy

\[
\Delta_c\le12.0\ \text{Å},
\]

inclusively, and every non-medoid unit must satisfy the same site-pair predicate with the selected medoid. A value greater than 12.0 Å is invalid, not merely boundary-sensitive. Missing centers do not by themselves fail the diameter cap; residue compatibility and mappability remain explicit.

The deterministic similarity used only to choose between otherwise valid assignments is the mean of available terms

\[
S(i,j)=\operatorname{mean}_{NA}\left(
\max\left(0,1-\frac{d(i,j)}{12}\right),
J(R(i),R(j))
\right),
\]

with range `[0,1]` and value 0 if both terms are unavailable. `S` chooses among compatible groups; it does not replace `D_cross` or the 12 Å cap.

## Residue consensus and boundary audit

Let `n_c` be the number of formal tool representatives in cluster `c`, with range `1–5` in this five-tool release. For canonical residue `r`, let `s_c(r)` be its representative support count, range `1–n_c`. The envelope and core are

\[
E_c=\{r:s_c(r)\ge1\},
\qquad
C_c=\{r:s_c(r)\ge\tau_c\},
\]

where

\[
\tau_c=
\begin{cases}
1, & n_c=1,\\
\max(2,\lceil n_c/2\rceil), & n_c>1.
\end{cases}
\]

`E_c` and `C_c` may be empty. The core/envelope ratio `q_c=|C_c|/|E_c|` lies in `[0,1]` and is `NA` when `E_c` is empty.

For `n_c≥2`, `boundary_sensitive` is true when any of the following strict conditions holds:

\[
\Delta_c>9.0,
\quad \sigma_c>4.0,
\quad q_c<0.20,
\quad \widetilde{J}_c<0.15.
\]

`σ_c` is the root-mean-square distance, in `[0,∞)` Å, from available formal-representative centers to the medoid center; it is `NA` when the medoid center is missing. `J̃_c` is the median of available formal-representative pairwise IoUs, in `[0,1]`, and is `NA` if none is available. An unavailable audit term does not trigger its condition. For `n_c<2`, the boundary flag is false. This audit is static and never removes a candidate.

## Static module scores

### Consensus module `C_cons`

Let `V_c` be the set of tools contributing a formal representative to cluster `c`, let `V_all` be the set of all five configured tools, and let the configured tool weights be

\[
w^{tool}=\{1.00,1.00,1.00,0.75,1.15\}
\]

for DoGSiteScorer, DoGSite3, CavityPlus, CASTpFold, and SiteMap, respectively. Then

\[
C_{cons,c}=100
\frac{\sum_{u\in V_c}w^{tool}_u}{\sum_{u\in V_{all}}w^{tool}_u}
\times
\begin{cases}
0.5,& |V_c|=1,\\
1,& |V_c|>1.
\end{cases}
\]

`C_cons` lies in `[0,100]`. A valid candidate has at least one formal representative. A nonpositive total configured tool weight is invalid configuration rather than a missing score. This is a static consensus measure, not a post-hoc evidence count.

### Geometry and ligandability modules

Let `g_{cu}` and `p_{cu}` be the available geometry and ligandability values of formal representative `u` in cluster `c`. Let `Top_m(X)` return the `min(m,|X|)` largest values of nonempty finite multiset `X`. Then

\[
G_{geo,c}=\operatorname{mean}(Top_3(\{g_{cu}:g_{cu}\text{ available}\})),
\]

\[
P_{lig,c}=\operatorname{mean}(Top_2(\{p_{cu}:p_{cu}\text{ available}\})).
\]

Each output lies between the minimum and maximum of its selected finite input values. The public schema requires finite numbers when present but does not impose a numeric bound, so the formal code range is finite real numbers rather than a guaranteed `[0,100]`. If no corresponding value is available, the module is `NA`. These are static detector-derived modules.

### Static evidence-quality module `Q_evidence`

Let `n_c` be formal-representative count, `f^z_c` the fraction with complete centers, `f^R_c` the fraction with nonempty residue sets, and `I^{SM}_c` indicate that SiteMap contributes. Their ranges are `n_c∈{1,…,5}`, `f^z_c,f^R_c∈[0,1]`, and `I^{SM}_c∈{0,1}`. Then

\[
Q_{evidence,c}=\operatorname{clamp}_{[0,100]}
\left(30+9n_c+18f^z_c+14f^R_c+4I^{SM}_c\right).
\]

The result lies in `[0,100]`. Zero representatives would produce `NA`, but such a candidate is invalid and is not released. Missing centers or residues lower their fractions; they are not imputed. This static module measures representation completeness only.

### Oligomer-relevance module `O_rel_formal`

The interface set `I_t` contains standard-amino-acid residues with an inter-chain heavy-atom contact at distance at most 5.0 Å. Let `P_c=R(c)` be the cluster envelope after canonicalization. Define

\[
f_c=\frac{|P_c\cap I_t|}{|P_c|},
\qquad
r_c=\frac{|P_c\cap I_t|}{|I_t|}.
\]

`f_c` (interface fraction) and `r_c` (interface recall) lie in `[0,1]`. `f_c` is `NA` when `P_c` is empty; `r_c` is `NA` when `I_t` is empty. Recall is reported but is not a direct `O_rel_formal` subscore.

Let `k_c` be the number of chains represented in the unique full cluster residue identifiers, range `0–∞`. If `n_{ch}` chains carry proportions `p_h` of the unique cluster residues, normalized chain entropy is

\[
H_c=
\begin{cases}
0,& n_{ch}\le1,\\
-\dfrac{\sum_{h=1}^{n_{ch}}p_h\log p_h}{\log n_{ch}},& n_{ch}>1.
\end{cases}
\]

`H_c∈[0,1]`; empty and single-chain residue sets give 0.

Let `d^I_c` be the minimum Euclidean distance from the cluster medoid center to any interface-atom coordinate, range `[0,∞)` Å and `NA` if either is unavailable. Its score is the inclusive piecewise contract

\[
D_c=
\begin{cases}
100,& d^I_c\le4,\\
90-5(d^I_c-4),& 4<d^I_c\le8,\\
70-6(d^I_c-8),& 8<d^I_c\le16,\\
\max(10,22-1.2(d^I_c-16)),& d^I_c>16,\\
NA,& d^I_c\text{ unavailable}.
\end{cases}
\]

The implemented boundary at exactly 4 Å is 100; the next branch begins strictly above 4 Å at values just below 90. `D_c` otherwise lies in `[10,100]` when present.

The other three subscores are

\[
F_c=\operatorname{clamp}_{[0,100]}(25+130f_c^*),
\qquad
f_c^*=\begin{cases}f_c,&f_c\text{ available},\\0,&f_c\text{ is NA},\end{cases}
\]

\[
L_c=
\begin{cases}
75+25H_c,& k_c\ge2,\\
65,& k_c<2\land f_c\ge0.35,\\
30,& \text{otherwise},
\end{cases}
\]

and

\[
X_t=\operatorname{clamp}_{[0,100]}\left(
45+10\min(4,K_t-1)+80\min(0.5,e_t)
\right),
\qquad
e_t=\begin{cases}|I_t|/N^{prot}_t,&N^{prot}_t>0,\\0,&N^{prot}_t=0.\end{cases}
\]

Here `F_c∈[0,100]` is the interface-fraction score; `L_c∈[30,100]` is chain-context support; `K_t` is the number of protein chains in the prepared structure; `N^{prot}_t` is its number of unique standard-amino-acid residues; `e_t∈[0,1]` is interface extent; and `X_t∈[0,100]` is the interface-extent/context score. The `f_c≥0.35` test is inclusive. Only `D_c` can be missing in the normal multi-chain/interface branch.

For a structure with at least two protein chains and a nonempty detected interface,

\[
O_{rel,c}=\operatorname{WMean}_{NA}
\left((F_c,L_c,D_c,X_t),(0.42,0.25,0.20,0.13)\right).
\]

For fewer than two protein chains, `O_rel,c=35`. For at least two chains but no detected interface residues, `O_rel,c=45`. `O_rel,c`, serialized as `O_rel_formal`, lies in `[0,100]`. `WMean_NA` removes unavailable terms and renormalizes their positive weights. The module is static and uses only the prepared assembly proxy.

## Missing-aware `OIPS-P_static`

Let the five static modules be

\[
M=\{C_{cons},G_{geo},P_{lig},O_{rel},Q_{evidence}\}
\]

with weights

\[
(w_C,w_G,w_P,w_O,w_Q)=(0.22,0.18,0.24,0.24,0.12).
\]

For candidate `c`, let `M_{m,c}` be the numeric value of module `m`, let `w_m` be its declared weight, and let `A_c⊆M` be the modules with finite values. Then

\[
OIPS\text{-}P_{static,c}=
\frac{\sum_{m\in A_c}w_mM_{m,c}}
     {\sum_{m\in A_c}w_m}.
\]

Missing modules are excluded from both numerator and denominator; they are not zero-filled. The result is `NA` if `A_c` is empty. When every available module lies in `[0,100]`, the score also lies in `[0,100]`; because the public schema does not bound `G_geo` or `P_lig`, the formal code range is finite real values. The released data use the declared detector score scales.

Within target `t`, ranks are descending:

\[
rank_{tc}=1+|\{j:OIPS_j>OIPS_c\}|+
\frac{|\{j:OIPS_j=OIPS_c\}|-1}{2}.
\]

Here `j` ranges over candidates of target `t`, and `OIPS_j` abbreviates `OIPS-P_static` for candidate `j`. This is average tie ranking. The range is `[1,N^{cand}_t]`, potentially fractional, where `N^{cand}_t` is the target's candidate count. An unavailable static score would have an unavailable rank and fails the released bundle contract. Post-hoc evidence never changes this rank or breaks ties.

## One-at-a-time module-weight sensitivity

For module `k` and direction `s∈{-1,+1}`, define the perturbation multiplier `a_s` as 0.8 for a decrease and 1.2 for an increase. The scenario weights are

\[
w_m^{(k,s)}=
\begin{cases}
a_s w_m,&m=k,\\
w_m,&m\ne k.
\end{cases}
\]

For candidate `c`, the perturbed score retains the baseline availability set `A_c` and renormalizes the available scenario weights:

\[
OIPS_c^{(k,s)}=
\frac{\sum_{m\in A_c}w_m^{(k,s)}M_{m,c}}
     {\sum_{m\in A_c}w_m^{(k,s)}}.
\]

Candidate membership, module values, missingness, and post-static labels are fixed. Perturbed within-target ranks use the same descending average-tie definition as the baseline. Let `B_t(1)` and `S_t^{(k,s)}(1)` be the baseline and scenario Top-1 identifier sets, and let `B_t(3)` and `S_t^{(k,s)}(3)` contain candidates with rank at most 3. Scenario summaries include

\[
Retention_{Top1}^{(k,s)}=\sum_t 1[B_t(1)=S_t^{(k,s)}(1)],
\]

\[
Jaccard_{Top3}^{(k,s)}=
\frac{1}{|T|}\sum_t
\frac{|B_t(3)\cap S_t^{(k,s)}(3)|}{|B_t(3)\cup S_t^{(k,s)}(3)|},
\]

and the median across targets of the Spearman correlation between baseline and scenario candidate-rank vectors. Reference and first-supported scenario ranks are the minimum perturbed ranks among the fixed `R_auto` and fixed `R_auto/A_auto` candidate sets, respectively. Their Top-k counts use inclusive rank thresholds.

## Reference and MD overlap

For any nonempty candidate residue set `P` and evidence-contact set `E`, define

\[
n_{PE}=|P\cap E|,
\quad
precision(P,E)=\frac{n_{PE}}{|P|},
\quad
recall(P,E)=\frac{n_{PE}}{|E|},
\quad
IoU(P,E)=\frac{n_{PE}}{|P\cup E|}.
\]

The count is a nonnegative integer and the three ratios lie in `[0,1]`. If either set is empty, the count is 0 and all three ratios are `NA`. Candidate residues and evidence residues use the same canonicalization. These are post-hoc quantities.

### Reference association

The reference-contact set contains standard-amino-acid residues having an `ATOM` coordinate within 4.5 Å, inclusively, of any atom in the selected ligand. `DCC_c` is Euclidean distance from candidate medoid center to ligand-atom centroid, in `[0,∞)` Å and `NA` when either center is unavailable. Let `r_c` be reference contact recall and `J_c` reference residue IoU. The automated reference rule is

\[
R_{auto,c}=
1[DCC_c\le6\land r_c\ge0.10]
\lor 1[DCC_c\le10\land(r_c\ge0.20\lor J_c\ge0.15)]
\lor 1[DCC_c\le12\land r_c\ge0.50\land J_c\ge0.15].
\]

All bounds are inclusive. Missing `DCC` or recall makes the rule false. Missing IoU makes only IoU-dependent branches false. An unresolved ligand selection or clear QC exclusion prevents the final `R_auto` evidence label even if this geometric predicate is true.

### MD mapping and concordance

For MD run `q`, `E_q` is the released persistent-contact set and `z_q` is its optional center. The same overlap definitions produce candidate coverage (`precision`), MD-contact coverage (`recall`), and Jaccard (`IoU`). Center distance `d^q_c=‖z(c)-z_q‖_2` lies in `[0,∞)` Å or is `NA`.

When `E_q` is nonempty, the best mapped candidate maximizes IoU, then overlap count, then minimizes center distance, then identifier. Let subscript `1` denote the immutable static-top candidate and `b` the best MD-mapped candidate. The ordered call contract is:

1. `concordant` if `J_1≥0.15` or (`precision_1≥0.20` and `d^q_1≤6` Å).
2. Otherwise `partially_concordant` if `J_1≥0.05`.
3. Otherwise `static_dynamic_conflict` if `b≠1` and (`J_b≥0.05` or `d^q_b≤6` Å).
4. Otherwise `boundary_shift` if `d^q_1≤8` Å.
5. Otherwise `static_dynamic_conflict`.

All stated thresholds are inclusive; `NA` comparisons are false. If `E_q` is empty, an apo simulation is `apo_only_context`; another context is `insufficient_MD_evidence`. A target without a released run is `MD_not_available`. `D_dyn_run_score` is reported as a finite real value or `NA` but is absent from every static formula.

## Top-k and MRR endpoints

For target `t`, let `r^R_t` be the smallest static rank among `R_auto` candidates. Let `E_R` be the set of targets with a resolved, nonempty reference ligand. When no `R_auto` candidate exists for an evaluable target, set `r^R_t=∞` for endpoint arithmetic. Reference Top-k and MRR are

\[
TopK_R(k)=\frac{1}{|E_R|}\sum_{t\in E_R}1[r^R_t\le k],
\qquad k\in\{1,3,5\},
\]

\[
MRR_R=\frac{1}{|E_R|}\sum_{t\in E_R}
\begin{cases}
1/r^R_t,&r^R_t<\infty,\\
0,&r^R_t=\infty.
\end{cases}
\]

`TopK_R` and `MRR_R` lie in `[0,1]` and are `NA` if `E_R` is empty. The `≤k` threshold is inclusive.

Let `r^S_t` be the smallest static rank among candidates labeled `R_auto` or `A_auto`, and let `E_S={t:r^S_t<∞}`. First-supported endpoints are

\[
TopK_S(k)=\frac{1}{|E_S|}\sum_{t\in E_S}1[r^S_t\le k],
\qquad k\in\{1,3\},
\]

\[
MRR_S=\frac{1}{|E_S|}\sum_{t\in E_S}\frac{1}{r^S_t}.
\]

They lie in `[0,1]` and are `NA` if `E_S` is empty. `E_R` and `E_S` are intentionally non-equivalent denominators.

## Native single-tool complete-case comparison

For tool `u`, let `U_{tu}` be the same-tool units for target `t`, ordered by the retained native keys `display_order`, then `sitemap_rank`, then `row_id`, with missing order values last. Each unit is mapped to its cluster-v2 identifier through the exhaustive source crosswalk. If multiple ordered units map to one cluster, only the earliest is retained. The resulting regional rank `q_{tuc}` is the one-based position of mapped cluster `c` in this deduplicated native sequence.

The single-tool first-reference rank is

\[
r^{u}_t=\min\{q_{tuc}:R_{auto,tc}=1\},
\]

and is infinite when no mapped region passes `R_auto`. The complete-case target set is

\[
E_{CC}=\{t: \forall u,\ |U^{mapped}_{tu}|\ge1\}.
\]

For this release, `|E_CC|=14`. For `k∈{1,3,5}`, each tool's complete-case endpoint and MRR are

\[
TopK_u(k)=\frac{1}{|E_{CC}|}\sum_{t\in E_{CC}}1[r^u_t\le k],
\qquad
MRR_u=\frac{1}{|E_{CC}|}\sum_{t\in E_{CC}}
\begin{cases}1/r^u_t,&r^u_t<\infty,\\0,&r^u_t=\infty.\end{cases}
\]

OIPS-P uses its immutable cluster-v2 static rank in the same equations and on the same target set. A no-hit target contributes zero, not a missing value.

For a finite reference rank and `N^{cand}_t` candidates, the reported rank percentile is

\[
Pctl_t=
\begin{cases}
1,&N^{cand}_t\le1,\\
1-\dfrac{r^R_t-1}{N^{cand}_t-1},&N^{cand}_t>1.
\end{cases}
\]

It lies in `[0,1]` and is `NA` without a finite reference rank.

## Bootstrap intervals

Let `\theta(D)` be one of the seven metrics `Reference_Top1`, `Reference_Top3`, `Reference_Top5`, `Reference_MRR`, `First_supported_Top1`, `First_supported_Top3`, or `First_supported_MRR` computed from target table `D`.

For target resampling, each replicate `b=1,…,10000` draws `N=21` target rows independently with replacement:

\[
D_b^{target}=\{T_{b1},\ldots,T_{bN}\},
\qquad T_{bj}\sim Uniform(D).
\]

For family-clustered resampling, let `F={f_1,…,f_L}` be the `L=13` family labels and `D_f` all targets in family `f`. Each replicate draws `L` labels with replacement and concatenates their full target groups:

\[
D_b^{family}=D_{F_{b1}}\uplus\cdots\uplus D_{F_{bL}},
\qquad F_{bj}\sim Uniform(F).
\]

Each method initializes an independent NumPy random generator with seed `20260710`. For method `m`, the interval is

\[
CI_m(\theta)=
[Q_{0.025}(\{\theta(D_b^m)\}),
 Q_{0.975}(\{\theta(D_b^m)\})].
\]

`Q_p` is the nan-aware empirical quantile used by NumPy. Metric and interval endpoints lie in `[0,1]`; a replicate with no evaluable denominator yields `NA` and is ignored by `Q_p`. The reported point estimate is `θ(D)`, not the bootstrap mean. All resampling and interval calculations are post-hoc.

## Family sensitivity

For each observed family `f∈F`, leave-one-family-out data are

\[
D_{-f}=\{t\in D:family(t)\ne f\}.
\]

The complete endpoint vector is recomputed as `θ(D_{-f})`. `excluded_target_N=|D|-|D_{-f}|` and `remaining_target_N=|D_{-f}|` are nonnegative integers. Ratio metrics retain `[0,1]`; ranks lie in `[1,∞)` when available and counts are nonnegative. Empty endpoint denominators produce `NA`. This analysis does not fit parameters or impute a removed family.

## `O_rel` ablation

For candidate `c`, let `A_c^{-O}` contain the available modules among `C_cons`, `G_geo`, `P_lig`, and `Q_evidence`, with the original weights `(0.22,0.18,0.24,0.12)`. The ablated score is

\[
Score^{-O}_c=
\frac{\sum_{m\in A_c^{-O}}w_mM_{m,c}}
{\sum_{m\in A_c^{-O}}w_m},
\]

If none is available, the score is `NA`; otherwise its range follows the available modules and is `[0,100]` when they are all in that range. `O_rel_formal` alone is removed. All other modules, candidate membership, labels, and original static ranks are held fixed.

`Without_O_rel_rank` is the descending average-tie rank of `Score^{-O}` within target. For any endpoint rank `r`, the reported rank change is

\[
\Delta r=r^{-O}-r^{full},
\]

with finite real range `[-(N^{cand}_t-1),N^{cand}_t-1]` when both ranks exist and `NA` otherwise. Positive values indicate worsening after removal. Top-3 inclusion change is

\[
\Delta Top3=1[r^{-O}\le3]-1[r^{full}\le3]\in\{-1,0,1\}.
\]

`top_cluster_identity_changed` is Boolean. Per-target `rank_Spearman_rho` is the Spearman correlation between the full and ablated candidate-rank vectors, range `[-1,1]`; it is `NA` if either vector has fewer than two distinct finite ranks. Category summaries take the arithmetic mean of available target correlations and identity-change indicators and recompute reference Top-1, Top-3, and MRR on reference-evaluable targets. This is a post-hoc sensitivity analysis and does not alter the released static ranking.
