# twinSEM
twin modeling on package OpenMx 

### what's the data?
- three groups: MZ, DZ/bio sibs, Adoptees
- rater items: conflict (con), involvement (inv), regard (preg)
  - running most bivariates w/ con/inv only since regard is wack
- each item measured by age (14/17), kid relationship with mom/dad, rater (self, parent, parent other parent)
  - ie con14ms is conflict with mom at 14 rated by kid, inv17frm is mom involvement with kid at 17 rated by father
  - items have large skew/kurtosis so they've been log transformed after control regressions
- bpd is outcome; at age 14 (bpd14) age 17 (bpd17)
  - nicely normally distributed so no log transform

control regressions:
- regressed out age and sex and saved resids 
  - whether log(item) ~ age + sex or item ~ age + sex depends on above
- items with 'r' after name (ie con14msr) means they've had age/sex regressed out

### model fun stuff:
- possibile to fit ACE (two/three), ADE (two/three), ACDE (three) 
  - most scripts contain models with all these structures + submodels
  - submodels ACDE: ACE, ADE, CDE (not feasible), AE, CE, DE, E

all fitted models thus far:
- saturated model w/ submodels
- univariate AC(D)E (cholesky/VC)
- moderated univariate ACE (cholesky)
- sib interaction AC(D)E (VC)
- bivariate AC(D)E three (cholesky/VC)
- rater bias model ACE (cholesky/VC)

to be added:
- common pathway ACE for three groups
- common pathway ACE with rater bias paths
- bivariate sib interaction

note: named functions 'path' are just cholesky decompositions



