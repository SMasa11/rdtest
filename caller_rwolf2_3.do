rwolf2 (rdrobust z1 x) (rdrobust z2 x) (rdrobust z3 x), indepvars(RD_Estimate,RD_Estimate,RD_Estimate) seed(7628775) nodots usevalid
qui clear
qui set obs 1
qui gen rw_pval_1 = e(rw_z1_RD_Estimate)
qui gen rw_pval_2 = e(rw_z2_RD_Estimate)
qui gen rw_pval_3 = e(rw_z3_RD_Estimate)