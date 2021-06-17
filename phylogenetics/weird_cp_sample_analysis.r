# Analysis to make sense of the lowland samples with upland chloroplast
#   haplotypes

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD LIBRARIES ###
library(data.table)

### INPUT DATA ###
res_tab_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt'
res_tab <- fread(res_tab_file)

meta_file <- '/home/grabowsky/tools/workflows/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

#weird_lib_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_first_draft_analysis/LIB_list_of_gulf_atlantic_admixed_cluster_in_upland.csv'
weird_lib_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_first_draft_analysis/List_of_upland_cp_LIBs.csv' 
weird_libs <- fread(weird_lib_file)

##########

weird_samps <- samp_meta[which(samp_meta$LIB %in% weird_libs$LIB), VCF_NAME]

res_weird_inds <- which(res_tab$samp_name %in% weird_samps)

res_tab[res_weird_inds, ploidy]

res_tab[res_weird_inds, .N, by = subgrp_v2]
res_tab[res_weird_inds, .N, by = ploidy]

res_tab[res_weird_inds, full_admix_k3_MW]

res_tab[res_weird_inds, list(subgrp_v2, full_admix_k3_MW)]

samp_meta[which(samp_meta$LIB %in% weird_libs$LIB), .N, by = STATE]
samp_meta[which(samp_meta$LIB %in% weird_libs$LIB), .N, by = UNI_ACC]
# 163
table(table(samp_meta[which(samp_meta$LIB %in% weird_libs$LIB), UNI_ACC]))
#  1  2  3  4  5  6 
# 75 55 26  4  1  2 

as.data.table(sort(weird_samps))

# The first set of samples...
######
 1:  J019.A - Atlantic, FL; ATL_02
 2: J019.L1 - Atlantic, FL
 3: J019.L2 - Atlantic, FL
 4: J019.L3 - Atlantic, FL
 5:  J177.A - Gulf(75%)/ATL admix, FL; GULF_12
 6: J179.L2 - Gulf(72%)/ATL admix; FL; GULF_12
 7: J179.L3 - Gulf(85%)/ATL admix; FL; GULF_12
 8:  J189.A - Gulf(75%)/ATL admix; FL; GULF_09
 9:  J209.A - Gulf, TX; GULF_06 - NOTE: GULF_06 includes 8X samples with 
## some MW ancestry; this clade also includes a weird 6X sample from FL that
## GRIN says is P. amarum...
10:  J209.B - Gulf, TX; GULF_06 
11:  J228.A - Gulf, TX; GULF_06
12:  J240.A - Gulf, TX; GULF_06
13:  J294.A - Gulf, TX: GULF_06
14:  J303.A - Gulf, TX: GULF_06
15:  J305.A - Gulf, TX: GULF_06
16:  J312.A - Gulf, TX: GULF_07, same NJ small clade as GULF_06
17:  J313.A - Gulf, TX: GULF_07, same NJ small clade as GULF_06
18:  J461.A - Gulf (94%) / ATL, MS; GULF_03
19:  J461.B - Gulf (95%) / ATL, MS; GULF_03
20:  J461.C - Gulf (97%) / ATL, MS; GULF_03
21:  J496.A - Gulf (95%) / MW, LA; GULF_09
22:  J496.C - Gulf (95%) / MW, LA; GULF_09
23:  J497.A - Gulf (95%) / MW, LA; GULF_09
24:  J497.B - Gulf (95%) / MW, LA; GULF_09
25:  J500.A - Gulf (97%) / ATL, MS; GULF_03
26:  J502.A - Gulf (99%) / ATL, MS; GULF_03
27:  J502.C - Gulf (97%) / ATL, MS; GULF_03
28:  J589_A - ATL (98%) / Gulf, VA; ATL_10; in subclade with P. amarum from NJ
29:  J589.B - ATL (98%) / Gulf, VA; ATL_10;
30:  J589.C - ATL (98%) / Gulf, VA; ATL_10
31:  J591.A - ATL (98%) / Gulf, VA; ATL_07; in same subclade with P. amarum as
#  J589
32:  J591.B - ATL (98%) / Gulf, VA;
33:  J593.B - ATL, NC; ATL_05
34:  J595.B - ATL, NC; ATL_05
35:  J595.C - ATL, NC; ATL_05
36:  J598.B - ATL, ML; ATL_08; this is sister to High Tide clutivar
37:  J613.A - ATL (71%) / Gulf, NC; ATL_07; is P. amarum !?!
38:  J613.C - ATL (68%) / Gulf, NC; ATL_07; is P. amarum !?!
39:  J621.B - ATL, ML; ATL_08; another sample from this accession is in 
#       a cluster with P. amarum from NJ!
40:  J621.C - ATL (94%)/ Gulf; ML; in cluster with P. amarum from NJ!
41:  J676.A - ATL, VA; ATL_08
42:  J676.C - ATL, VA; ATL_08
43:  J677.B - ATL, NC; ATL_05
44:  J682.B - ATL, NJ; ATL_10

######

weird_gulf_inds <- intersect(res_weird_inds, grep('GULF_', res_tab$subgrp_v2))
weird_gulf_ord <- order(res_tab[weird_gulf_inds, samp_name])
res_tab[weird_gulf_inds[weird_gulf_ord], list(samp_name, ploidy)]
as.data.table(sort(res_tab[weird_gulf_inds, samp_name])

tot_8X_gulf_inds <- intersect(grep('GULF_', res_tab$subgrp_v2), 
  which(res_tab$ploidy == '8X'))

setdiff(tot_8X_gulf_inds, weird_gulf_inds)
# 7 of 17 not part of upland cp clade

    samp_name ploidy
 1:    J177.A     4X - Gulf(75%)/ATL admix, FL; GULF_12
 2:   J179.L2     4X - Gulf(72%)/ATL admix; FL; GULF_12
 3:   J179.L3     4X - Gulf(85%)/ATL admix; FL; GULF_12
 4:    J188.A     4X - Gulf(93%)/ATL admix; FL; GULF_08
 5:   J188.L1     8X - Gulf(88%)/ATL admix; FL; GULF_09 but same asscession as
# other J188's
 6:   J188.L2     4X - Gulf(93%)/ATL admix; FL; GULF_08
 7:   J188.L3     4X - Gulf(95%)/ATL admix; FL; GULF_08
 8:    J189.A     4X - Gulf(75%)/ATL admix; FL; GULF_09
 9:    J209.A     4X  - Gulf, TX; GULF_06 - NOTE: GULF_06 includes 8X samples 
## with some MW ancestry; this clade also includes a weird 6X sample from 
## FL that GRIN says is P. amarum...
10:    J209.B     4X - Gulf, TX; GULF_06
11:    J213_A     8X - Gulf(97%)/MW admix; TX; GULF_06
12:    J221_A     8X - Gulf(99%)/MW admix; TX: GULF_06
13:    J228.A     4X - Gulf, TX; GULF_06 
14:    J236_A     8X - Gulf(97%)/MW admix; TX; GULF_06
15:    J240.A     4X - Gulf, TX; GULF_06
16:    J294.A     4X - Gulf, TX: GULF_06
17:    J302.A     8X - Gulf(97%)/MW admix; TX; GULF_06
18:    J303.A     4X - Gulf, TX: GULF_06
19:    J304_A     8X - GULF(99%)/MW admix; TX; GULF_06
20:    J305.A     4X - Gulf, TX: GULF_06
21:    J306.A     4X
22:    J306.B     4X
23:    J306.C     4X
24:    J312.A     4X
25:    J313.A     4X
26:    J328_A     8X
27:    J337.A     4X
28:    J461.A     4X
29:    J461.B     4X
30:    J461.C     4X
31:    J462_B     8X
32:    J462_C     8X
33:    J465.A     4X
34:    J465.B     4X
35:    J465.C     4X
36:    J466_D     4X
37:    J466.A     4X
38:    J466.B     4X
39:    J466.C     4X
40:    J477.A     4X
41:    J477.C     4X
42:    J482.A     4X
43:    J482.B     4X
44:    J483.C     4X
45:    J495_C     8X
46:    J496.A     4X
47:    J496.C     4X
48:    J497.A     4X
49:    J497.B     4X
50:    J499_B     4X
51:    J499.A     4X
52:    J499.C     4X
53:    J500.A     4X
54:    J501_B     4X
55:    J501.C     4X
56:    J502.A     4X
57:    J502.C     4X
58:    J514.A     4X
59:    J535.A     4X




