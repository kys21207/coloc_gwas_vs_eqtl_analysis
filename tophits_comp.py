import gwaslab as gl

dta1 = gl.Sumstats("ALL-GWAS_ALL-GWAS_ENSG00000139874.regenie_block_1163.txt", fmt="regenie")
dta2 = gl.Sumstats("coloc_inputs/batch1/gwasA_block_1163.txt", fmt="regenie")
#dta2= gl.Sumstats("coloc_inputs/batch1/ALL-GWAS_ALL-GWAS_ENSG00000139874.regenie_block_1163.txt",
#             snpid="ID",
#             chrom="CHROM",
#             pos="GENPOS",
#             ref="ALLELE0",
#             alt="ALLELE1",
#             beta="BETA",
#             se="SE",
#             mlog10p="LOG10P",
#             sep=" ")


dta1.basic_check()
dta2.basic_check()

dta1.fill_data(to_fill=["P"], extreme=True)
dta2.fill_data(to_fill=["P"], extreme=True)


# auto extract lead SNPs and compare the effect sizes
a = gl.compare_effect(path1= dta1,
                      path2= dta2,
                      sig_level=1,
                      legend_title=r'$ P < 1 x 10^{-1}$ in:',
                      save="summary_plots/ms_all_comers_gwas_eqtl_compare_plot.png"
                      )

#a = gl.plot_miami2(path1= dta1,
#                   path2= dta2,
#                   build = "38",
#                   skip=2,
#                   cut1=20,
#                   cut2=15,
#                   id1="SNPID",
#                   id2="SNPID",
#                   anno1="GENENAME",
#                   anno2="GENENAME",
##                   additional_line1=[1e-14],
##                  anno_set1=["rs3798519"],
##                  pinpoint1=[["rs3798519","rs35560038"],["rs7933262","rs8098510"]],
##                  pinpoint_color1=["purple","black"],
##                  highlight1=["rs8098510"],
##                  highlight2=[["rs8098510","rs3798519"], ["rs1491850"]],
#                   jagged=True,
#                   verbose1=False,
#                   verbose2=False,
#                   save="/home/jovyan/session_data/summary_plots/ms_regenie_miami_plot.png"
#)

#a=dta1.plot_mqq(mode="m",
#              region=(16,10844539,12182186),      
#              build = "38",
#        stratified = True,
#        sig_level = 1e-6,
#        suggestive_sig_line = True,
#        suggestive_sig_line_color = "red",
#        anno = "GENENAME",
#        anno=True,
#                    save="summary_plots/ms_all_comers_manhattan_plot.png"

#                    anno_set=["7:156793450_G_GA","7:157038803_A_G"]
#                   )

