from treeswift import read_tree_newick
import pandas as pd

## suboscines dataset

wastral_castles = read_tree_newick("../biological/suboscines-harvey/treepl_castles_T400F.wastral.rooted.ingroup.tre")
castles_astral = read_tree_newick("../biological/suboscines-harvey/treepl_castles_T400F.astral.rooted.ingroup.tre")
caml_astral = read_tree_newick("../biological/suboscines-harvey/treepl_examl_T400F.astral.rooted.ingroup.tre")
caml_caml = read_tree_newick("../biological/suboscines-harvey/treepl_concat_T400F.examl.rooted.ingroup.tre")
caml_castles = read_tree_newick("../biological/suboscines-harvey/treepl_castles_T400F.examl.rooted.ingroup.tre")

# old present day : 64.7143

max_pd=40.4552 #

dict1 = caml_castles.lineages_through_time(color="green",show_plot=False,present_day=max_pd)
dict2 = caml_caml.lineages_through_time(color="lightgreen",show_plot=False,present_day=max_pd)
dict3 = castles_astral.lineages_through_time(color="blue",show_plot=False,present_day=max_pd)
dict4 = caml_astral.lineages_through_time(color="lightblue",show_plot=False,present_day=max_pd)
dict5 = wastral_castles.lineages_through_time(color="red",present_day=max_pd,export_filename='ltt.pdf')


df = pd.DataFrame(columns=['method', 'topology', 'time', 'lineages'])

for key, value in dict1.items():
    df.loc[len(df)] = ['TreePL+CASTLES-Pro', 'CAML', key, value]
for key, value in dict2.items():
    df.loc[len(df)] = ['TreePL+Concat(RAxML)', 'CAML', key, value]

for key, value in dict3.items():
    df.loc[len(df)] = ['TreePL+CASTLES-Pro', 'ASTRAL', key, value]
for key, value in dict4.items():
    df.loc[len(df)] = ['TreePL+Concat(RAxML)', 'ASTRAL', key, value]

for key, value in dict5.items():
    df.loc[len(df)] = ['TreePL+CASTLES-Pro', 'wASTRAL', key, value]

df.to_csv('suboscines_ltt_ingroup.csv')

#df.to_csv('caml_ltt_ingroup.csv')
'''
print(wastral_castles.treeness())
print(castles_astral.treeness(), caml_astral.treeness())
print(caml_castles.treeness(), caml_caml.treeness())
'''

## stiller dataset
'''
castles_mdcat= read_tree_newick("../biological/avian-stiller/mdcat_median_40Kl_astral4_stiller_nolabel.rooted.tre")
caml_mdcat = read_tree_newick("../biological/avian-stiller/mdcat_median_40Kl_caml_stiller.rooted.tre")

dict1 = castles_mdcat.lineages_through_time(color="green",show_plot=False,present_day=127.28139999999999)
dict2 = caml_mdcat.lineages_through_time(color="lightgreen",export_filename='ltt_stiller.pdf',present_day=127.28139999999999)

df = pd.DataFrame(columns=['method', 'time', 'lineages'])

for key, value in dict1.items():
    df.loc[len(df)] = ['MD-Cat+CASTLES-Pro', key, value]
for key, value in dict2.items():
    df.loc[len(df)] = ['MD-Cat+Concat(RAxML)', key, value]

df.to_csv('ltt_stiller.csv')
'''

'''
castles_treepl= read_tree_newick("../biological/avian-stiller/treepl_minmax_castlespro_stiller.rooted.tre")
caml_treepl = read_tree_newick("../biological/avian-stiller/treepl_minmax_astral_63K_concat_bl.rooted.tre")
castles_wlogdate= read_tree_newick("../biological/avian-stiller/wlogdate_median_astral4_stiller.rooted.no_label.tre")
caml_wlogdate = read_tree_newick("../biological/avian-stiller/wlogdate_median_astral_63K_concat_bl.rooted.no_label.tre")

print(castles_mdcat.treeness(), caml_mdcat.treeness())
print(castles_treepl.treeness(), caml_treepl.treeness())
print(castles_wlogdate.treeness(), caml_wlogdate.treeness())
'''
