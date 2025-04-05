from treeswift import read_tree_newick
import pandas as pd

## suboscines dataset

wastral_castles = read_tree_newick("../biological/suboscines-harvey/treepl_castles_T400F.wastral.rooted.ingroup.tre")
castles_astral = read_tree_newick("../biological/suboscines-harvey/treepl_castles_T400F.astral.rooted.ingroup.tre")
caml_astral = read_tree_newick("../biological/suboscines-harvey/treepl_examl_T400F.astral.rooted.ingroup.tre")
caml_caml = read_tree_newick("../biological/suboscines-harvey/treepl_concat_T400F.examl.rooted.ingroup.tre")
caml_castles = read_tree_newick("../biological/suboscines-harvey/treepl_castles_T400F.examl.rooted.ingroup.tre")

print(wastral_castles.treeness())
print(castles_astral.treeness(), caml_astral.treeness())
print(caml_castles.treeness(), caml_caml.treeness())

## stiller dataset
'''
castles_mdcat= read_tree_newick("../biological/avian-stiller/mdcat_median_40Kl_astral4_stiller_nolabel.rooted.tre")
caml_mdcat = read_tree_newick("../biological/avian-stiller/mdcat_median_40Kl_caml_stiller.rooted.tre")
castles_treepl= read_tree_newick("../biological/avian-stiller/treepl_median_castlespro_stiller.rooted.tre")
caml_treepl = read_tree_newick("../biological/avian-stiller/treepl_median_astral_63K_concat_bl.rooted.tre")
castles_wlogdate= read_tree_newick("../biological/avian-stiller/wlogdate_median_astral4_stiller.rooted.no_label.tre")
caml_wlogdate = read_tree_newick("../biological/avian-stiller/wlogdate_median_astral_63K_concat_bl.rooted.no_label.tre")

print(castles_mdcat.treeness(), caml_mdcat.treeness())
print(castles_treepl.treeness(), caml_treepl.treeness())
print(castles_wlogdate.treeness(), caml_wlogdate.treeness())
'''
