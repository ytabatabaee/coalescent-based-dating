import argparse
import matplotlib.pyplot as plt
import dendropy
import math
import pandas as pd
import seaborn as sns
import numpy as np


def plot_x_y_line(ax):
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]),
            np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, '--', alpha=0.75, zorder=0, color='black')
    ax.set_aspect('equal')


def plot_correlations(df, name1, name2):
    plt.cla()
    sns.lmplot(data=df, x="log10(l1)", y="log10(l2)", palette='Dark2', hue="Node Type",
               scatter_kws={'alpha': 0.8, 'linewidth': 0}, order=2, ci=0) #
    plt.grid(linestyle='--', linewidth=0.5)
    ax = plt.gca()
    ax.set_xlabel('CASTLES-Pro')
    ax.set_ylabel('CAML')
    ax.text(-1.3, -1.5, 'y=x')
    plot_x_y_line(ax)
    plt.savefig(name1+'_'+name2+'_correlations.pdf', bbox_inches='tight')


def calculate_node_age(t):
    for node in t.preorder_node_iter():
        if node.parent_node:
            node.age = node.edge.length + node.parent_node.age # time
        else: # root node
            node.age = 0

    tree_height = t.find_node_with_taxon_label('Cariama Y101002').age # use a taxa name from the tree
    print(tree_height)

    # backward time (present leafs 0, root oldest)
    for node in t.preorder_node_iter():
        node_type = 'terminal' if node.is_leaf() else 'internal'
        node.age = round(abs(tree_height - node.age), 3)
        node.edge.length = node.age

    return t


def compare_node_age(tree_path1, tree_path2):
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(path=tree_path1, schema='newick', taxon_namespace=tns, rooting='force-rooted')
    t2 = dendropy.Tree.get(path=tree_path2, schema='newick', taxon_namespace=tns, rooting='force-rooted')

    t1 = calculate_node_age(t1)
    t2 = calculate_node_age(t2)

    df_branches = pd.DataFrame(columns=['Taxon', "Node Type", "l1", "l2", 'log10(l1)', 'log10(l2)', 'l1/l2'])

    length_diffs = dendropy.calculate.treecompare._get_length_diffs(t1, t2)
    idx = 0

    for node in t1.postorder_node_iter():
        (l1, l2) = length_diffs[idx]
        node_label = node.taxon.label if node.is_leaf() else '  '
        node_type = 'terminal' if node.is_leaf() else 'internal'
        if node_type == 'internal':
            df_branches.loc[len(df_branches.index)] = [node_label, node_type, l1, l2,
                                                       math.log10(l1) if l1 > 0 else np.nan,
                                                       math.log10(l2) if l2 > 0 else np.nan,
                                                       l1/l2 if l2 > 0 else np.nan]
        idx += 1

    print(df_branches[['Taxon', "Node Type", "l1", "l2", "l1/l2"]].to_string())
    name1 = args.tree1.split('/')[-1].split('.')[0].upper()
    name2 = args.tree2.split('/')[-1].split('.')[0].upper()
    df_branches[['Taxon', "Node Type", "l1", "l2"]].to_csv(name1 + '_' + name2 + '_node_ages.csv')
    df_branches = df_branches.dropna()
    print(len(df_branches))

    if args.plot:
        plot_correlations(df_branches, name1, name2)

    df_branches['l1'] = df_branches['l1'].apply(lambda x: x if x > 0 else 1e-6)
    df_branches['l2'] = df_branches['l2'].apply(lambda x: x if x > 0 else 1e-6)
    print('\nBias:', np.mean(df_branches['l1'] - df_branches['l2']))
    print('Mean absolute error:', np.mean(np.abs(df_branches['l1'] - df_branches['l2'])))
    print('Root mean square error (RMSE):', np.sqrt(np.mean((df_branches['l1'] - df_branches['l2'])**2)))
    print('Mean logarithmic error:', np.mean(np.abs(np.log10(df_branches['l1']) - np.log10(df_branches['l2']))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compare node ages of two trees")
    parser.add_argument("-t1", "--tree1", type=str,  required=True,
                        help="tree file with branch lengths in newick format")
    parser.add_argument("-t2", "--tree2", type=str, required=True,
                        help="tree file with branch lengths in newick format")
    parser.add_argument("-p", "--plot", default=False, required=False,
                        action='store_true',
                        help="plot correlations between node ages of input trees")
    args = parser.parse_args()
    compare_node_age(args.tree1, args.tree2)
