import argparse
import matplotlib.pyplot as plt
import dendropy
import math
import pandas as pd
import seaborn as sns
import numpy as np
import re


def extract_node_metadata(node):
    metadata = {
        't': None,
        't_lower': None,
        't_upper': None,
        'mu': None,
        'mu_lower': None,
        'mu_upper': None,
        'bl': None,
        'bl_lower': None,
        'bl_upper': None
    }

    # ---- Extract ALL metadata from node.annotations ----
    if node.annotations:
        for key in metadata.keys():
            value = node.annotations.get_value(key)
            if value is not None:
                try:
                    metadata[key] = float(value)
                except (TypeError, ValueError):
                    pass  # leave as None if conversion fails

    # ---- Branch length ----
    if node.edge and node.edge.length is not None:
        try:
            metadata['bl'] = float(node.edge.length)
        except (TypeError, ValueError):
            pass

    return metadata


def is_within(point, lower, upper):
    if lower is None or upper is None or point is None:
        return None
    return lower <= point <= upper


def compute_CI_containment(df):

    coal = df[df["Method"] == "MD-Cat+CoalBL"]
    con  = df[df["Method"] == "MD-Cat+ConBL"]

    merged = coal.merge(
        con,
        on=["Calibrations", 'Condition','Replicate','clade_id'],
        suffixes=('_coal','_con')
    )

    merged['coal_in_con_CI'] = (
        (merged['t_coal'] >= merged['t_lower_con']) &
        (merged['t_coal'] <= merged['t_upper_con'])
    )

    return merged


def compute_CI_zscores(df, param):
    """
    Compute Z-scores of CoalBL point estimates relative to ConBL CIs,
    safely handling missing or zero-width CIs.
    """
    coal = df[df["Method"] == "MD-Cat+CoalBL"]
    con  = df[df["Method"] == "MD-Cat+ConBL"]

    merged = coal.merge(
        con,
        on=["Calibrations", "Condition", "Replicate", "clade_id"],
        suffixes=("_coal", "_con")
    )

    # Midpoint and half-width of ConBL CI
    midpoint = (merged[param+'_lower_con'] + merged[param+'_upper_con']) / 2
    half_width = (merged[param+'_upper_con'] - merged[param+'_lower_con']) / 2

    # Safe z-score: only compute if half_width > 0
    merged['z_score'] = merged.apply(
        lambda row: (row[param+'_coal'] - (row[param+'_lower_con'] + row[param+'_upper_con'])/2) / ((row[param+'_upper_con'] - row[param+'_lower_con'])/2)
        if pd.notna(row[param+'_coal']) and pd.notna(row[param+'_lower_con']) and pd.notna(row[param+'_upper_con']) and (row[param+'_upper_con'] > row[param+'_lower_con'])
        else np.nan,
        axis=1
    )

    # Flag overlap (inside CI)
    merged['coal_in_con_CI'] = (merged[param+'_coal'] >= merged[param+'_lower_con']) & \
                               (merged[param+'_coal'] <= merged[param+'_upper_con'])

    # Overlap ratio
    valid = merged.dropna(subset=['coal_in_con_CI'])
    overlap_ratio = valid['coal_in_con_CI'].sum() / len(valid)

    # Safe z-score stats
    z_valid = merged['z_score'].dropna()
    mean_z = z_valid.mean()
    std_z = z_valid.std()
    n_1 = (z_valid.abs() > 1).sum()
    n_2 = (z_valid.abs() > 2).sum()

    print(len(coal), len(con), len(merged))

    print(param, param+'_lower', param+'_upper')
    print(f"{valid['coal_in_con_CI'].sum()}/{len(valid)} nodes ({overlap_ratio:.2%}) inside ConBL CI | "
          f"Z-score: mean={mean_z:.3f}, std={std_z:.3f} | "
          f"{n_1} nodes >1 CI width away, {n_2} nodes >2 CI widths away")
    z = merged['z_score'].dropna()
    print("Median Z:", z.median())
    print("Mean Z:", z.mean())
    print("95th percentile:", z.quantile(0.95))
    print("99th percentile:", z.quantile(0.99))
    print("Max Z:", z.max())
    width = merged[param+'_upper_con'] - merged[param+'_lower_con']
    print((width < 1e-10).sum())

    return merged, overlap_ratio


def summarize_overlap(merged):
    valid = merged.dropna(subset=['coal_in_con_CI'])

    total = len(valid)
    inside = valid['coal_in_con_CI'].sum()

    if total == 0:
        print("No valid comparisons.")
        return None

    print(f"{inside}/{total} nodes ({inside/total:.2%}) "
          "have CoalBL inside ConBL CI")

    return inside / total


def summarize_zscores(merged):
    # Drop nodes with missing z-score
    valid = merged.dropna(subset=['z_score'])

    total = len(valid)
    if total == 0:
        print("No valid nodes for z-score calculation.")
        return None

    mean_z = valid['z_score'].mean()
    std_z = valid['z_score'].std()

    # Count nodes inside ConBL CI
    inside = valid['coal_in_con_CI'].sum()

    # Count nodes exceeding 1 or 2 CI widths
    over1 = (valid['z_score'].abs() > 1).sum()
    over2 = (valid['z_score'].abs() > 2).sum()

    print(
        f"{inside}/{total} nodes ({inside / total:.2%}) inside ConBL CI | "
        f"Z-score: mean={mean_z:.3f}, std={std_z:.3f} | "
        f"{over1} nodes >1 CI width away, {over2} nodes >2 CI widths away"
    )


if __name__ == "__main__":
    model_conditions = [
        "outgroup.0.species.0.15.genes.0.15",
        "outgroup.1.species.1.5.genes.1.5",
        "outgroup.0.species.5.genes.5",
        "outgroup.1.species.5.genes.5",
        "outgroup.0.species.1.5.genes.1.5",
        "outgroup.1.species.0.15.genes.0.15"
    ]
    rows = []

    calibs = ['', 'n3_', 'n5_']
    calib_num = [0, 3, 5]

    ru_txts = ['']  # add 'root_unfixed_' if needed
    ru_states = [True]

    methods = ['MD-Cat+CoalBL', 'MD-Cat+ConBL']

    dataset_path = '/u/syt3/scratch/dating-data/MVroot/'
    gt_type = 'estimatedgenetre.gtr'
    footer = '_s_tree.trees.rooted.labeled'

    df = pd.DataFrame(columns=[
        "Calibrations",
        "Condition",
        "root_fixed",
        "Method",
        "Replicate",
        'gtee',
        'ad',
        "clade_id",
        "node_type",
        "t",
        "t_lower",
        "t_upper",
        "mu",
        "mu_lower",
        "mu_upper",
        "bl",
        "bl_lower",
        "bl_upper"
    ])

    for c in range(len(calibs)):
        calib = calibs[c]
        calib_value = calib_num[c]

        for u in range(len(ru_txts)):

            if calib == '' and not ru_states[u]:
                continue

            for condition in model_conditions:

                print(f"Processing: calib={calib_value}, condition={condition}")

                for r in range(1, 101):
                    replicate = str(r).zfill(3)

                    tns = dendropy.TaxonNamespace()

                    ad_path = dataset_path + condition + '/' + replicate + '/ad.txt'
                    gtee_path = dataset_path + condition + '/' + replicate + '/gtee_gtr.txt'
                    try:
                        with open(ad_path) as f:
                            ad = float(f.read())
                        with open(gtee_path) as g:
                            gtee = float(g.read())
                    except:
                        print('no ad or gtee file')

                    for method in methods:

                        try:
                            # Build tree filename
                            if method == "MD-Cat+ConBL":
                                tree_prefix = f"mdcat_CI_{calib}{ru_txts[u]}RAxML_result.concat_align"
                                tree_path = (
                                    dataset_path +
                                    condition + "/" +
                                    replicate + "/" +
                                    tree_prefix +
                                    footer
                                )
                            else:
                                tree_prefix = f"mdcat_CI_{calib}{ru_txts[u]}castlespro_"
                                tree_path = (
                                    dataset_path +
                                    condition + "/" +
                                    replicate + "/" +
                                    tree_prefix +
                                    gt_type +
                                    footer
                                )

                            # Load tree

                            with open(tree_path) as f:
                                newick = f.read()

                            # Convert CI to proper metadata comment
                            newick = re.sub(
                                r':\s*([0-9.eE+-]+)\s*\[\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\]',
                                r':\1[&bl_lower=\2,bl_upper=\3]',
                                newick
                            )

                            newick = re.sub(
                                r"'(\w+)\[([^\]]+)\]'",
                                r"\1[&\2]",
                                newick
                            )

                            tree = dendropy.Tree.get(
                                data=newick,
                                schema="newick",
                                taxon_namespace=tns,
                                preserve_underscores=True,
                                extract_comment_metadata=True
                            )

                            tree.encode_bipartitions()

                            # Extract node metadata
                            for node in tree.postorder_node_iter():

                                if node.parent_node is None:
                                    continue  # skip root (no meaningful bipartition)

                                node_type = 'terminal' if node.is_leaf() else 'internal'

                                # Unique clade identifier (bitstring of descendant taxa)
                                clade_id = node.edge.bipartition.leafset_as_bitstring()
                                meta = extract_node_metadata(node)

                                rows.append([
                                    calib_value,
                                    condition,
                                    ru_states[u],
                                    method,
                                    replicate,
                                    gtee,
                                    ad,
                                    clade_id,
                                    node_type,
                                    meta['t'],
                                    meta['t_lower'],
                                    meta['t_upper'],
                                    meta['mu'],
                                    meta['mu_lower'],
                                    meta['mu_upper'],
                                    meta['bl'],
                                    meta['bl_lower'],
                                    meta['bl_upper']
                                ])

                        except Exception as e:
                            print(f"Skipping {condition} rep {replicate} ({method})")
                            print(e)

    print("Extraction complete.")

    df = pd.DataFrame(rows, columns=[
        "Calibrations",
        "Condition",
        "root_fixed",
        "Method",
        "Replicate",
        'gtee',
        'ad',
        "clade_id",
        "node_type",
        "t",
        "t_lower",
        "t_upper",
        "mu",
        "mu_lower",
        "mu_upper",
        "bl",
        "bl_lower",
        "bl_upper"
    ])

    df.to_csv("mvroot_MDcat_CI.csv", index=False)
    print("Computing containment...")

    #merged = compute_CI_containment(df)
    # summarize_overlap(merged)

    # Save results
    #merged.to_csv("coal_vs_con_CI_containment.csv", index=False)

    merged_with_z, overlap_ratio = compute_CI_zscores(df, param='bl')
    merged_with_z, overlap_ratio = compute_CI_zscores(df, param='t')
    merged_with_z, overlap_ratio = compute_CI_zscores(df, param='mu')
    #summarize_zscores(merged_with_z)
    merged_with_z.to_csv("coal_vs_con_CI_zscores.csv", index=False)
