"""
This script uses the jaccard indices of the jaccard.py script. 
The jaccard distances are correlated with the phylogenetic distances.
A correlation plot with a regression line is produced as a panel plot.
They latter need to be provided in the subfolder "/phylogenetic trees/"
@Author: Maximilian Sprang, Github: Muedi
"""
# %%
from itertools import combinations
from matplotlib import colors
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, kendalltau, pearsonr
import pandas as pd
import numpy as np
import seaborn as sns

# %%from matplotlib import colors
names = {
    "ec":"E. cloacae",
    "ef":"E. faecium",
    "ab":"A. baumannii",
    "pa":"P. aeruginosa",
    "kp":"K. pneumoniae",
    "sa":"S. aureus"
}
species = ["ef", "sa", "kp", "ab", "pa", "ec"]
correls = pd.DataFrame(
    index=["pearson", "spearman", "kendall"], columns=species)
# %%

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 12))
for s, ax in zip(species, axs.ravel()):
    print(s)
    name = names[s]
    file_meta = "data paper/metadata_{}_def_brief.tsv".format(s)
    meta = pd.read_csv(file_meta, compression='infer', sep="\t", engine='c')
    meta = meta.set_index("ID")
    # dorp ambigous
    meta = meta.loc[meta.Ambiguous.ne("ambiguous").values]
    # data = data.loc[:, meta.index]
    # remove MLST group -
    meta = meta.loc[meta.MLST.ne("-").values]
    ab_jaccard_distance = pd.read_csv("output/_{}_jacc_dissimilaity_MLST.csv".format(s), index_col=0)
    ab_jaccard_distance = ab_jaccard_distance.add_prefix("{}_".format(s)).T
    ab_jaccard_distance = ab_jaccard_distance.add_prefix("{}_".format(s)).T
    ab_phylogenic_distance = pd.read_csv(
        "phylogenetic trees/{}_dist_matrix.tsv".format(s), sep="\t", index_col=0)

    common_mlst = [mlst for mlst in ab_jaccard_distance.columns.tolist()
                if mlst in ab_phylogenic_distance.columns.tolist()]
    # get the same index and overlapping mlsts
    ab_phylogenic_distance = ab_phylogenic_distance.loc[common_mlst,
                                                        common_mlst]
    ab_jaccard_distance = ab_jaccard_distance.loc[common_mlst,
                                                common_mlst]
    # deframe, keep indices
    ab_phylogenic_distance_flat = pd.Series(squareform(ab_phylogenic_distance),
                                            index=pd.MultiIndex.from_tuples(tuple(combinations(ab_phylogenic_distance.index, 2)))
                                            )
    ab_jaccard_distance_flat = pd.Series(squareform(ab_jaccard_distance),
                                    index=pd.MultiIndex.from_tuples(tuple(combinations(ab_jaccard_distance.index, 2)))
                                    )
    # compute correlations
    _pearsonr = pearsonr(ab_jaccard_distance_flat,
                         ab_phylogenic_distance_flat)
    _spearmanr = spearmanr(ab_jaccard_distance_flat,
                           ab_phylogenic_distance_flat)
    _kendalltau = kendalltau(ab_jaccard_distance_flat,
                             ab_phylogenic_distance_flat)
    correls[s] = [_pearsonr[0], _spearmanr[0], _kendalltau[0]]

    # the block below does some cosmetics on the crispr names
    # If multiple are present they are split into a list of strings
    # Ambigous strains need to be dropped as in the jaccard script.
    # get df with true false for crispr type with types as index and mlst as cols (test)
    test = meta.drop("Ambiguous", axis=1)
    test["crispr_clean"] = test.CRISPRtype.str.split(",")
    test = test.explode("crispr_clean")
    test["idx"] = test.groupby('MLST')["crispr_clean"].cumcount()
    test = test.pivot_table(index="crispr_clean", columns="MLST", values="idx")
    test = test.notnull()
    # build numeric map
    mlst_map = pd.Series(["No_crispr"] * len(set(meta.MLST)), index=set(meta.MLST))
    mlst_map.loc[test.columns] = [
        test[mlst].loc[test[mlst].eq(True)].index.values for mlst in test.columns]
    mlst_map = mlst_map.astype(str)
    unique_combis = list(set(mlst_map))
    unique_combis.remove("No_crispr")
    unique_combis = ["No_crispr"] + unique_combis
    unique_combis = {unique_combis[i]:i for i in range(len(unique_combis))}
    # build colormap
    colors_hex = ["#000000","#FFBF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C","#92896B"]
    colors_rgb = [colors.to_rgb(color_hex) for color_hex in colors_hex]


    # plot preprocessing
    df_ml = pd.DataFrame({"mlsts": ab_jaccard_distance_flat.index.to_list(),
                          "Jaccard distance": ab_jaccard_distance_flat,
                          "Phylogenetic distance": ab_phylogenic_distance_flat,
                          })
    df_ml["Crispr-Status"] = ["Both Crispr-Cas" if (mlst_map[x[0][3:]] != "No_crispr" and mlst_map[x[1][3:]] != "No_crispr") else "One Crispr-Cas" if (
        mlst_map[x[0][3:]] != "No_crispr" or mlst_map[x[1][3:]] != "No_crispr") else "None Crispr-Cas" for x in df_ml.mlsts.to_list()]
    # TODO: use mlsts list of groups to index mlst_map and get info if one or both of the mlsts contains a crispr group
    # plot
    sns_palette = {"One Crispr-Cas": colors_hex[1],
                   "None Crispr-Cas": colors_hex[2],
                   "Both Crispr-Cas": colors_hex[3]}
    # sns_rgb = {"One Crispr": colors_rgb[1], 
    #                "None Crispr": colors_rgb[2],
    #                "Both Crispr": colors_rgb[3]}
    #fig, ax = plt.subplots(figsize=(5, 5))
    sns.scatterplot(data=df_ml,
                    x="Phylogenetic distance",
                    y="Jaccard distance",
                    hue="Crispr-Status",
                    markers=['o'],
                    palette=sns_palette,
                    alpha=0.7,
                    ax=ax, 
                    legend=None
                    )
    sns.regplot(data=df_ml,
                x="Phylogenetic distance",
                y="Jaccard distance",
                marker='None',
                #order=3,
                ax=ax
                )
    ax.set_ylabel("Gene profile distance")
    ax.set_xlabel("Phylogenetic distance")
    ax.set_title("$\it{}$ (corr = {})".format(name, np.round(_pearsonr[0], 2)))
    # fig.savefig("output/{}_jacc-vs-phylo-lmplot-linear.svg".format(s))

handles = [Patch(color=sns_palette[key], label=key) for key in sns_palette.keys()]
axs[1, 2].legend(handles=handles, loc="lower right")
fig.savefig("output/panel_jacc-vs-phylo-lmplot-linear.png", dpi=300)
fig.savefig("output/panel_jacc-vs-phylo-lmplot-linear.tiff", dpi=300)
fig.savefig("output/panel_jacc-vs-phylo-lmplot-linear.pdf")

print("finished")