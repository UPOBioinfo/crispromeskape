"""
This script computes the jaccard distance between the MLST groups of each species.
The MLST groups are condensed from all strains that are marked as the respectiuve group,
for each gene the majority vote of all corresponding strands appoints the genes status.
o fthese distances and the phylogenetic distances of the strains a hierarchical clustering
is performed and plotted as a dendrogram.
They latter need to be provided in the subfolder "/phylogenetic trees/"
@Author: Maximilian Sprang, Github: Muedi
"""# %%
# jaccard index/distance
# from itertools import combinations
from scipy.spatial.distance import squareform
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import scipy.stats as st
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import colors
from scipy.stats import spearmanr, kendalltau, pearsonr
# import seaborn as sns

# cores for jaccard index 
n_cores = 6

def drop_prefix(self, prefix):
    self.columns = self.columns.str.lstrip(prefix)
    return self

pd.core.frame.DataFrame.drop_prefix = drop_prefix

filenames = ["data paper/" + f for f in listdir("data paper/")
             if isfile(join("data paper/", f))]

                 # build colormap
colors_hex = ["#000000","#FFBF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C","#92896B"]
colors_rgb = [colors.to_rgb(color_hex) for color_hex in colors_hex]

unique_combis = ["No_crispr", "I + II", "I + III", "II + III", "I + IV", "IV", "III", "II", "I"]
unique_combis = {unique_combis[i]:i for i in range(len(unique_combis))}
my_palette = colors.ListedColormap(colors_rgb[1:len(unique_combis)+1])

# %%
names = {
    "_ec":"E. cloacae",
    "_ef":"E. faecium",
    "_ab":"A. baumannii",
    "_pa":"P. aeruginosa",
    "_kp":"K. pneumoniae",
    "_sa":"S. aureus"
}
species = ["_ef", "_sa", "_kp", "_ab", "_pa", "_ec"]
# list_keys = list(load_dict.keys())
correls = pd.DataFrame(index=["pearson", "spearman", "kendall"], columns=species)
fig, axs = plt.subplots(nrows=2, ncols=6, figsize=(20, 20))
for idx, s in enumerate(species):
    name = names[s]
    file_data = [
        file_data for file_data in filenames if s in file_data and "gnmatrix" in file_data][0]
    data = pd.read_csv(file_data, compression='infer', sep="\t", memory_map=True, engine='c')
    data = data.set_index("#genes/strains")
    #data = data.T # long format
    file_meta = [
        file_meta for file_meta in filenames if s in file_meta and "metadata" in file_meta][0]
    meta = pd.read_csv(file_meta, compression='infer', sep="\t", engine='c')
    meta = meta.set_index("ID")
    # dorp ambigous
    meta = meta.loc[meta.Ambiguous.ne("ambiguous").values]
    # data = data.loc[:, meta.index]
    # remove MLST group -
    meta = meta.loc[meta.MLST.ne("-").values]
    data = data.loc[:, meta.index]
    
    if not isfile("output/{}_jacc_dissimilaity_MLST.csv".format(s)):
        condensed_data = pd.DataFrame(index=set(meta.MLST), columns=data.index)
        # GET majority vote 
        for mlst in set(meta.MLST):
            # majority vote
            c0 = data.loc[:, meta.MLST.eq(mlst)].values
            condensed_data.loc[mlst] = st.mode(c0, axis=1)[0].flatten()


        #pairwise distance
        data_MLST = pairwise_distances(
            condensed_data.values, metric='jaccard', n_jobs=4)

        # for combi in combis_mlst:
        #     jacc_dis = distance.jaccard(majority_c0, majority_c1)
        #     data_MLST[combi[0]][combi[1]] = jacc_dis
        #     data_MLST[combi[1]][combi[0]] = jacc_dis
        #     break
        data_MLST_nonsquare = squareform(data_MLST)
        data_MLST = pd.DataFrame(
            data_MLST, index=condensed_data.index, columns=condensed_data.index)
        data_MLST.to_csv("output/{}_jacc_dissimilaity_MLST.csv".format(s))
        print("{} finished".format(s))
    else:
        data_MLST = pd.read_csv(
            "output/{}_jacc_dissimilaity_MLST.csv".format(s), index_col=0)
        data_MLST_nonsquare = squareform(data_MLST)

    # dendro
    # get mapping of mlst and crispr, if True the mlst group contains at least one Crispr strain
    # mlst_map = meta.groupby('MLST').count()['CRISPRtype'] >= 1
    # mlst_map = mlst_map.loc[data_MLST.index]
    # mlst_map = mlst_map.astype(int) #change to numeric for palette
    # my_palette = plt.cm.get_cmap("Accent", 2) # two colors crispr or not
    # # my_color = mlst_map.codes
# %%
    # color mapping for single crispr types
    #get df with true false for crispr type with types as index and mlst as cols
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
    #simplify crispr types
    # mlst_map_blanko = [x if x == "No_crispr" else "" for x in mlst_map]
    # mlst_map_blanko = pd.Series(mlst_map_blanko, index=mlst_map.index)
    for mlst in mlst_map.index:
        if mlst_map[mlst] == "No_crispr":
            pass
        elif "'I-" in mlst_map[mlst] and "'II-" in mlst_map[mlst]:
            mlst_map[mlst] = "I + II"
        elif "'I-" in mlst_map[mlst] and "'III-" in mlst_map[mlst]:
            mlst_map[mlst] = "I + III"
        elif "'II-" in mlst_map[mlst] and "'III-" in mlst_map[mlst]:
            mlst_map[mlst] = "II + III"
        elif "'I-" in mlst_map[mlst] and "'IV-" in mlst_map[mlst]:
            mlst_map[mlst] = "I + IV"
        elif "'IV-" in mlst_map[mlst]:
            mlst_map[mlst] = "IV"
        elif "'III-" in mlst_map[mlst]:
            mlst_map[mlst] = "III"
        elif "'II-" in mlst_map[mlst]:
            mlst_map[mlst] = "II"
        elif "'I-" in mlst_map[mlst]:
            mlst_map[mlst] = "I"
        else:
            print("No hardcoded case!")
            pass
    # mlst_map = mlst_map_blanko
    unique_combis_local = list(set(mlst_map))
    unique_combis_local.remove("No_crispr")
    unique_combis_local = ["No_crispr"] + unique_combis_local
    unique_combis_local = {unique_combis_local[i]:i for i in range(len(unique_combis_local))}


    # my_palette = colors.ListedColormap(colors_rgb[1:len(unique_combis)+1])

    _linkages = ['ward']
    labellist = ["{:>35}".format(l) for l in data_MLST.index.to_list()]
    for _linkage in _linkages:
        linked = linkage(data_MLST_nonsquare, method=_linkage)
        with plt.rc_context({'lines.linewidth': 2.5}):
            # fig, ax = plt.subplots(figsize=(5, 20))
            ax = axs[1, idx]
            # plt.title("Dendrogram Jaccard dissimilarity of {} with {} linkage".format(
            #     s, _linkage))
            dendrogram(linked,
                    orientation='left',
                    leaf_rotation=0,
                    distance_sort='descending',
                    show_leaf_counts=True, 
                    labels=labellist, 
                    ax=ax)

            # Apply the right color to each label
            # ax = plt.gca()
            ylbls = ax.get_ymajorticklabels()
            for lbl in ylbls:
                label = lbl.get_text()
                # lbl.set_backgroundcolor(my_palette(mlst_map.loc[label]))
                # feed label to map to dict, to get number for palette
                # lbl.set_backgroundcolor(my_palette(
                #     unique_combis[mlst_map.loc[label.strip()]]))
                lbl.set_size(3)
                lbl.set_bbox(dict(pad=0.0,
                            facecolor=my_palette(unique_combis[mlst_map.loc[label.strip()]]),
                            edgecolor=my_palette(unique_combis[mlst_map.loc[label.strip()]]))
                            )
            # legend_elements = [Patch(facecolor=my_palette(0), label='No Crispr (<1)'),
            #                    Patch(facecolor=my_palette(1), label='Has Crispr (>=1)')]
            legend_elements = [Patch(facecolor=my_palette(unique_combis[key]), label=key) for key in unique_combis_local.keys()]
            # ax.legend(handles=legend_elements, loc='upper left')
            #plt.show()
            fig.savefig("output/dendro_{}_{}-linkage.svg".format(s, _linkage), dpi=300)

    # get phylogenetic trees and compute correlation
    ab_phylogenic_distance = pd.read_csv(
        "phylogenetic trees/{}_dist_matrix.tsv".format(s[1:]), sep="\t", index_col=0)
    ab_phylogenic_distance = ab_phylogenic_distance.drop_prefix("{}_".format(s[1:])).T
    ab_phylogenic_distance = ab_phylogenic_distance.drop_prefix("{}_".format(s[1:])).T
    # cast everything to int
    ab_phylogenic_distance.index = ab_phylogenic_distance.index.astype(int)
    ab_phylogenic_distance.columns = ab_phylogenic_distance.columns.astype(int)
    data_MLST.index = data_MLST.index.astype(int)
    data_MLST.columns = data_MLST.columns.astype(int)

    common_mlst = [mlst for mlst in data_MLST.columns.tolist()
                   if mlst in ab_phylogenic_distance.columns.tolist()]
    common_mlst = list(map(int, common_mlst))
    # get the same index and overlapping mlsts
    ab_phylogenic_distance = ab_phylogenic_distance.loc[common_mlst,
                                                        common_mlst]
    data_MLST = data_MLST.loc[common_mlst,common_mlst]
    # flatten, flat mlst data for jaccard is already available (data_MLST_nonsquare)
    ab_phylogenic_distance_flat = squareform(ab_phylogenic_distance)
    data_MLST_nonsquare = squareform(data_MLST)

    _pearsonr = pearsonr(data_MLST_nonsquare,
                           ab_phylogenic_distance_flat)
    _spearmanr = spearmanr(data_MLST_nonsquare,
                             ab_phylogenic_distance_flat)
    _kendalltau = kendalltau(data_MLST_nonsquare,
                               ab_phylogenic_distance_flat)
    correls[s] = [_pearsonr[0],_spearmanr[0],_kendalltau[0]]
    
    # tree for phylogenic distances
    _linkages = ['ward']
    labellist = ["{:>35}".format(l)
                 for l in ab_phylogenic_distance.index.to_list()]
    for _linkage in _linkages:
        linked = linkage(ab_phylogenic_distance_flat, method=_linkage)
        with plt.rc_context({'lines.linewidth': 2.5}):
            # fig, ax = plt.subplots(figsize=(5, 20))
            ax = axs[0, idx]
            # plt.title("Dendrogram phylogenic distances of {} with {} linkage".format(
                # s, _linkage))
            dendrogram(linked,
                    orientation='left',
                    leaf_rotation=0,
                    distance_sort='descending',
                    show_leaf_counts=True,
                    labels=labellist, 
                    ax=ax)

            # Apply the right color to each label
            # ax = plt.gca()
            ylbls = ax.get_ymajorticklabels()
            for lbl in ylbls:
                label = lbl.get_text()
                # lbl.set_backgroundcolor(my_palette(mlst_map.loc[label]))
                # feed label to map to dict, to get number for palette
                # lbl.set_backgroundcolor(my_palette(
                #     unique_combis[mlst_map.loc[label.strip()]]))
                lbl.set_size(3)
                lbl.set_bbox(dict(pad=0.0,
                            facecolor=my_palette(unique_combis[mlst_map.loc[label.strip()]]),
                            edgecolor=my_palette(unique_combis[mlst_map.loc[label.strip()]]))
                            )            # legend_elements = [Patch(facecolor=my_palette(0), label='No Crispr (<1)'),
            #                    Patch(facecolor=my_palette(1), label='Has Crispr (>=1)')]
            legend_elements = [Patch(facecolor=my_palette(
                unique_combis[key]), label=key) for key in unique_combis_local.keys()]
            # ax.legend(handles=legend_elements, loc='upper left')
            #plt.show()

for ax, col in zip(axs[1], species):
    ax.set_title("$\it{}$".format(names[col]), size='xx-large')

for ax, row in zip(axs[:,0], ["Phylogeny (MLST)", "Gene profile"]):
    ax.set_ylabel(row, rotation=90, size='xx-large')

legend_elements = [Patch(facecolor=my_palette(
                unique_combis[key]), label=key) for key in unique_combis.keys()]
axs[0, 0].legend(handles=legend_elements, loc='lower left')

fig.savefig("output/fig_supplTrees.png".format(s, _linkage), dpi=300)
fig.savefig("output/fig_supplTrees.svg".format(s, _linkage), dpi=300)
fig.savefig("output/fig_supplTrees.tiff".format(s, _linkage), dpi=300)
# save correls
correls.to_csv("output/correlations.csv")
