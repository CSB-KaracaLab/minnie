# Copyright 2020 Deniz Dogan, Ezgi Karaca
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Plotting and visualization functions.
"""

import datetime
import itertools
import logging
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.legend_handler import HandlerPatch

pd.options.mode.chained_assignment = None

# Setup logger
# _private name to prevent collision/confusion with parent logger
#logging.getLogger(__name__).addHandler(logging.NullHandler())

now = datetime.datetime.now()
timestamp = now.strftime("%H.%M.%S")


def filter_todnaall(project_ids, between, spp, per, itypes):
    """
    To get interactions  at the protein-DNA interface or all
    """
    pathx = os.getcwd()
    if spp == "specific":
        spec = "spec"
        cx_spp = "complex_specific"
    elif spp == "common":
        spec = "common"
        cx_spp = "common"

    df_collec = pd.DataFrame(columns=["time", "len", "uniq", "group", "category", "hue"])
    for Name in project_ids:
        logging.info(f'Processing {Name} ...')
        ffile_path = f'{pathx}/{Name}/04_compare_complex/{str(per)}_freq_filtered/{str(per)}_freq_perres/{cx_spp}/{Name}_{itypes}_compared_{spec}_perres.csv'
        ffile = pd.read_csv(ffile_path)

        selected_ffile = pd.DataFrame()
        if between == "protein-dna":
            selected_ffile1 = ffile[(ffile.prot_or_dna == between)]
            selected_ffile1["category"] = itypes
            selected_ffile1["hue"] = itypes
            selected_ffile1["group"] = Name
            selected_ffile = selected_ffile.append(selected_ffile1)

        elif between == "all":
            selected_ffile = ffile
        else:
            raise ValueError

        d = pd.DataFrame(columns=["time", "len", "uniq", "group", "category", "hue"])
        for time in selected_ffile.time.unique():
            b = selected_ffile[selected_ffile.time == time]
            c = pd.DataFrame([time, len(b.donor_acceptor), len(b.donor_acceptor.unique()), Name, itypes, itypes]).T
            c.columns = ["time", "len", "uniq", "group", "category", "hue"]
            d = d.append(c)

        df_collec = df_collec.append(d)

    return df_collec


def filter_todraw(project_ids, chainIDs, spp, per, itypes):
    """
    To get filter data via chainIDs
    """

    pathx = os.getcwd()
    if spp == "specific":
        spec = "spec"
        cx_spp = "complex_specific"
    elif spp == "common":
        spec = "common"
        cx_spp = "common"

    fg = chainIDs[0].split("+")
    sg = chainIDs[1].split("+")
    chaingroups = list(itertools.product(fg, sg))

    df_collec = pd.DataFrame(columns=["time", "len", "uniq", "group", "category", "hue"])
    for Name in project_ids:
        logging.info(f'Processing {Name} ...')
        ffile_path = f'{pathx}/{Name}/04_compare_complex/{str(per)}_freq_filtered/{str(per)}_freq_perres/{cx_spp}/{Name}_{itypes}_compared_{spec}_perres.csv'

        ffile = pd.read_csv(ffile_path)

        selected_ffile = pd.DataFrame()
        for element in chaingroups:
            selected_ffile1 = ffile[(ffile.donor_chain == element[0]) & (ffile.acceptor_chain == element[1]) |
                                    (ffile.donor_chain == element[1]) & (ffile.acceptor_chain == element[0])]
            selected_ffile1["category"] = itypes
            selected_ffile1["hue"] = itypes
            selected_ffile1["group"] = Name
            selected_ffile = selected_ffile.append(selected_ffile1)

        d = pd.DataFrame(columns=["time", "len", "uniq", "group", "category", "hue"])
        for time in selected_ffile.time.unique():
            b = selected_ffile[selected_ffile.time == time]
            c = pd.DataFrame([time, len(b.donor_acceptor), len(b.donor_acceptor.unique()), Name, itypes, itypes]).T
            c.columns = ["time", "len", "uniq", "group", "category", "hue"]
            d = d.append(c)

        df_collec = df_collec.append(d)

    return df_collec


def draw_fig(df_hbond_collec, itypes, fName, sName, fcolor, scolor, filename, spp):
    if spp == "specific":
        spec = "spec"
        cx_spp = "complex_specific"
    elif spp == "common":
        spec = "common"
        cx_spp = "common"

    pathx = os.getcwd()
    pathy = f'{pathx}/graphs'
    if not os.path.exists(pathy):
        os.makedirs(pathy, exist_ok=True)

    colors = [fcolor, scolor]
    import matplotlib.patches as patches
    if itypes == "hbonds":
        itypes = "Hydrogen"
    else:
        itypes = itypes.capitalize()

    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(4, 4)

    rect = patches.Rectangle((0, 0), 0, 0, linewidth=1, edgecolor='black', facecolor='none')
    ax = fig.add_subplot(gs[:, :])
    ax.add_patch(rect)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel("Number of " + itypes + " Bonds", size=18)
    ax.yaxis.set_label_coords(-0.11, 0.5)

    ax_box1 = fig.add_subplot(gs[:, 0:2])

    boxprops = dict(linestyle='-', linewidth=1.5, edgecolor='#7f7f7f')
    medianprops = {"linestyle": '-', "linewidth": "1.5", "color": "#333333"}
    capprops = dict(linestyle='-', linewidth=1.5, color='#7f7f7f')
    whiskerprops = dict(linestyle='-', linewidth=1.5, color='#7f7f7f')

    order = [fName, sName]
    try:
        ax_box1.set_ylim(min(df_hbond_collec["uniq"]) - 1, max(df_hbond_collec["uniq"]) + 1)

        sns.boxplot(data=df_hbond_collec, x="group", y="uniq", ax=ax_box1, palette=colors, order=order,
                    medianprops=medianprops, boxprops=boxprops, capprops=capprops, whiskerprops=whiskerprops)

        ax_box1.set_ylabel("")
        ax_box1.set_xlabel("")
        ax_box1.spines['right'].set_visible(False)

        ax_hist1 = fig.add_subplot(gs[:, 2], sharey=ax_box1)
        ax_hist1.set_title(fName)
        ax_hist2 = fig.add_subplot(gs[:, 3], sharex=ax_hist1, sharey=ax_box1)
        ax_hist2.set_title(sName)

        df_hbond_collec_wt = df_hbond_collec[df_hbond_collec.group == fName]
        sns.distplot(df_hbond_collec_wt["uniq"], ax=ax_hist1, vertical=True, kde=False, rug=False,
                     hist_kws=dict(edgecolor="black"),
                     color=fcolor)
        df_hbond_collec_mut = df_hbond_collec[df_hbond_collec.group == sName]
        ax_hist1.spines['right'].set_visible(False)
        sns.distplot(df_hbond_collec_mut["uniq"], ax=ax_hist2, vertical=True, kde=False, rug=False,
                     hist_kws=dict(edgecolor="black"),
                     color=scolor)
        ax_hist1.set_ylabel("")
        ax_hist2.set_ylabel("")
        plt.setp(ax_hist1.get_yticklabels(), visible=False)
        plt.setp(ax_hist2.get_yticklabels(), visible=False)
        plt.setp(ax_box1.get_yticklabels(), visible=True)
        ax_box1.tick_params(axis="y", which=u'both', length=0)
        texts = [fName, sName]

        class HandlerEllipse(HandlerPatch):
            def create_artists(self, legend, orig_handle,
                               xdescent, ydescent, width, height, fontsize, trans):
                center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
                p = mpatches.Rectangle(xy=center, width=width + xdescent,
                                       height=height + ydescent, linewidth=1, edgecolor='b')
                self.update_prop(p, orig_handle, legend)
                p.set_transform(trans)
                return [p]

        c = [mpatches.Circle(0.5, 0.5, facecolor=colors[i], edgecolor='black', linewidth=1) for i in range(len(texts))]
        plt.legend(c, texts, bbox_to_anchor=(1.55, 0.95), loc='center',
                   handler_map={mpatches.Rectangle: HandlerEllipse()})
        plt.locator_params(axis="both", integer=True, tight=True)

        logging.info("Drawing graph for {} bonds...".format(itypes))

        if filename:
            if os.path.exists(f'{pathy}/{filename}_{spp}.png'):
                logging.info(f'Overwriting {pathy}/{filename}_{spec}.png ...')
            plt.savefig(f'{pathy}/{filename}_{spp}.png', bbox_inches='tight', dpi=300)
            plt.close("all")

        elif os.path.exists(f'{pathy}/{fName}_{sName}_{itypes}_{spec}.png'):
            plt.savefig(f'{pathy}/{fName}_{sName}_{itypes}_{spec}_{timestamp}.png', bbox_inches='tight', dpi=300)
            plt.close("all")
        else:
            plt.savefig(f'{pathy}/{fName}_{sName}_{itypes}_{spec}.png', bbox_inches='tight', dpi=300)
            plt.close("all")
    except ValueError:
        logging.info(f'No drawable {itypes} interactions..')
