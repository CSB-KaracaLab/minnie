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
Filtering functions.
"""

import glob
import logging
import os
import sys

import pandas as pd

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


# inside pathx (MD)
def time_freq_filter(filex, project_id, per):
    pathx = os.getcwd()
    file = os.path.basename(filex)
    fname = project_id
    bondtype = file.split(".csv")[0].split('_merged_')[1]
    try:
        first = pd.read_csv(filex)
    except:
        sys.exit(1)

    os.makedirs(f'{pathx}/{project_id}/03_time_freq_filter', exist_ok=True)
    pathxx = f'{pathx}/{project_id}/03_time_freq_filter'
    os.chdir(pathxx)

    pathy = pathxx + "/" + str(per) + "_freq_filtered"
    os.makedirs(str(per) + "_freq_filtered", exist_ok=True)
    os.chdir(pathy)

    if first.empty:
        logging.info(f'{first} is empty ..')
    else:

        # fIRST
        logging.info('Finding percentages: {}'.format(fname))
        firstx = []
        for adx in first.donor_acceptor.unique():
            bbx = first[first["donor_acceptor"] == adx]
            bbx_filt = bbx.time.unique().size
            first_filt = first.time.unique().size
            firstx.append([adx, bbx_filt / first_filt * 100])

        firstxy = pd.DataFrame(firstx)
        firstxy.columns = ["donor_acceptor", "percentage"]

        logging.info('Writing to file percentage: {}'.format(fname))
        morefirstxy = firstxy[firstxy.percentage > float(per)]

        if len(morefirstxy.donor_acceptor) == 0:
            pathz = pathy + "/" + str(per) + "_freq"
            os.makedirs(str(per) + "_freq", exist_ok=True)
            os.chdir(pathz)

            morefirstxy = pd.DataFrame(columns=firstxy.columns)
            morefirstxy.to_csv(pathz + "/" + fname + "_" + f'{bondtype}' +
                               "_" + str(per) + "_freq.csv", index=None)
            os.chdir("..")

            os.makedirs(str(per) + "_freq_perres", exist_ok=True)
            pathq = pathy + "/" + str(per) + "_freq_perres"
            os.chdir(pathq)

            first_perres = pd.DataFrame(columns=first.columns)
            first_perres.to_csv(pathq + "/" + fname + "_" + f'{bondtype}' +
                                "_" + str(per) + "_freq_perres.csv",
                                index=None)
        else:

            pathz = pathy + "/" + str(per) + "_freq"
            os.makedirs(str(per) + "_freq", exist_ok=True)
            os.chdir(pathz)

            morefirstxy.to_csv(pathz + "/" + fname + "_" + f'{bondtype}' + "_"
                               + str(per) + "_freq.csv", index=None)
            logging.info('Writing to file list: {}'.format(fname))

            first_perres = pd.DataFrame()
            for da in morefirstxy.donor_acceptor.unique():
                df = first[first.donor_acceptor == da]
                first_perres = first_perres.append(df)

            first_perres.sort_values(by="time", inplace=True)
            first_perres.reset_index(drop=True)

            os.chdir("..")

            os.makedirs(str(per) + "_freq_perres", exist_ok=True)
            pathq = pathy + "/" + str(per) + "_freq_perres"
            os.chdir(pathq)

            first_perres.to_csv(pathq + "/" + fname + "_" + f'{bondtype}' +
                                "_" + str(per) + "_freq_perres.csv",
                                index=None)


def make_freq_folders(pathy, per):
    """
    Creates folders to write and read common and complex-specific bonds
    within 04_compare_cx_spp folder

    """
    import os
    os.chdir(pathy)
    pathz = pathy + "/" + str(per) + "_freq_filtered"
    os.makedirs(str(per) + "_freq_filtered", exist_ok=True)

    for fold in ["_freq", "_freq_perres"]:
        os.chdir(pathz)
        # to add freq
        pathq = pathz + "/" + str(per) + fold
        os.makedirs(str(per) + fold, exist_ok=True)
        os.chdir(pathq)

        # pathq_common = pathq + "/common"
        if not os.path.exists("common"):
            os.makedirs("common", exist_ok=True)

        os.chdir(pathq)
        # pathq_spp = pathq + "/complex_specific"
        os.makedirs("complex_specific", exist_ok=True)


def get_paths(pathy, per, fold, com_spp):
    import os
    os.chdir(pathy)
    pathtowrite = \
        pathy + "/" + per + "_" + "freq_filtered/" + per + fold + "/" + com_spp
    return pathtowrite


def compare_bonds(project_id, per):
    def find_specific_bonds(name, morefirstxy_df, moresecxy_df, patho):
        """
        To find and write bonds specific to name complex.
        Its bond network is saved in morefirstxy_df.
        moresecxy_df shows bond network of given second complex.
        morefirstxy_df will be compared with moresecxy_df.
        """
        logging.info("Running on {} - {}".format(bondtype, name))
        i = 0
        spp_first = pd.DataFrame(columns=morefirstxy_df.columns)
        common_first = pd.DataFrame(columns=morefirstxy_df.columns)
        for item in morefirstxy_df.donor_acceptor:
            item_swapped = item.split(":")[1] + ":" + item.split(":")[0]
            if item in moresecxy_df.donor_acceptor.unique():
                df = pd.DataFrame(morefirstxy_df.iloc[i, :]).transpose()
                common_first = common_first.append(df)
            elif item_swapped in moresecxy_df.donor_acceptor.unique():
                df = pd.DataFrame(morefirstxy_df.iloc[i, :]).transpose()
                common_first = common_first.append(df)
            else:
                df = pd.DataFrame(morefirstxy_df.iloc[i, :]).transpose()
                spp_first = spp_first.append(df)
            i = i + 1

        spp_first.sort_values(by="donor_acceptor", ascending=False)
        spp_first.reset_index(drop=True, inplace=True)

        fold = "_freq"
        com_spp = "complex_specific"
        pathq_spp = get_paths(patho, str(per), fold, com_spp)
        towr1 = f'{pathq_spp}/{name}_{bondtype}_compared_spec.csv'
        spp_first.to_csv(towr1, index=False)
        common_first.sort_values(by="donor_acceptor", ascending=False)
        common_first.reset_index(drop=True, inplace=True)
        com_spp = "common"
        pathq_common = get_paths(patho, str(per), fold, com_spp)
        towr2 = f'{pathq_common}/{name}_{bondtype}_compared_common.csv'
        common_first.to_csv(towr2, index=False)
        return spp_first, common_first

    # write  bond list specific to xx one
    def list_specific_bonds(name, first_df, spp_first_l, common_first, patho):

        fold = "_freq_perres"
        spp_first_perres = pd.DataFrame(columns=first_df.columns)
        it_over = [["complex_specific", spp_first_l], ["common", common_first]]
        for llist in it_over:
            spp_first_l = llist[1]
            for item in spp_first_l.donor_acceptor.unique():
                f = first_df[first_df.donor_acceptor == item]
                spp_first_perres = spp_first_perres.append(f)

            spp_first_perres.sort_values(by="time", ascending=False)
            spp_first_perres.reset_index(drop=True, inplace=True)

            pathl_spp = get_paths(patho, str(per), fold, llist[0])
            if llist[0] == "complex_specific":
                spp = "spec"
            else:
                spp = "common"
            logging.info("{} list to {}".format(llist[0], name))
            spp_first_perres.to_csv(pathl_spp + "/" + name + "_" +
                                    f'{bondtype}_compared_{spp}_perres.csv',
                                    index=False)

    pathx = os.getcwd()
    fname = project_id[0]
    sname = project_id[1]

    file_lists_freq_fname = glob.glob(
        f'{pathx}/{fname}/03_time_freq_filter/{str(per)}_freq_filtered/'
        f'{str(per)}_freq/*csv')
    file_lists_freq_sname = glob.glob(
        f'{pathx}/{sname}/03_time_freq_filter/{str(per)}_freq_filtered/'
        f'{str(per)}_freq/*csv')

    file_lists_freq = file_lists_freq_fname + file_lists_freq_sname

    tocompare = {}
    for filex in file_lists_freq:
        file = os.path.basename(filex)

        if f'{fname}_' in filex:
            name = fname
        else:
            name = sname
        try:
            bondtype = file.split(f'{name}_')[1].split('_')[0]
        except:
            pass
        if bondtype == "ring":
            bondtype = "ring_stacking"
        first = pd.read_csv(filex)
        if bondtype in tocompare.keys():
            tocompare[bondtype].update({name: first})
        else:
            tocompare.update({bondtype: {name: first}})

    for bondtype in tocompare.keys():
        os.chdir(pathx)
        pathy = f'{pathx}/{fname}/04_compare_complex'
        os.makedirs(f'{pathx}/{fname}/04_compare_complex', exist_ok=True)

        pathz = f'{pathx}/{sname}/04_compare_complex'
        os.makedirs(f'{pathx}/{sname}/04_compare_complex', exist_ok=True)

        # FIRST complex
        make_freq_folders(pathy, per)
        morefirstxy = tocompare[bondtype][fname]

        fold = "_freq_perres"
        patha = f'{pathx}/{fname}/03_time_freq_filter/' \
                f'{str(per)}_freq_filtered/{str(per)}{fold}'
        first = pd.read_csv(patha + "/" + fname + "_" + f'{bondtype}'
                            + "_" + str(per) + fold + ".csv")

        # SECOND complex
        make_freq_folders(pathz, per)
        moresecxy = tocompare[bondtype][sname]
        logging.info("sname : {}".format(sname))

        fold = "_freq_perres"
        pathb = f'{pathx}/{sname}/03_time_freq_filter/' \
                f'{str(per)}_freq_filtered/{str(per)}{fold}'
        sec = pd.read_csv(pathb + "/" + sname + "_" + f'{bondtype}' +
                          "_" + str(per) + fold + ".csv")

        # find and write bonds specific to first one
        spp_first, common_first = \
            find_specific_bonds(fname, morefirstxy, moresecxy, pathy)

        # find and write bonds specific to second one
        spp_sec, common_sec = \
            find_specific_bonds(sname, moresecxy, morefirstxy, pathz)

        # write  bond list specific to first one
        list_specific_bonds(fname, first, spp_first, common_first, pathy)
        # write  bond list specific to second one
        list_specific_bonds(sname, sec, spp_sec, common_sec, pathz)
