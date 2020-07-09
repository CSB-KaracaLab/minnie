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
Analysis functions.
"""

import logging
import pathlib

import interfacea as ia
import pandas as pd

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


def split_pdbs(pdbfile, project_id):
    """Splits a multi-model PDB file into individual PDB files.

    Arguments
    ---------
        pdbfile : pathlib.Path
            input PDB file to process.
        project_id : str
            identifier for this complex.
    """

    logging.info(f'Splitting {project_id} ...')

    curdir = pathlib.Path('.').resolve(strict=True)
    output_dir = curdir / project_id / '01_frames'
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(pdbfile, 'r') as pdblines:
        model = []
        model_no = 0
        for line in pdblines:
            if line.startswith(('ATOM', 'HETATM', 'TER')):
                model.append(line)
            elif line.startswith('ENDMDL'):
                fpath = output_dir / f'md_{model_no}.pdb'
                with open(fpath, 'w') as outfile:
                    model.append('END')
                    print(''.join(model), file=outfile)
                model = []
                model_no += 1

    logging.info(f'Read {model_no} models from input file.')


def comb_int(pdbfile, project_id, itype, include_intra=False, **kwargs):
    """Analyzes the interactions in one or more PDB files.

    Arguments
    ---------
        ....
    """

    curdir = pathlib.Path('.').resolve(strict=True)
    output_dir = curdir / str(
        project_id[0]) / '02_interfacea_results' / f'{itype}'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup interfacea
    ia.set_log_level('minimal')
    mol = ia.read(str(pdbfile))  # convert from Path to str.
    analyzer = ia.InteractionAnalyzer(mol)

    func = getattr(analyzer, f'get_{itype}', None)
    if func is None:
        raise ValueError(
            f'Unknown analysis function: get_{itype}'
        )

    kwargs['include_intra'] = include_intra
    func(**kwargs)

    df_table = analyzer.itable._table
    df_table.columns = [
        'itype', 'donor_chain', 'acceptor_chain', 'donor_resnm',
        'acceptor_resnm', 'donor_resid', 'acceptor_resid',
        'donor_atom', 'acceptor_atom'
    ]

    if len(df_table) == 0:
        logging.info(
            f'No interactions of type "{itype}" found in input file'
        )
        return

    donor_list = df_table.apply(
        lambda x: x['donor_resnm'] + str(x['donor_resid']),
        axis=1
    )
    df_table.loc[:, 'donor'] = donor_list

    donorC_list = df_table.apply(
        lambda x, sep="_": x['donor'] + "_" + str(x['donor_chain']),
        axis=1
    )
    df_table.loc[:, 'donorC'] = donorC_list

    acceptor_list = df_table.apply(
        lambda x: x['acceptor_resnm'] + str(x['acceptor_resid']),
        axis=1
    )
    df_table.loc[:, 'acceptor'] = acceptor_list

    acceptorC_list = df_table.apply(
        lambda x, sep="_": x['acceptor'] + "_" + str(x['acceptor_chain']),
        axis=1
    )
    df_table.loc[:, 'acceptorC'] = acceptorC_list

    donor_acceptor_list = df_table.apply(
        lambda x, sep=":": x['donorC'] + sep + str(x['acceptorC']),
        axis=1
    )
    df_table.loc[:, 'donor_acceptor'] = donor_acceptor_list

    chain_type = df_table.apply(
        lambda x: (
            "intra" if (x["acceptor_chain"] == x["donor_chain"])
            else "inter"
        ),
        axis=1
    )
    df_table.loc[:, 'chain_type'] = chain_type

    # You could probably move this to a data.py module.
    # Made them sets since you are doing 'x in y' operations.
    protein_residues = {
        'ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR',
        'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER',
        'TRP', 'VAL'
    }
    dna_residues = {'DA', 'DC', 'DG', 'DT'}

    prot_or_dna = []
    for i in range(len(df_table.donor_resnm)):
        d_is_dna = df_table.donor_resnm[i] in dna_residues
        d_is_prot = df_table.donor_resnm[i] in protein_residues
        a_is_dna = df_table.acceptor_resnm[i] in dna_residues
        a_is_prot = df_table.acceptor_resnm[i] in protein_residues

        if d_is_dna and a_is_dna:
            cmplx_type = 'dna-dna'
        elif d_is_prot and a_is_prot:
            cmplx_type = 'protein-protein'
        elif (d_is_dna and a_is_prot) or (d_is_prot and a_is_dna):
            cmplx_type = 'protein-dna'
        else:
            cmplx_type = 'other'  # because why not?

        prot_or_dna.append(cmplx_type)
    df_table["prot_or_dna"] = prot_or_dna

    # renumber index of hbond dataframe
    df_table.index = list(range(len(df_table)))

    if itype == 'hbonds':
        non_spp_atoms = {
            "O2P", "O1P", "N", "O", "OC1", "OC2",
            "O4'", "O5'", "O3'", "H", "HA"
        }
        specificity = []
        for i in df_table.index:
            specificity_acc = df_table["acceptor_atom"][i] in non_spp_atoms
            specificity_donor = df_table["acceptor_atom"][i] in non_spp_atoms
            if specificity_acc or specificity_donor:
                specificity.append("non-specific")
            else:
                specificity.append("specific")
        df_table['specificity'] = specificity

    frame_no = pdbfile.stem.split('md_')[1]
    df_table["time"] = [frame_no] * len(df_table)

    logging.info(f'Writing "{itype}" bonds to disk...')
    df_table.to_csv(
        str(output_dir / f'md_{frame_no}.pdb_{itype}_all.csv'),
        index=False
    )


def combine_interfacea_results(project_id, clean=False):
    curdir = pathlib.Path('.').resolve(strict=True)
    output_dir = curdir / str(project_id[0]) / '02_interfacea_results'
    output_dir = output_dir.resolve(strict=True)

    for child in output_dir.iterdir():
        if child.is_dir():
            logging.info(f'processing csv files in {child} ..')

            dfx = pd.DataFrame()
            for csvfile in child.rglob('*all.csv'):  # recursive find!
                df = pd.read_csv(csvfile)
                dfx = dfx.append(df)  # very slow!

            logging.info(f'writing combined analysis to {child}')
            dfx.to_csv(
                str(child / f'{project_id[0]}_merged_{child.stem}.csv'),
                index=False
            )

    # Remove intermediate files
    if clean:
        logging.info('removing intermediate files..')
        for csvfile in output_dir.rglob('*all.csv'):
            csvfile.unlink()
