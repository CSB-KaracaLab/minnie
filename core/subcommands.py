#!/usr/bin/env python

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
Command-line subcommand executors for minnie top-level script.
"""

import logging
import os
import sys

from .parallel import (
    create_mp_pool,
    parallelize
)

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


def splitpdbs(args):
    """Executes splitpdbs subcommand.

    Arguments
    ---------
        args: argparse.Namespace
            Result of parsing user arguments with argparse.
    """

    from .analysis import split_pdbs

    # Validate files exist
    args.pdbs = [
        f.resolve(strict=True) for f in args.pdbs
    ]

    # Validate project ids
    if args.project_ids:
        assert len(args.pdbs) == len(args.project_ids), \
            'Number of trajectories does not match number of ids: ' + \
            f'{len(args.pdbs)} != {len(args.project_ids)}'
    else:
        logging.info(f'{"Auto-assigning project ids from input file names"}')
        args.project_ids = [
            f.stem for f in args.pdbs
        ]

    for pdbfile, project_id in zip(args.pdbs, args.project_ids):
        split_pdbs(pdbfile, project_id)


def findbonds(args):
    """Executes findbonds command.

    Arguments
    ---------
        args: argparse.Namespace
            Result of parsing user arguments with argparse.
    """

    from .analysis import (
        comb_int,
        combine_interfacea_results
    )

    _itypes = [
        'hbonds',
        'ionic',
        'hydrophobic',
        'ring_stacking'
    ]
    if args.pdbfile:  # single structure
        pdblist = [args.pdbfile.resolve(strict=True)]
    elif args.folder:
        folder = args.folder.resolve(strict=True)
        pdblist = [
            f for f in folder.rglob('*.pdb')
        ]
    else:
        raise NotImplementedError(
            f'input namespace requires "pdbfile" or "folder" parameter'
        )

    if not args.project_id:
        logging.info(f'Auto-assigning project id from input file names')
        args.project_id = (args.pdbfile or args.folder).stem

    if 'all' in args.itypes:
        args.types = _itypes
    else:
        args.types = args.itypes

    mp_pool = create_mp_pool(args.nproc)
    for itype in args.types:
        parallelize(
            mp_pool,
            comb_int,
            pdblist,
            project_id=args.project_id,
            itypes=itype,
            include_intra=args.intra
        )
    if len(pdblist) != 1:
        combine_interfacea_results(
            args.project_id,
            args.clean
        )

    mp_pool.close()


def timefilter(self):
    """Apply critical interaction filter"""
    from .filtering import (
        time_freq_filter
    )

    if not self.project_id:
        logging.info(f'{"Auto-assigning project id from input file names"}')
        self.project_id = (self.pdbfile or self.folder).stem

    path_back = os.getcwd()

    if not self.files:
        print(f'\n{"where is the file(s) ??"}\n')
    elif not self.per:
        print(f'\n{"Please specify a cutoff value to filter out bonds !!"}\n')

    if self.per:
        for filex in self.files:
            os.chdir(path_back)
            time_freq_filter(filex, self.project_id, self.per)
    else:
        print("Please give a cutoff value")


def comparecx(self):
    """Calculate common and distinct interactions in two cases"""
    from .filtering import (
        compare_bonds
    )

    compare_bonds(self.project_ids, self.per)


def graph(self):
    """aaaaand graphs!"""
    from .graphs import (
        filter_todnaall,
        filter_todraw,
        draw_fig
    )
    if self.between is not None and self.chainIDs is not None:
        print("\nPlease specify either chainIDs or betweenness.")
        sys.exit(1)
    elif self.between is None and self.chainIDs is None:
        print("\nPlease specify either chainIDs or betweenness.")
        sys.exit(1)
    elif self.itypes == "all":
        for itypesx in ["hbonds", "ionic", "hydrophobic", "ring_stacking"]:
            if self.between:
                print(itypesx)
                df_collec = filter_todnaall(self.project_ids,
                                            self.between, self.spp,
                                            self.per, str(itypesx))
                draw_fig(df_collec, str(itypesx), self.project_ids[0],
                         self.project_ids[1],
                         self.colors[0], self.colors[1], self.filename,
                         self.spp)
            elif self.chainIDs:
                df_collec = filter_todraw(self.project_ids, self.chainIDs,
                                          self.spp, self.per,
                                          str(itypesx))
                draw_fig(df_collec, str(itypesx), self.project_ids[0],
                         self.project_ids[1],
                         self.colors[0], self.colors[1], self.filename,
                         self.spp)
    elif self.between == "protein-dna":
        df_collec = filter_todnaall(self.project_ids, self.between,
                                    self.spp, self.per, self.itypes)
        draw_fig(df_collec, self.itypes, self.project_ids[0],
                 self.project_ids[1],
                 self.colors[0], self.colors[1], self.filename,
                 self.spp)
    elif self.between == "all":
        df_collec = filter_todnaall(self.project_ids, self.between,
                                    self.spp, self.per, self.itypes)
        draw_fig(df_collec, self.itypes, self.project_ids[0],
                 self.project_ids[1],
                 self.colors[0], self.colors[1], self.filename,
                 self.spp)
    else:
        df_collec = filter_todraw(self.project_ids, self.chainIDs,
                                  self.spp, self.per, self.itypes)
        draw_fig(df_collec, self.itypes, self.project_ids[0],
                 self.project_ids[1],
                 self.colors[0], self.colors[1], self.filename,
                 self.spp)
