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
minnie: a structural ensemble analysis package to deduce the fingerprints of
binding
"""

import argparse
import logging
import pathlib
import sys

from core.gui import (
    guivar
)
from core.subcommands import (  # eventually this should be minnie.core
    splitpdbs,
    findbonds,
    timefilter,
    comparecx,
    graph
)

# Setup top-level logger
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y/%m/%d %H:%M:%S'
)

SUBCOMMANDS = {
    'splitpdbs': splitpdbs,
    'findbonds': findbonds,
    'timefilter': timefilter,
    'comparecx': comparecx,
    'graph': graph
}

if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='minnie', description=__doc__)
    parser.add_argument(
        '--nproc',
        default=None,
        type=int,
        help='Number of processes to use for parallel execution.'
    )

    # Setup subparsers
    subparsers = parser.add_subparsers(dest="subcommand", help="")

    # minnie splitpdb
    splitpdb = subparsers.add_parser(
        'splitpdbs',
        help='Split a trajectory into single frames',
        usage=guivar.splitpdbs[0]
    )
    splitpdb.add_argument(
        '--pdbs',
        nargs='+',
        type=pathlib.Path,
        dest='pdbs',
        help=argparse.SUPPRESS
    )
    splitpdb.add_argument(
        '-i'  # simpler than -cn
        '--id',  # simpler than --complexName
        nargs='+',
        type=str,
        dest='project_ids',
        help=argparse.SUPPRESS,
    )

    # minnie findbonds
    findbonds = subparsers.add_parser(
        'findbonds',
        help='Calculates interactions between and/or within monomers',
        description="Executes findbonds command",
        usage=guivar.findbonds[0]
    )

    findbonds.add_argument(
        '-i'
        '--id',
        nargs='+',
        type=str,
        dest='project_id',
        help=argparse.SUPPRESS
    )
    inputopts = findbonds.add_mutually_exclusive_group(required=True)
    inputopts.add_argument(
        '-f',
        '--pdbfile',
        type=pathlib.Path,
        help=argparse.SUPPRESS,
    )
    inputopts.add_argument(
        '-d',
        '--folder',
        type=pathlib.Path,
        help=argparse.SUPPRESS,
    )
    findbonds.add_argument(
        '--itypes',
        choices=[
            'hbonds',
            'ionic',
            'hydrophobic',
            'ring_stacking',
            'all'
        ],
        nargs='+',
        default=['hbonds'],
        type=str,
        dest='itypes',
        help=argparse.SUPPRESS,
    )
    findbonds.add_argument(
        '--intra',
        action='store_true',
        help=argparse.SUPPRESS,
    )
    findbonds.add_argument(
        '--clean',
        action='store_true',
        help=argparse.SUPPRESS,
    )

    timefilter = subparsers.add_parser(
        'timefilter',
        help='Apply critical interaction filter',
        usage=guivar.timefilter[0]
    )
    timefilter.add_argument(
        '-f',
        '--files',
        nargs='+',
        type=pathlib.Path,
        dest='files',
        help=argparse.SUPPRESS,
    )
    timefilter.add_argument(
        '-i'  # simpler than -cn
        '--id',  # simpler than --complexName
        type=str,
        dest='project_id',
        help=argparse.SUPPRESS,
    )
    timefilter.add_argument(
        '--per',  # simpler than --complexName
        dest='per',
        type=int,
        required=True,
        help=argparse.SUPPRESS,
    )
    comparecx = subparsers.add_parser(
        'comparecx',
        help='Calculate common and distinct interactions in two cases',
        usage=guivar.comparecx[0]
    )

    comparecx.add_argument(
        '-i'  # simpler than -cn
        '--id',  # simpler than --complexName
        nargs='+',
        type=str,
        dest='project_ids',
        help=argparse.SUPPRESS,
    )
    comparecx.add_argument(
        '--per',  # simpler than --complexName
        dest='per',
        type=int,
        help=argparse.SUPPRESS,
    )
    graph = subparsers.add_parser(
        'graph',
        help='aaaaand graphs!',
        usage=guivar.graph[0]
    )

    graph.add_argument(
        '-i'  # simpler than -cn
        '--id',  # simpler than --complexName
        nargs='+',
        type=str,
        dest='project_ids',
        help=argparse.SUPPRESS,
    )
    graph.add_argument(
        '--per',  # simpler than --complexName
        dest='per',
        type=int,
        help=argparse.SUPPRESS,
    )
    graph.add_argument(
        '-b'  # simpler than -cn
        '--between',  # simpler than --complexName
        dest='between',
        help=argparse.SUPPRESS,
    )
    graph.add_argument(
        '-c'  # simpler than -cn
        '--chainIDs',  # simpler than --complexName
        nargs='+',
        dest='chainIDs',
        help=argparse.SUPPRESS,
    )
    graph.add_argument(
        '--filename',  # simpler than --complexName
        dest='filename',
        help=argparse.SUPPRESS,
    )
    graph.add_argument(
        '--colors',  # simpler than --complexName
        dest='colors',
        default=['#D9B4CC', '#6F81A6'],
        help=argparse.SUPPRESS,
    )
    graph.add_argument(
        '--itypes',  # simpler than --complexName
        dest='itypes',
        help=argparse.SUPPRESS,
    )
    graph.add_argument(
        '-s',  # simpler than --complexName
        dest='spp',
        default="specific",
        help=argparse.SUPPRESS,
    )

    # Parse unknown args to function
    # from https://stackoverflow.com/a/37367814
    parsed, unknown = parser.parse_known_args()
    for argname in unknown:
        if argname.startswith(("-", "--")):
            findbonds.add_argument(argname)

    args = parser.parse_args()
    logging.info(f'{" ".join(sys.argv)}')

    func = SUBCOMMANDS.get(args.subcommand)
    if func is None:
        raise ValueError(f'Unknown subcommand: {args.subcommand}')

    func(args)  # execute subcommand
