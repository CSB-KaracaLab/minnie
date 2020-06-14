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

import pathos
from pathos.multiprocessing import ProcessingPool as Pool

from core.subcommands import (  # eventually this should be minnie.core
    splitpdbs,
    findbonds
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
}

# # Decorator to mark subcommands.
# def subcommand(parent=SUBPARSERS):

#     def decorator(func):
#         parser = SUBPARSERS.add_parser(
#             func.__name__,
#             help= func.__doc__ ,
#             conflict_handler='resolve'
#         )
#         parserg = parser.add_argument_group("parameters to give")
#         parserg.set_defaults(func=func)
#         return parserg

#     return decorator


# # Decorator to add options to subcommands
# def option(*args, **kwargs):
#     def decorator(parserg):
#         parserg.add_argument(*list(args), **kwargs)
#         return parserg
#     return decorator


# ###### ---- Subcommands ---- ######


# # ---- find bonds ---- #
# @option('--help','-h',action='store_true')
# @option('--systematic','-s',choices=["True","False"], help="Use this option if you want to analyze multiple pdbs in a folder",
#         default="True")
# @option('--pdbs','-p',nargs="*", help="Give single *.pdb or give folder path ")
# @option('--complexName','-cn', help="Project ID of your complex(s)")
# @option('-i', choices=["hbonds","ionic","hydrophobic","ring_stacking","all"],default="hbonds",
#                      dest="intType", help="choose interaction type")
# @option('-intra',"--includeIntra", choices=["True","False"],default=False,
#                      help="includes intramolecular residue pairs in the search. By default, False.")
# @option("-d", default=2.5, dest="hbond_distance",help=" cutoff value to find hbonds")
# @subcommand()
# def findbonds(options):
#     """Calculates interactions between and/or within monomers"""
#     if options.help:
#         print("Calculates interactions between and/or within monomers\n"
#             f'\n\033[1mUsage: minnie findbonds \n'
#               f'                        -cn, --complexName     <string>     \n '
#               f'                                               Project ID of your complex\n\n'

#               f'                        -p, --pdbs             [<.pdb>/<path>] (singleframe.pdb)   \n'
#               f'                                               Give single *.pdb or give folder path \n\n'

#               f'                        -i                     [<hbonds>/<ionic>/<hydrophobic>/<ring_stacking>/<all>] (hbonds)    \n'
#               f'                                               Calculates which types of interactions \n\n'

#               f'                        -d                      <float> (2.5)                 \n'
#               f'                                               Cut-off to define a hydrogen bond\n\n'

#               f'                        -intra, --includeIntra [<"True">/<"False">] ("False")  \n'
#               f'                                               What do you want to analyze, all or only inter-monomer contacts? \033[0m \n\n\n\n'

#               f'\n\033[1mUsage example:\033[0m\n\n'
#               " Single frame    - minnie findbonds -cn sox4 -p sox4/02_frames/md_0.pdb -i hbonds  -s False  \n"
#               " Multiple frames - minnie findbonds -cn sox4 -p sox4/02_frames/* -i hbonds \n"
#               " Multiple frames - minnie findbonds -cn sox4 -p sox4/02_frames/* -i all \n"

#               )
#     elif not options.pdbs:
#         print(f'where is pdb??')
#     elif not options.complexName:
#         print(f'Please specify complex name(s)')


#     elif (options.systematic) == "True":
#         pdb_list = options.pdbs
#         if (options.intType == "all"):
#             for intType in ["hbonds", "ionic", "hydrophobic", "ring_stacking"]:
#                 pool = Pool(pathos.multiprocessing.cpu_count() - 2)
#                 pool.map(analysis.comb_int, pdb_list, len(pdb_list) * [str(options.complexName)],
#                          len(pdb_list) * [str(intType)], len(pdb_list) * [str(options.includeIntra)],
#                          len(pdb_list) * [str(options.hbond_distance)])
#                 #pool.close()

#         else:
#             pool = pathos.multiprocessing.ProcessingPool(pathos.multiprocessing.cpu_count() - 2)
#             pool.map(analysis.comb_int, pdb_list, len(pdb_list) * [str(options.complexName)],
#                      len(pdb_list) * [str(options.intType)],len(pdb_list) * [str(options.includeIntra)],
#                      len(pdb_list) * [str(options.hbond_distance)] )
#             pool.close()
#         analysis.combine_interfacea_results(options.complexName)
#     elif (options.systematic) == "False":
#         if (options.intType == "all"):
#             for intType in ["hbonds", "ionic", "hydrophobic", "ring_stacking"]:
#                 analysis.comb_int(options.pdbs[0], options.complexName, intType, options.includeIntra, options.hbond_distance)
#         else:
#             analysis.comb_int(options.pdbs[0],options.complexName,options.intType,options.includeIntra,options.hbond_distance)


#         analysis.combine_interfacea_results(options.complexName)





# # ---- Time filtering ---- #
# @option('--help','-h',action='store_true')
# @option('--files','-f',nargs="*", help="Files")
# @option('--complexName','-cn',help="Project ID of your complex(s)")
# @option('--per', help="observation frequency (in %) to classify an interaction as critical?", type = int)
# @subcommand()
# def timefilter(options):
#     """Apply critical interaction filter"""
#     path_back = os.getcwd()
#     if options.help:
#         print(f'\n\033[1mUsage: minnie timefilter \n'
#               f'                        -cn, --complexName     <string>      \n '
#               f'                                               Project ID of your complex\n\n'

#               f'                         -f, --files            [<.csv>]               \n'
#               f'                                                Files of your complex\n\n'
#               f'                         --per                  <float>                \n'
#               f'                                                Observation frequency to classify an interaction as critical\033[0m \n\n\n\n'

#             f'\n\033[1mUsage example:\033[0m\n\n'
#               " Single file    - minnie timefilter -f sox4/03_interfacea_results/hbonds/sox4_merged_hbonds.csv -cn sox4 --per 25  \n"
#               " Multiple files - minnie timefilter -f sox4/03_interfacea_results/*/sox4_merged_*.csv -cn sox4 --per 25  \n"
#               )
#     elif not options.files:
#         print(f'\nwhere is the file(s) ??\n')
#     elif not options.complexName:
#         print(f'\nPlease specify a complex name(s) !!\n')
#     elif not options.per:
#         print(f'\nPlease specify a cutoff value to filter out bonds !!\n')

#     if (options.per):
#         for filex in options.files:
#             os.chdir(path_back)
#             filtering.time_freq_filter(filex,options.complexName,options.per)
#         #pool = pathos.multiprocessing.ProcessingPool(pathos.multiprocessing.cpu_count() - 2)
#         #pool.map(filtering.time_freq_filter, options.files,
#         #         len(options.files)*[options.complexName],
#         #         len(options.files)*[options.per])
#         #pool.close()


# # ---- Compare interaction networks between two complex ---- #
# @option('--help','-h',action='store_true')
# @option('--per', help="observation frequency (in %) to classify an interaction as critical?", type = int)
# @option('--complexName','-cn',nargs=2, help="Project ID of your complex(s)")
# @subcommand()
# def compareCX(options):
#     """Calculate common and distinct interactions in two cases"""
#     if options.help:
#         print(f'\n\033[1mUsage: minnie compareCX  \n'
#               f'                        -cn, --complexName     <string> <string>     \n '
#               f'                                               Project ID of your complex(s)\n\n'
#               f'                        --per                  <float>                \n'
#               f'                                                Observation frequency to classify an interaction as critical\033[0m \n\n\n\n'

#               f'\n\033[1mUsage example:\033[0m\n\n'
#               " minnie compareCX -cn sox4 sox18 --per 25  \n")
#     else:
#         filtering.compare_bonds(options.complexName,options.per)





# # ---- draw graphs! ---- #
# @option('--help','-h',action='store_true')
# @option('--per', help="observation frequency (in %) to classify an interaction as critical?", type = int)
# @option('--colors',nargs="*", help="Color IDs of the complexes (for plotting)",default=['#D9B4CC', '#6F81A6'])
# @option('--chainIDs','-c',nargs=2,dest="chains", help="Give chainID(s)")
# @option('--complexName','-cn',nargs=2, help="Project ID of your complex(s)")
# @option('-s', choices=["specific","common"],default="specific", dest="spp")
# @option('-i', choices=["hbonds","ionic","hydrophobic","ring_stacking","all"],default="hbonds",
#                      dest="intType", help="choose interaction type")
# @option('-b',"--between", choices=["protein-dna","all"], dest="between")
# @option('--filename', dest="filename",help="filename to use while writing graph",default="")
# @subcommand()
# def graph(options):
#     """aaaaand graphs!"""
#     try:
# 	    if options.help:
# 	        print(f'\n\033[1mUsage: minnie graph  \n'
# 	              f'                        -cn, --complexName     <string> <string>     \n'
# 	              f'                                               Project IDs of your complex(s)\n\n'

# 	              f'                        --per                  <float>                \n'
# 	              f'                                               Observation frequency to classify an interaction as critical \n\n'


# 	              f'                        -b, --between          [<protein-dna>/<all>] (all)   \n'
# 	              f'                                               Between protein-dna or keep all \n\n'

# 	              f'                        -c, --chainIDs         <string> <string>       \n'
# 	              f'                                               Give ChainIDs to proceed\n\n'

# 	              f'                        --filename             <string>                           \n'
# 	              f'                                               Give a name to output file (optional)\n\n'

# 	              f'                        --colors               [<hex colors>] ("#D9B4CC", "#6F81A6")     \n'
# 	              f'                                               Color IDs of the complexes (optional)\n\n'


# 	              f'                        -i                     [<hbonds>/<ionic>/<hydrophobic>/<ring_stacking>/<all>] (hbonds)    \n'
# 	              f'                                               Calculates which types of interactions \n\n'

# 	              f'                        -s                     [<specific>/<common>] (specific)                           \n'
# 	              f'                                               Complex-specific or common interactions\033[0m \n\n\n\n'

# 	            f'Please do not use "--between" and "--chainIDs" options at the same time\n\n'

# 	            f'\n\033[1mUsage example:\033[0m\n\n'
# 	              " minnie graph -cn 'sox4' 'sox18' --per 25 -i hbonds -s specific  -c A+B C --colors '#D9B4CC' '#6F81A6' \n"
# 	              " minnie graph -cn 'sox4' 'sox18' --per 25 -i ionic -c A+B C  \n"
# 	              " minnie graph -cn 'sox4' 'sox18' --per 25 -i ionic -b protein-dna \n"
# 	              " minnie graph -cn 'sox4' 'sox18' --per 25 -i ionic -b protein-dna --filename sox4_sox18_protein_dna \n")
# 	    elif (options.between) and (options.chains):
# 	        raise Exception()
# 	    elif options.intType == "all":
# 	        print(options.between)
# 	        for intTypex in ["hbonds","ionic","hydrophobic","ring_stacking"]:

# 	            if  (options.between):
# 	                print(intTypex)
# 	                df_collec=graphs.filter_todnaall(options.complexName,options.between,options.spp,options.per,str(intTypex))
# 	                graphs.draw_fig(df_collec, str(intTypex), options.complexName[0], options.complexName[1],
# 	                                options.colors[0], options.colors[1],options.filename, options.spp)

# 	            elif (options.chains):
# 	                df_collec=graphs.filter_todraw(options.complexName,options.chains,options.spp,options.per,str(intTypex))
# 	                graphs.draw_fig(df_collec, str(intTypex), options.complexName[0], options.complexName[1],
# 	                                        options.colors[0], options.colors[1],options.filename, options.spp)

# 	    elif options.between == "protein-dna":
# 	        df_collec=graphs.filter_todnaall(options.complexName,options.between,options.spp,options.per,options.intType)
# 	        graphs.draw_fig(df_collec, options.intType, options.complexName[0], options.complexName[1],
# 	                        options.colors[0], options.colors[1],options.filename, options.spp)

# 	    elif options.between == "all":
# 	        df_collec=graphs.filter_todnaall(options.complexName,options.between,options.spp,options.per,options.intType)
# 	        graphs.draw_fig(df_collec, options.intType, options.complexName[0], options.complexName[1],
# 	                        options.colors[0], options.colors[1],options.filename, options.spp)
# 	    else:
# 	        df_collec=graphs.filter_todraw(options.complexName,options.chains,options.spp,options.per,options.intType)
# 	        graphs.draw_fig(df_collec, options.intType, options.complexName[0], options.complexName[1],
# 	                        options.colors[0], options.colors[1],options.filename, options.spp)
#     except TypeError:
#         print(f'\nPlease check given parameters''')

#     except Exception:
#         print(f'\nPlease, either use -b or -c option''')



# # ---- Clean  ---- #
# @option('--help','-h',action='store_true')
# @option('--complexName','-cn',nargs=2, help="Project ID of your complex(s)")
# @subcommand()
# def clean(options):
#     """To remove unnecessary folders"""
#     if options.help:
#         print(f'\n\033[1mUsage: minnie clean \n'
#               f'                        -cn, --complexName     <string> <string>     \n '
#               f'                                               Project ID of your complex(s)\n\n'

#               f'\n\033[1mUsage example:\033[0m\n\n'
#               " minnie clean -cn sox4 sox18  \n")
#     else:
#         try:
#             core.clean.cleanx(options.complexName[0])
#         except FileNotFoundError:
#             print(f'Nothing to clean inside {options.complexName[0]}')

#         try:
#             core.clean.cleanx(options.complexName[1])
#         except FileNotFoundError:
#             print(f'Nothing to clean inside {options.complexName[1]}')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='minnie', description=__doc__)
    parser.add_argument(
        '--nproc',
        default=None,
        type=int,
        help='Number of processes to use for parallel execution.'
    )

    # Setup subparsers
    subparsers = parser.add_subparsers(dest="subcommand")

    # minnie splitpdb
    splitpdb = subparsers.add_parser(
        'splitpdbs',
        help='Split a trajectory into single frames'
    )
    splitpdb.add_argument(
        'pdbs',
        nargs='+',
        type=pathlib.Path,
        help='Input trajectory file(s) in PDB format'
    )
    splitpdb.add_argument(
        '-i'  # simpler than -cn
        '--id',  # simpler than --complexName
        nargs='+',
        type=str,
        dest='project_ids',
        help='Project ID for your analyses. One per trajectory.',
    )

    # minnie findbonds
    findbonds = subparsers.add_parser(
        'findbonds',
        help='Calculates interactions between and/or within monomers'
    )
    findbonds.add_argument(
        '-i'
        '--id',
        nargs='+',
        type=str,
        dest='project_id',
        help='Project ID for your analysis.',
    )
    inputopts = findbonds.add_mutually_exclusive_group(required=True)
    inputopts.add_argument(
        '-f',
        '--pdbfile',
        type=pathlib.Path,
        help='Input PDB file in PDB format.'
    )
    inputopts.add_argument(
        '-d',
        '--folder',
        type=pathlib.Path,
        help='Input directory with PDB files.'
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
        help='Interaction types to analyze.',
    )
    findbonds.add_argument(
        '--intra',
        action='store_true',
        help='Include intra-monomer interactions.',
    )
    findbonds.add_argument(
        '--clean',
        action='store_true',
        help='Remove intermediate files upon completion.',
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
