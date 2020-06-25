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


class guivar:
    splitpdbs = ['%(prog)s\n'
                 f'                        -i, --id     <string> <string>     \n '
                 f'                                     Project ID for your analyses. One per trajectory.\n\n'

                 f'                        --pdbs       [<traj.pdb>] [<traj.pdb>]   \n'
                 f'                                      Trajectory ID of your complex(s)\033[0m \n\n\n\n'

                 f'\n\033[1mUsage example:\033[0m\n\n'
                 " minnie splitpdbs -i sox4 sox18 --pdbs sox4.pdb sox18.pdb \n"
                 " minnie splitpdbs -i sox4  --pdbs sox4.pdb  \n"]

    findbonds = ['%(prog)s\n'
                 f'                        -i, --id      <string>     \n '
                 f'                                      Project ID for your analysis.\n\n'

                 f'                        -f, --pdbfile [<.pdb>] (singleframe.pdb)   \n'
                 f'                                      Input PDB file in PDB format. \n\n'

                 f'                        -d, --folder  [<path>]                  \n'
                 f'                                      Input directory with PDB files.\n\n'

                 f'                        --itypes      [<hbonds>/<ionic>/<hydrophobic>/<ring_stacking>/<all>] (hbonds)'
                 f'\n                                      Calculates which types of interactions \n\n'

                 f'                        --intra        [<"True">/<"False">] ("False")  \n'
                 f'                                      Include intra-monomer interactions. \033[0m \n\n\n\n'

                 f'                        --clean\n'
                 f'                                      Remove intermediate files upon completion. \033[0m \n'
                 f'\n\n\n'


                 f'\n\033[1mUsage example:\033[0m\n\n'
                 " Single frame      - minnie findbonds -i sox4 -f 'sox4/01_frames/md_0.pdb' --itypes all \n"
                 " Multiple frames   - minnie findbonds -i sox4 -d 'sox4/01_frames/' --clean \n"
                 " Multiple frames   - minnie findbonds -i sox4 -d 'sox4/01_frames/' --intra \n"]
    timefilter = ['%(prog)s\n'
                  f'                        -i, --id     <string> <string>     \n '
                  f'                                      sProject ID of your complex\n\n'

                  f'                         -f, --files            [<.csv>]               \n'
                  f'                                                Files of your complex\n\n'
                  f'                         --per                  <float>                \n'
                  f'                                                Observation frequency to classify an interaction as critical\033[0m \n\n\n\n'

                  f'\n\033[1mUsage example:\033[0m\n\n'
                  " Single file    - minnie timefilter -f sox4/02_interfacea_results/hbonds/sox4_merged_hbonds.csv -i sox4 --per 25  \n"
                  " Multiple files - minnie timefilter -f sox4/02_interfacea_results/*/sox4_merged_*.csv -i sox4 --per 25  \n"]
    comparecx = ['%(prog)s\n'
                 f'                        -i, --id     <string> <string>     \n '
                 f'                                      sProject ID of your complex\n\n'


                 f'                         --per       <float>                \n'
                 f'                                     Observation frequency to classify an interaction as critical\033[0m \n\n\n\n'

                 f'\n\033[1mUsage example:\033[0m\n\n'
                 " minnie comparecx -i sox4 sox18 --per 25  \n"]

    graph = ['%(prog)s\n'
             f'                        -i, --id               <string> <string>     \n'
             f'                                               Project IDs of your complex(s)\n\n'

             f'                        --per                  <float>                \n'
             f'                                               Observation frequency to classify an interaction as critical \n\n'


             f'                        -b, --between          [<protein-dna>/<all>]   \n'
             f'                                               Between protein-dna or keep all \n\n'

             f'                        -c, --chainIDs         <string> <string>       \n'
             f'                                               Give ChainIDs to proceed\n\n'

             f'                        --filename             <string>                           \n'
             f'                                               Give a name to output file (optional)\n\n'

             f'                        --colors               [<hex colors>] ("#D9B4CC", "#6F81A6")     \n'
             f'                                               Color IDs of the complexes (optional)\n\n'


             f'                        --itypes               [<hbonds>/<ionic>/<hydrophobic>/<ring_stacking>/<all>] (hbonds)    \n'
             f'                                               Calculates which types of interactions \n\n'

             f'                        -s                     [<specific>/<common>] (specific)                           \n'
             f'                                               Complex-specific or common interactions\033[0m \n\n\n\n'

             f'Please do not use "--between" and "--chainIDs" options at the same time\n\n'

             #"\n\033[1m Usage example: \033[0m\n\n"
             "\nUsage example:  \n\n"
             " minnie graph -i 'sox4' 'sox18' --per 25 --itypes hbonds -s specific  -c A+B C --colors '#D9B4CC' '#6F81A6' \n"
             " minnie graph -i 'sox4' 'sox18' --per 25 --itypes ionic -c A+B C  \n"
             " minnie graph -i 'sox4' 'sox18' --per 25 --itypes ionic -b protein-dna \n"
             " minnie graph -i 'sox4' 'sox18' --per 25 --itypes ionic -b protein-dna --filename sox4_sox18_protein_dna \n"]
