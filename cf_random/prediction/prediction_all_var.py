#!/bin/env python3
# -*- coding: utf-8 -*-
"""Helpers to run ColabFold prediction batches across varied MSAs.

Small utilities used by higher-level workflows to orchestrate
multiple ColabFold runs with different MSA settings.
"""

import logging
import os
import random
from pathlib import (
    Path,
)

import numpy as np
from colabfold.batch import (
    get_queries,
    run,
)
from colabfold.utils import (
    setup_logging,
)


class CF_MSA_MAX:
    def __init__(self, search_dir, output_dir, pdb_name, rseed, num_seeds, model_type) -> str:
        print(search_dir)

        setup_logging(Path(output_dir) / "log.txt")
        logger = logging.getLogger(__name__)

        queries, is_complex = get_queries(search_dir)

        run(
            queries=queries,
            result_dir=output_dir,
            num_models=5,
            is_complex=is_complex,
            model_type=model_type,
            num_seeds=int(num_seeds),
            random_seed=int(rseed),
            data_dir=Path("."),
        )


class CF_MSA_VAR:
    def __init__(self, pdb1_name, search_dir, output_dir, rseed, num_seeds, model_type):
        #### shallow MSA section
        #### Global viarlable
        max_msa = 1
        ext_msa = 2
        pre_random_seed = np.array(rseed)  ## needed to remove future
        random_seed = "".join(map(str, pre_random_seed))

        self.pdb1_name = pdb1_name

        max_msa = 1
        ext_msa = 2

        TMscores_random = []

        for multi in (1, 2, 2, 2, 2, 2, 2):
            max_msa = max_msa * multi
            ext_msa = ext_msa * multi

            output_dir_var = (
                output_dir + str(random_seed) + "_max_" + str(max_msa) + "_ext_" + str(ext_msa)
            )

            setup_logging(Path(output_dir_var) / "log.txt")
            logger = logging.getLogger(__name__)

            queries, is_complex = get_queries(search_dir)

            run(
                queries=queries,
                result_dir=output_dir_var,
                num_models=5,
                is_complex=is_complex,
                model_type=model_type,
                num_seeds=int(num_seeds),
                random_seed=int(random_seed),
                max_seq=int(max_msa),
                max_extra_seq=int(ext_msa),
                data_dir=Path("."),
            )

        fin_pred_dir = pdb1_name + "_predicted_models_rand_" + str(random_seed) + "_max_*"
        gen_dir = "predictions_all/" + pdb1_name

        if not os.path.exists(gen_dir):
            os.makedirs(gen_dir)
            mv_command = "mv " + fin_pred_dir + " predictions_all/" + pdb1_name
            print(mv_command)
            os.system(mv_command)
        else:
            mv_command = "mv " + fin_pred_dir + " predictions_all/" + pdb1_name
            print(mv_command)
            os.system(mv_command)


class PredictionAll:
    def __init__(self, pdb1_name, search_dir, search_multi_dir, nMSA, model_type):

        num_seeds = 5 + nMSA

        pre_random_seed = np.random.randint(0, 16, 1)
        random_seed = "".join(map(str, pre_random_seed))
        output_dir = pdb1_name + "_predicted_models_full_rand_" + str(random_seed)

        ##### Perform predction with full-length MSA
        MSA_full = CF_MSA_MAX(search_dir, output_dir, pdb1_name, random_seed, num_seeds, model_type)
        pwd = os.getcwd() + "/"

        # Directory section
        gen_dir = "predictions_all/" + pdb1_name

        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)

        pred_dir = pdb1_name + "_predicted_models_full_rand_" + str(random_seed) + "/"
        mv_folder_cmd = "mv " + pred_dir + " predictions_all/" + pdb1_name
        print(mv_folder_cmd)
        os.system(mv_folder_cmd)

        ##### check out varied-MSA with (msa-max: 1, 2, 4, 8, 16, 32, 64) (msa-extra: 2, 4, 8, 16, 32, 64, 128)
        output_dir = pdb1_name + "_predicted_models_rand_"
        random_seed = random.sample(range(100), 1)
        if model_type != "alphafold2_multimer_v3":
            MSA_var = CF_MSA_VAR(
                pdb1_name, search_dir, output_dir, random_seed, num_seeds, model_type
            )
        else:
            MSA_var = CF_MSA_VAR(
                pdb1_name, search_multi_dir, output_dir, random_seed, num_seeds, model_type
            )
