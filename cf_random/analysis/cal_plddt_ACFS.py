# 0_000_scores_rank_001_alphafold2_ptm_model_4_seed_000.json
pattrn = re.compile(r".*?_scores_rank_(?P<rank>\d+)_alphafold2.*")

# default if pattern doesn't work
rank = "000"


def read_plddt(jsonfile):
    """Reads pLDDT scores from an AlphaFold prediction JSON file.

    Args:
        jsonfile (str): Path to the JSON file containing prediction scores.

    Returns:
        numpy.ndarray: Array of pLDDT scores as float64.
    """
    import json
    import numpy as np

    with open(jsonfile) as json_file:
        data = json.load(json_file)

    plddt_scores = np.array(data["plddt"], dtype="float64")

    return plddt_scores


def fract_good(score):
    """Calculates the average pLDDT score from an array of scores.

    Args:
        score (numpy.ndarray): Array of pLDDT scores.

    Returns:
        float: Average pLDDT score rounded to 2 decimal places.
    """
    import numpy as np

    avg_plddt = round(np.average(score), 2)
    return avg_plddt


class plddt_cal:
    """Calculates average pLDDT scores for protein models across different categories.

    This class processes AlphaFold prediction JSON files to extract pLDDT scores
    and computes average scores for different MSA categories and model types.
    """

    def __init__(self, sub_list, category, pdb_name, nMSA, nENS, model_type):
        """Initializes pLDDT calculation for given subdirectories and parameters.

        Args:
            sub_list (list): List of subdirectory paths to process.
            category (str): MSA category ('full-MSA', 'additional-MSA', 'random-MSA').
            pdb_name (str): Name of the PDB structure.
            nMSA (int): Number of MSA sequences.
            nENS (int): Number of ensemble models.
            model_type (str): Type of AlphaFold model.
        """
        import sys
        import numpy as np

        if len(sub_list) == 0:
            sys.exit(1)

        print("working...")
        print(sub_list)

        values_all, out_dict_all, cnt = self._process_subdirs(sub_list)

        if category == "full-MSA":
            cnt = int(cnt / 5)
        # For additional-MSA and others, cnt remains as is

        print(cnt)
        print(values_all)

        # Reshape based on category and model_type
        if category == "full-MSA":
            values_all_resh = values_all.reshape(nMSA + 5, 5)
        elif category == "additional-MSA" and model_type == "alphafold2_multimer_v3":
            values_all_resh = values_all.reshape((nENS + 20), 5)
        elif category == "additional-MSA":
            values_all_resh = values_all.reshape((nENS + 20), 5)
        elif category == "random-MSA" and model_type != "alphafold2_multimer_v3":
            values_all_resh = values_all.reshape((nMSA + 5) * 7, 5)
        elif category == "random-MSA":
            values_all_resh = values_all.reshape((nMSA + 5) * 7, 5)

        print("                ")
        print("Calculated pLDDT")
        print(values_all_resh)
        np.savetxt("plddt_" + category + "_" + pdb_name + ".csv", values_all_resh, fmt="%2.3f")

    def _process_subdirs(self, sub_list):
        """Processes subdirectories to extract pLDDT scores.

        Args:
            sub_list (list): List of subdirectory paths.

        Returns:
            tuple: (values_all, out_dict_all, cnt) where values_all is numpy array
                   of scores, out_dict_all is dict of key-value pairs, cnt is count.
        """
        import glob
        from pathlib import Path
        import numpy as np

        out_dict_all = {}
        values_all = np.array([])
        cnt = 0

        for subdir in sub_list:
            if Path(subdir).is_dir():
                subdir_name = Path(subdir).name
                jsonfiles = glob.glob(str(subdir) + "/*_scores*json")

                for jsonfile in jsonfiles:
                    plddt_score = read_plddt(jsonfile)
                    values = fract_good(plddt_score)
                    values_all = np.append(values_all, values)
                    jsonfilepath = Path(jsonfile)
                    jsonfilename = jsonfilepath.stem
                    match = pattrn.match(jsonfilename)
                    if match:
                        rank = match.group("rank")

                    key_pair = subdir_name + ":" + rank
                    if key_pair not in out_dict_all:
                        out_dict_all[key_pair] = values

                    cnt += 1

        return values_all, out_dict_all, cnt
