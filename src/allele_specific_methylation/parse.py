import json
import yaml
from pathlib import Path
from typing import Literal


def parse_combine_vcf_config(config_fn: str | Path, file_type: Literal["yaml", "json"]):
    """Parse a configuration file for combining VCFs for multiple samples

    :param config_fn: path to the configuration file
    :param file_type: type of the configuration file, either "yaml" or "json"
    :return: the parsed configuration as a dictionary
    """
    match file_type:
        case "yaml":
            with open(config_fn, "r") as yaml_file:
                loaded_config = yaml.safe_load(yaml_file)

        case "json":
            with open(config_fn, "r") as json_file:
                loaded_config = json.load(json_file)

        case _:
            raise NotImplementedError(
                f"Handling for config of type {file_type} is not implemented"
            )

    sample_configs = {}
    for sample_id in loaded_config:
        config_dict = {
            "short_read_snv": None,
            "short_read_indel": None,
            "long_read": None,
        }
        if any(vcf_type not in loaded_config[sample_id] for vcf_type in config_dict):
            raise ValueError(
                f"Config file must contain 'short_read_snv', 'short_read_indel', and 'long_read' VCFs for {sample_id}"
            )

        for vcf_type in config_dict:
            if "path" in loaded_config[sample_id][vcf_type]:
                config_dict[vcf_type] = {
                    "path": loaded_config[sample_id][vcf_type]["path"]
                }
            else:
                raise ValueError(
                    f"A path must be provided for {vcf_type} in "
                )

            if loaded_config[sample_id][vcf_type]["rename"] is True:
                rename_dict = {}
                if "normal_library_id" in loaded_config[sample_id][vcf_type]:
                    library_id = loaded_config[sample_id][vcf_type]["normal_library_id"]
                    rename_dict[library_id] = "NORMAL"
                if "tumor_library_id" in loaded_config[sample_id][vcf_type]:
                    library_id = loaded_config[sample_id][vcf_type]["tumor_library_id"]
                    rename_dict[library_id] = "TUMOR"

                config_dict[vcf_type]["rename"] = rename_dict
            else:
                config_dict[vcf_type]["rename"] = None

        sample_configs[sample_id] = config_dict

    return sample_configs
