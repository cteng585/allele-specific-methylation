import json
import re
import warnings
from pathlib import Path
from typing import Literal

import yaml

from allele_specific_methylation.bioapps.classes import RequestHandler


def bioapps_libraries(
    request_handle: RequestHandler,
    sample_ids: list[str],
) -> list[str]:
    """POG IDs are not the primary identifiers in the BioApps database

    Multiple different "participant_study_identifiers" can map to the same "patient_identifer".
    Use the "patient_analysis" endpoint to get the "patient_identifier" values for a list of
    POG IDs

    :return:
    """
    patient_analysis_params = {
        "participant_study_identifier": ",".join(
            re.search(r"^POG[0-9]*", study_id).group(0) for study_id in sample_ids
        ),
    }
    response = request_handle.get(
        "patient_analysis",
        params=patient_analysis_params,
    )

    if not response.ok:
        msg = (
            f"Failed to retrieve patient analysis information from BioApps API: "
            f"{response.status_code} {response.reason}"
        )
        raise ConnectionError(msg)

    source_libraries = {}
    for record in response.json():
        patient_id = record.get("patient_identifier", None)
        patient_libraries = {}
        study_id = set()
        patient_sources = []
        for source in record["sources"]:
            patient_sources.append(source["original_source_name"])
            study_id.add(source["participant_study_identifier"])
            pathology = source["pathology"]

            match pathology:
                case "Normal":
                    library_label = "NORMAL"
                case "Diseased" | "Malignant":
                    library_label = "TUMOR"
                case _:
                    msg = f"Unknown or unhandled pathology {pathology}"
                    warnings.warn(msg, RuntimeWarning)
                    continue

            for library in source["libraries"]:
                patient_libraries[library["name"]] = library_label

        if len(study_id) > 1:
            msg = (
                f"Multiple study identifiers found for patient {patient_id}: "
                f"{', '.join(study_id)}"
                f"Using the first one"
            )
            warnings.warn(msg, RuntimeWarning)

        study_id = study_id.pop()
        for source_name in patient_sources:
            source_libraries[source_name] = {
                "patient_id": patient_id,
                "study_id": study_id,
                "libraries": patient_libraries,
            }

    return source_libraries


def generate_config(
    analysis_dir: str | Path,
    config_type: Literal["yaml", "json"] = "yaml",
    config_path: str | Path = "config.yaml",
):
    analysis_dir = Path(analysis_dir.removesuffix("/"))
    source_directories = [
        path_obj.stem for path_obj in Path(analysis_dir).iterdir()
    ]
    request_handle = RequestHandler()
    source_libraries = bioapps_libraries(request_handle, source_directories)

    config = {}
    for source_directory in source_directories:
        patient_id = source_libraries[source_directory]["patient_id"]
        study_id = source_libraries[source_directory]["study_id"]
        libraries = source_libraries[source_directory]["libraries"]

        # account for cases where a study identifier is a constant created by the clair variant calling workflow
        libraries["SAMPLE"] = "TUMOR"

        # TODO: might want to change the structure of this to:
        # config[source_directory] = {
        #     "patient_id": patient_id,
        #     "study_id": study_id,
        #     "libraries": {...},
        #     "paths": {
        #         "short_read_snv": "...",
        #         "short_read_indel": "...",
        #         "long_read": "...",
        #     },
        # }
        config[source_directory] = {
            "patient_id": patient_id,
            "study_id": study_id,
            "short_read_snv": {
                "path": f"{analysis_dir}/{source_directory}/short_read.snvs.vcf",
                "rename": True,
                "libraries": libraries,
            },
            "short_read_indel": {
                "path": f"{analysis_dir}/{source_directory}/short_read.indels.vcf",
                "rename": True,
                "libraries": libraries,
            },
            "long_read": {
                "path": f"{analysis_dir}/{source_directory}/long_read.vcf.gz",
                "rename": True,
                "libraries": libraries,
            }
        }

    match config_type:
        case "yaml":
            yaml.Dumper.ignore_aliases = lambda *args: True
            with open(config_path, "w") as outfile:
                yaml.dump(config, outfile, default_flow_style=False)

        case "json":
            with open(config_path, "w") as outfile:
                json.dump(config, outfile, indent=4)

    return
