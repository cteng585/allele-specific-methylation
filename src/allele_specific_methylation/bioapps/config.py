import json
import re
import warnings
from pathlib import Path
from typing import Literal

import yaml

from allele_specific_methylation.bioapps.classes import RequestHandler


def patient_identifier_mappings(
    request_handle: RequestHandler,
    participant_study_identifiers: list[str],
) -> list[str]:
    """POG IDs are not the primary identifiers in the BioApps database

    Multiple different "participant_study_identifiers" can map to the same "patient_identifer".
    Use the "patient_analysis" endpoint to get the "patient_identifier" values for a list of
    POG IDs

    :return:
    """
    patient_analysis_params = {
        "participant_study_identifier": ",".join(
            re.search(r"^POG[0-9]*", study_id).group(0) for study_id in participant_study_identifiers
        ),
    }
    response = request_handle.get(
        "patient_analysis",
        params=patient_analysis_params,
    )
    patient_identifiers = [record["patient_identifier"] for record in response.json()]

    patient_params = {
        "patient_identifier": ",".join(patient_identifiers),
        "relations": "sources",
        "sources_columns": "participant_study_identifier",
    }
    response = request_handle.get(
        "patient",
        params=patient_params,
    )

    # patient identifier is the primary ID in the BioApps db for a given patient
    # participant_study_identifier is ordinarily the POG ID but may also be a PATH ID
    patient_mappings = {}
    for record in response.json():
        patient_identifier = record["patient_identifier"]
        patient_study_ids = [
            source.get("participant_study_identifier")
            for source in record.get("sources", [])
            if
            source.get("participant_study_identifier") and re.search(
                r"POG", source.get("participant_study_identifier")
            )
        ]

        for patient_study_id in patient_study_ids:
            patient_mappings[patient_identifier] = patient_study_id

    return patient_mappings


def bioapps_libraries(
    request_handle: RequestHandler,
    patient_mappings: dict[str, str],
):
    library_params = {
        "patient": ",".join(patient_mappings.keys())
    }
    response = request_handle.get(
        "library/info",
        params=library_params,
    )

    source_libraries = {}
    for record in response.json():
        library_name = record.get("name", None)
        if not library_name:
            continue

        source = record.get("source", None)
        if not source:
            continue

        if source.get("pathology"):
            pathology = source["pathology"]
        elif source.get("pathology_type"):
            pathology = source["pathology_type"]
        else:
            msg = f"No pathology information found for library {library_name}"
            warnings.warn(msg, RuntimeWarning)
            continue

        source_name = source["original_source_name"]
        patient_identifier = source["patient"]["patient_identifier"]

        match pathology:
            case "Normal":
                library_label = "NORMAL"
            case "Diseased" | "Malignant":
                library_label = "TUMOR"
            case _:
                msg = f"Unknown or unhandled pathology {pathology}"
                warnings.warn(msg, RuntimeWarning)
                continue

        if source_name not in source_libraries:
            source_libraries[source_name] = {
                "patient_id": patient_identifier,
                "libraries": {
                    library_name: library_label
                },
            }
        elif patient_identifier != source_libraries[source_name]["patient_id"]:
            msg = (
                f"Source {source_name} has multiple patient IDs: "
                f"{patient_identifier}, {source_libraries[source_name]['patient_id']}"
            )
            raise KeyError(msg)
        else:
            source_libraries[source_name]["libraries"][library_name] = library_label

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
    bioapps_id_to_pog = patient_identifier_mappings(request_handle, source_directories)
    source_libraries = bioapps_libraries(request_handle, bioapps_id_to_pog)

    config = {}
    for source_directory in source_directories:
        patient_id = source_libraries[source_directory]["patient_id"]
        study_id = bioapps_id_to_pog[patient_id]
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
