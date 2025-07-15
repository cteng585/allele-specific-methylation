import json
import re
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
        "participant_study_identifier": ",".join(participant_study_identifiers),
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


def patient_libraries(
    request_handle: RequestHandler,
    participant_study_identifiers: list[str],
    patient_mappings: dict[str, str],
):
    library_params = {
        "patient": ",".join(patient_mappings.keys())
    }
    response = request_handle.get(
        "library/info",
        params=library_params,
    )

    participant_libraries = {
        participant_study_id: [] for participant_study_id in participant_study_identifiers
    }
    for record in response.json():
        library_id = record.get("name", None)
        if not library_id:
            continue

        source = record.get("source")
        if source and source.get("pathology"):
            patient_identifier = source["patient"]["patient_identifier"]
            library_label = record["source"]["pathology"]

            match library_label:
                case "Normal":
                    library_label = "NORMAL"
                case "Diseased":
                    library_label = "TUMOR"
                case _:
                    raise Exception(f"Unknown library label: {library_label}")

            protocol_name = record["protocol"]["name"]
            if re.search(r"nanopore|genome", protocol_name, re.IGNORECASE):
                participant_libraries[patient_mappings[patient_identifier]].append(
                    (library_id, library_label)
                )

    # check that at least 3 libraries (one for short read normal, one for short read tumor, one for long read tumor) are present for each participant
    wrong_library_num = []
    for study_identifier in participant_libraries:
        if len(participant_libraries[study_identifier]) < 3:
            wrong_library_num.append(study_identifier)

    if wrong_library_num:
        raise ValueError(
            f"Not enough libraries found for the following POG IDs: {', '.join(wrong_library_num)}. "
            "At least one short read normal, one short read tumor, and one long read tumor library are required."
        )

    return participant_libraries


def generate_config(
    analysis_dir: str | Path,
    config_type: Literal["yaml", "json"] = "yaml",
    config_path: str | Path = "config.yaml",
):
    analysis_dir = Path(analysis_dir.removesuffix("/"))
    participant_study_identifiers = [
        path_obj.stem for path_obj in Path(analysis_dir).iterdir()
    ]
    request_handle = RequestHandler()
    bioapps_id_to_pog = patient_identifier_mappings(request_handle, participant_study_identifiers)
    pog_libraries = patient_libraries(request_handle, participant_study_identifiers, bioapps_id_to_pog)

    config = {}
    for study_identifier in pog_libraries:
        libraries = {
            library[0]: library[1] for library in pog_libraries[study_identifier]
        }

        # account for cases where a study identifier is a constant created by the clair variant calling workflow
        libraries["SAMPLE"] = "TUMOR"

        config[study_identifier] = {
            "short_read_snv": {
                "path": f"{analysis_dir}/{study_identifier}/short_read.snvs.vcf",
                "rename": True,
                "libraries": libraries,
            },
            "short_read_indel": {
                "path": f"{analysis_dir}/{study_identifier}/short_read.indels.vcf",
                "rename": True,
                "libraries": libraries,
            },
            "long_read": {
                "path": f"{analysis_dir}/{study_identifier}/long_read.vcf.gz",
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

        case _:
            raise NotImplementedError(f"Unknown config type: {config_type}")
