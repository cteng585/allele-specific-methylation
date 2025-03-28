from typing import Literal, Optional, Union

from pydantic import BaseModel, ConfigDict


class VCFMetadata(BaseModel):
    """
    TODO: add docstring
    """

    ID: str
    Number: Union[int, Literal["A", "R", "G", "."]]
    Description: str

    model_config = ConfigDict(extra="allow")


class VCFContigField(VCFMetadata):
    """
    Represents a contig field in a VCF file
    """

    Number: None = None
    Description: None = None
    MetadataType: str = "ContigField"
    length: int
    assembly: Optional[str]


class VCFFormatField(VCFMetadata):
    """
    Represents a FORMAT field in a VCF file
    """

    MetadataType: str = "FormatField"
    Type: Literal["Integer", "Float", "Character", "String"]


class VCFInfoField(VCFMetadata):
    """
    Represents an INFO field in a VCF file
    """

    MetadataType: str = "InfoField"
    Type: Literal["Integer", "Float", "Flag", "Character", "String"]
    Source: Optional[str] = None
    Version: Optional[str] = None
