import re
from typing import Literal, Union
from typing_extensions import Annotated

from pydantic import BaseModel, ConfigDict


class VCFFormatField(BaseModel):
    """
    ##FORMAT=<ID=ID,Number=number,Type=type,Description="description">
    """
    ID: str
    Number: Union[int, Literal["A", "R", "G", "."]]
    Type: Literal["Integer", "Float", "Character", "String"]
    Description: str
    
    model_config = ConfigDict(extra='allow')
