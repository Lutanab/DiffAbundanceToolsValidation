from enum import Enum
from pydantic import BaseModel
from pydantic import Field
from pydantic import model_validator
from pydantic import validator
from typing import Text
from typing import List
from typing import Union
from typing_extensions import Annotated
from typing_extensions import Self


Proportion = Annotated[float, Field(ge=0, le=1)]


class RBool(Enum):
    TRUE = "TRUE"
    FALSE = "FALSE"

    def __bool__(self):
        return self == RBool.TRUE


class CountTransform(Enum):
    CLR = "clr"
    LOG10P = "log10p"
    LOG10 = "log10"
    COMPOSITIONAL = "compositional"


class MicrobiomeDistance(Enum):
    BRAY = "bray"
    JSD = "jsd"


class RawGlmParams(BaseModel):
    lfc_range: float
    p_val_treshold: Proportion
    transforms: List[CountTransform]


class ManualRawGlmParams(BaseModel):
    lfc_range: float
    transforms: List[CountTransform]


class PcaParams(BaseModel):
    transform: CountTransform
    distance: MicrobiomeDistance


class MetricsParams(BaseModel):
    lfc_range: float
    p_val_treshold: Proportion


class GenerationParams(BaseModel):
    da_taxa_percentage: Proportion
    prevalence_upper_limit: Proportion
    prevalence_lower_limit: Proportion
    n_samples: int
    tar_samples_percentage: Proportion
    mannual_params: RBool

    def taxa_params_eq(self: 'GenerationParams', other: 'GenerationParams'):
        if not isinstance(other, GenerationParams):
            return False
        return (self.da_taxa_percentage == other.da_taxa_percentage) and\
            (self.prevalence_upper_limit == other.prevalence_upper_limit) and\
            (self.prevalence_lower_limit == other.prevalence_lower_limit)
            

    @validator("n_samples")
    def check_positive_number(cls, value):
        if value < 0:
            raise ValueError("n_samples must be positive integer")
        return value

    @model_validator(mode='after')
    def validate_upper_lower_limits(self):
        if self.prevalence_lower_limit > self.prevalence_lower_limit:
            raise ValueError("prevalence_upper_limit must be >= prevalence_lower_limit")
        return self


class SnakemakeConfig(BaseModel):
    series_number: int
    series_name: Text
    fixed_da_taxa: RBool
    data_generator: Text
    source: Text
    generation_params_list: List[GenerationParams]
    raw_glm_params_list: List[RawGlmParams]
    manual_raw_glm_params_list: List[ManualRawGlmParams]
    pca_params_list: List[PcaParams]
    metrics_params_list: List[MetricsParams]

    @model_validator(mode='after')
    def validate_taxa_params_consistency(self):
        if self.fixed_da_taxa == RBool.TRUE and len(self.generation_params_list) > 1:
            generation_params_1 = self.generation_params_list[0]
            for gen_params in self.generation_params_list[1:]:
                if not generation_params_1.taxa_params_eq(gen_params):
                    raise ValueError("DA taxa generation params inconsistent")
        return self


