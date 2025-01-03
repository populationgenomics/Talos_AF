"""
A home for all data models used in Talos
"""

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field

NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
CHROM_ORDER = list(map(str, range(1, 23))) + NON_HOM_CHROM

# some kind of version tracking
CURRENT_VERSION = '0.0.1'

# ratios for use in AB testing
MAX_WT = 0.15
MIN_HET = 0.25
MAX_HET = 0.75
MIN_HOM = 0.85


class FileTypes(Enum):
    """
    enumeration of permitted input file types
    """

    HAIL_TABLE = '.ht'
    MATRIX_TABLE = '.mt'
    PED = 'ped'
    VCF = '.vcf'
    VCF_GZ = '.vcf.gz'
    VCF_BGZ = '.vcf.bgz'


class PhenoPacketHpo(BaseModel):
    """
    A representation of a HPO term
    """

    id: str
    label: str


class Coordinates(BaseModel):
    """
    A representation of genomic coordinates
    """

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def string_format(self) -> str:
        """
        forms a string representation: chr-pos-ref-alt
        """
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def __lt__(self, other) -> bool:
        """
        enables positional sorting
        """
        # this will return False for same chrom and position
        if self.chrom == other.chrom:
            return self.pos < other.pos
        # otherwise take the relative index from sorted chromosomes list
        if self.chrom in CHROM_ORDER and other.chrom in CHROM_ORDER:
            return CHROM_ORDER.index(self.chrom) < CHROM_ORDER.index(other.chrom)
        # if self is on a canonical chromosome, sort before HLA/Decoy etc.
        if self.chrom in CHROM_ORDER:
            return True
        return False


class VariantCommon(BaseModel):
    """
    the abstracted representation of a variant from any source
    """

    coordinates: Coordinates = Field(repr=True)
    info: dict[str, str | int | float | list[str] | list[float] | dict[str, str] | bool] = Field(default_factory=dict)
    het_samples: set[str] = Field(default_factory=set, exclude=True)
    hom_samples: set[str] = Field(default_factory=set, exclude=True)
    sample_support: list[str] = Field(default_factory=list, exclude=True)
    phased: dict = Field(default_factory=dict)

    def __str__(self):
        return repr(self)

    def __lt__(self, other):
        return self.coordinates < other.coordinates

    def __eq__(self, other):
        return self.coordinates == other.coordinates

    def check_ab_ratio(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for AB ratio checking - not implemented for SVs
        """
        return set()

    def get_sample_flags(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for flag checking - not implemented for SVs (yet)
        """
        return set()

    def check_read_depth(self, *args, **kwargs) -> set[str]:  # noqa: ARG002, ANN002, ANN003
        """
        dummy method for read depth checking - not implemented for SVs
        """
        return set()


class SmallVariant(VariantCommon):
    """
    a representation of a small variant
    note that transcript_consequences is not optional
    we require that something specific to SmallVariant(s) is mandatory
    this is in order to correctly deserialise the Small/Structural objects
    into the appropriate types. If it is optional, Pydantic can coerce
    everything as a SmallVariant
    """

    depths: dict[str, int] = Field(default_factory=dict, exclude=True)
    ab_ratios: dict[str, float] = Field(default_factory=dict, exclude=True)
    transcript_consequences: list[dict[str, str | float | int]]

    def get_sample_flags(self, sample: str) -> set[str]:
        """
        gets all report flags for this sample - currently only one flag
        """
        return self.check_ab_ratio(sample) | self.check_read_depth(sample)

    def check_read_depth(self, sample: str, threshold: int = 10, var_is_cat_1: bool = False) -> set[str]:
        """
        flag low read depth for this sample

        Args:
            sample (str): sample ID matching VCF
            threshold (int): cut-off for flagging as a failure
            var_is_cat_1 (bool): flag if this variant is a category 1

        Returns:
            return a flag if this sample has low read depth
        """
        if var_is_cat_1:
            return set()
        if self.depths[sample] < threshold:
            return {'Low Read Depth'}
        return set()

    def check_ab_ratio(self, sample: str) -> set[str]:
        """
        AB ratio test for this sample's variant call

        Args:
            sample (str): sample ID

        Returns:
            set[str]: empty, or indicating an AB ratio failure
        """
        het = sample in self.het_samples
        hom = sample in self.hom_samples
        variant_ab = self.ab_ratios.get(sample, 0.0)

        if (variant_ab <= MAX_WT) or (het and not MIN_HET <= variant_ab <= MAX_HET) or (hom and variant_ab <= MIN_HOM):
            return {'AB Ratio'}
        return set()


class AFCategoryC(BaseModel):
    """
    details relevant to a category C gene
    """

    # a list of specific individual variants
    only_variants: Optional[list[str]] = Field(default=None)
    # a list of specific individual protein consequences
    only_protein_consequences: Optional[list[str]] = Field(default=None)
    # a list of specific individual disorders, only ClinVar entries with these annotated disorders are relevant
    clinvar_disorders: Optional[list[str]] = Field(default=None)
    # only report on high impact/clinvar path on these transcripts
    only_transcripts: Optional[list[str]] = Field(default=None)


class AFGeneEntity(BaseModel):
    a: bool = Field(default=False)
    b: bool = Field(default=False)
    c: Optional[AFCategoryC] = Field(default=None)


class AFSpecification(BaseModel):
    """
    A representation of the AF specification
    """

    genes: dict[str, str] = Field(default_factory=dict)
    metadata: dict[str, str] = Field(default_factory=dict)


class MinimalVariant(BaseModel):
    """
    A minimal representation of a variant for reporting purposes
    """

    gene: str
    coordinates: Coordinates = Field(repr=True)
    info: dict[str, Any] = Field(default_factory=dict)
    genotype: str = Field(default_factory=str)


class ReportableVariant(BaseModel):
    """
    A variant passing MOI tests, to be reported
    2 components - the mandatory first variant, and optionally the second variant
    """

    sample: str
    var_1: MinimalVariant
    var_2: Optional[MinimalVariant] = None
    gene: str
    rule: str
    flags: set[str] = Field(default_factory=set)
    labels: set[str] = Field(default_factory=set)

    def __eq__(self, other):
        """
        makes reported variants comparable
        """
        return self.sample == other.sample and self.var_1.coordinates == other.var_1.coordinates

    def __lt__(self, other):
        return self.var_1.coordinates < other.var_1.coordinates


class ResultMeta(BaseModel):
    """
    metadata for a result set
    """

    version: str = Field(default_factory=str)
    input_file: str = Field(default_factory=str)


class MemberSex(Enum):
    MALE = 'male'
    FEMALE = 'female'
    UNKNOWN = 'unknown'


class FamilyMembers(BaseModel):
    affected: bool = Field(default=False)
    ext_id: str = Field(default_factory=str)
    sex: MemberSex = Field(default=MemberSex.UNKNOWN)


class ParticipantMeta(BaseModel):
    ext_id: str
    family_id: str
    members: dict[str, FamilyMembers] = Field(default_factory=dict)
    phenotypes: list[PhenoPacketHpo] = Field(default_factory=list)


class ParticipantResults(BaseModel):
    """
    A representation of a result set
    """

    variants: list[ReportableVariant] = Field(default_factory=list)
    metadata: ParticipantMeta = Field(default_factory=ParticipantMeta)


class ResultData(BaseModel):
    """
    A representation of a result set
    """

    results: dict[str, ParticipantResults] = Field(default_factory=dict)
    metadata: ResultMeta = Field(default_factory=ResultMeta)
    version: str = CURRENT_VERSION


class PedigreeMember(BaseModel):
    """
    This will be a more searchable implementation of the peds pedigree
    """

    family: str
    id: str
    mother: str | None = None
    father: str | None = None
    sex: str
    affected: str
    ext_id: str = 'Missing'
    hpo_terms: list[PhenoPacketHpo] = Field(default_factory=list)
