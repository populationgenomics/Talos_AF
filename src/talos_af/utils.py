"""
classes and methods shared across reanalysis components

HTTPX requests are backoff-wrapped using tenaticy
https://tenacity.readthedocs.io/en/latest/
"""

import httpx
import json
import re
import zoneinfo
from collections import defaultdict
from datetime import datetime
from itertools import chain, combinations_with_replacement, islice
from pathlib import Path
from typing import Any

import cyvcf2
from cloudpathlib.anypath import to_anypath
from tenacity import retry, stop_after_attempt, wait_exponential_jitter, retry_if_exception_type

from talos_af.config import config_retrieve
from talos_af.models import (
    Coordinates,
    FileTypes,
    SmallVariant,
)
from talos_af.logger import get_logger

HOMREF: int = 0
HETALT: int = 1
UNKNOWN: int = 2
HOMALT: int = 3

# in cyVCF2, these ints represent HOMREF, and UNKNOWN
BAD_GENOTYPES: set[int] = {HOMREF, UNKNOWN}
PHASE_SET_DEFAULT = -2147483648
NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
X_CHROMOSOME = {'X'}

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
TODAY = datetime.now(tz=TIMEZONE).strftime('%Y-%m-%d_%H:%M')

DATE_RE = re.compile(r'\d{4}-\d{2}-\d{2}')


def chunks(iterable, chunk_size):
    """
    Yield successive n-sized chunks from an iterable

    Args:
        iterable (): any iterable - tuple, str, list, set
        chunk_size (): size of intervals to return

    Returns:
        intervals of requested size across the collection
    """

    if isinstance(iterable, set):
        iterable = list(iterable)

    for i in range(0, len(iterable), chunk_size):
        yield iterable[i : (i + chunk_size)]


def generator_chunks(generator, size):
    """
    Iterates across a generator, returning specifically sized chunks

    Args:
        generator (): any generator or method implementing yield
        size (): size of iterator to return

    Returns:
        a subset of the generator results
    """
    iterator = iter(generator)
    for first in iterator:
        yield list(chain([first], islice(iterator, size - 1)))


def identify_file_type(file_path: str) -> FileTypes | Exception:
    """
    return type of the file, if present in FileTypes enum

    Args:
        file_path (str):

    Returns:
        A matching file type, or die
    """
    pl_filepath = Path(file_path)

    # pull all extensions (e.g. .vcf.bgz will be split into [.vcf, .bgz]
    if not (extensions := pl_filepath.suffixes):
        raise ValueError('cannot identify input type from extensions')

    if extensions[-1] == '.ht':
        return FileTypes.HAIL_TABLE
    if extensions[-1] == '.mt':
        return FileTypes.MATRIX_TABLE
    if extensions == ['.vcf']:
        return FileTypes.VCF
    if extensions == ['.vcf', '.gz']:
        return FileTypes.VCF_GZ
    if extensions == ['.vcf', '.bgz']:
        return FileTypes.VCF_BGZ
    raise TypeError(f'File cannot be definitively typed: {extensions}')


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential_jitter(initial=1, max=5, exp_base=2),
    retry=retry_if_exception_type(
        (
            httpx.ReadTimeout,
            httpx.ConnectError,
            httpx.TooManyRedirects,
            httpx.RequestError,
        ),
    ),
    reraise=True,
)
def get_json_response(url):
    """
    takes a request URL, checks for healthy response, returns the JSON
    For this purpose we only expect a dictionary return
    List use-case (activities endpoint) no longer supported

    for backoff handling I'm using the tenacity library

    Args:
        url (str): URL to retrieve JSON format data from

    Returns:
        the JSON response from the endpoint
    """
    response = httpx.get(url, headers={'Accept': 'application/json'}, timeout=60, follow_redirects=True)
    if response.is_success:
        return response.json()
    raise ValueError('The JSON response could not be parsed successfully')


def get_phase_data(samples, var) -> dict[str, dict[int, str]]:
    """
    read phase data from this variant

    Args:
        samples ():
        var ():
    """
    phased_dict: dict[str, dict[int, str]] = defaultdict(dict)

    # first set the numpy.ndarray to be a list of ints
    # then zip against ordered sample IDs
    # this might need to store the exact genotype too
    # i.e. 0|1 and 1|0 can be in the same phase-set
    # but are un-phased variants

    try:
        for sample, phase, genotype in zip(samples, map(int, var.format('PS')), var.genotypes):
            # cyvcf2.Variant holds two ints, and a bool
            allele_1, allele_2, phased = genotype
            if not phased:
                continue
            gt = f'{allele_1}|{allele_2}'
            # phase set is a number
            if phase != PHASE_SET_DEFAULT:
                phased_dict[sample][phase] = gt
    except KeyError as ke:
        get_logger().info('failed to find PS phase attributes')
        try:
            # retry using PGT & PID
            for sample, phase_gt, phase_id in zip(samples, var.format('PGT'), var.format('PID')):
                if phase_gt != '.' and phase_id != '.':
                    phased_dict[sample][phase_id] = phase_gt

        except KeyError as ke2:
            get_logger().info('also failed using PID and PGT')
            raise ke from ke2

    return dict(phased_dict)


def create_small_variant(var: cyvcf2.Variant, samples: list[str]):
    """
    takes a small variant and creates a Model from it

    Args:
        var ():
        samples ():
    """

    coordinates = Coordinates(chrom=var.CHROM.replace('chr', ''), pos=var.POS, ref=var.REF, alt=var.ALT[0])
    depths: dict[str, int] = dict(zip(samples, map(int, var.gt_depths)))
    info: dict[str, Any] = {x.lower(): y for x, y in var.INFO} | {'seqr_link': coordinates.string_format}

    # optionally - ignore some categories from this analysis
    if ignore_cats := config_retrieve(['ValidateMOI', 'ignore_categories'], []):
        info = {key: val for key, val in info.items() if key not in ignore_cats}

    het_samples, hom_samples = get_non_ref_samples(variant=var, samples=samples)

    # set the class attributes
    boolean_categories = [key for key in info if key.startswith('categoryboolean')]
    sample_categories = [key for key in info if key.startswith('categorysample')]
    support_categories = [key for key in info if key.startswith('categorysupport')]

    # overwrite with true booleans
    for cat in support_categories + boolean_categories:
        info[cat] = info.get(cat, 0) == 1

    # sample categories are a set of strings or 'missing'
    for sam_cat in sample_categories:
        if isinstance(info[sam_cat], str):
            info[sam_cat] = info[sam_cat].split(',') if info[sam_cat] != 'missing' else set()
        if isinstance(info[sam_cat], list):
            info[sam_cat] = set(info[sam_cat])

    phased = get_phase_data(samples, var)
    ab_ratios = dict(zip(samples, map(float, var.gt_alt_freqs)))
    transcript_consequences = extract_csq(csq_contents=info.pop('csq', ''))

    return SmallVariant(
        coordinates=coordinates,
        info=info,
        het_samples=het_samples,
        hom_samples=hom_samples,
        sample_support=support_categories,
        phased=phased,
        depths=depths,
        ab_ratios=ab_ratios,
        transcript_consequences=transcript_consequences,
    )


# CompHetDict structure: {sample: {variant_string: [variant, ...]}}
# sample: string, e,g, CGP12345
CompHetDict = dict[str, dict[str, list[SmallVariant]]]
GeneDict = dict[str, list[SmallVariant]]


def canonical_contigs_from_vcf(reader) -> set[str]:
    """
    read the header fields from the VCF handle
    return a set of all 'canonical' contigs

    Args:
        reader (cyvcf2.VCFReader):
    """

    # contig matching regex - remove all HLA/decoy/unknown
    contig_re = re.compile(r'^(chr)?[0-9XYMT]{1,2}$')

    return {
        contig['ID']
        for contig in reader.header_iter()
        if contig['HeaderType'] == 'CONTIG' and re.match(contig_re, contig['ID'])
    }


def gather_gene_dict_from_contig(
    contig: str,
    variant_source,
) -> GeneDict:
    """
    takes a cyvcf2.VCFReader instance, and a specified chromosome
    iterates over all variants in the region, and builds a lookup

    optionally takes a second VCF and incorporates into same dict

    Args:
        contig (): contig name from VCF header
        variant_source (): the VCF reader instance

    Returns:
        A lookup in the form
        {
            gene1: [var1, var2],
            gene2: [var3],
            ...
        }
    """

    # a dict to allow lookup of variants on this whole chromosome
    contig_variants = 0
    contig_dict = defaultdict(list)

    # iterate over all variants on this contig and store by unique key
    # if contig has no variants, prints an error and returns []
    for variant in variant_source(contig):
        small_variant = create_small_variant(
            var=variant,
            samples=variant_source.samples,
        )

        # update the variant count
        contig_variants += 1

        # update the gene index dictionary
        contig_dict[small_variant.info.get('gene_id')].append(small_variant)

    get_logger().info(f'Contig {contig} contained {contig_variants} variants, in {len(contig_dict)} genes')

    return contig_dict


def read_json_from_path(read_path: str | None = None, default: Any = None, return_model: Any = None) -> Any:
    """
    take a path to a JSON file, read into an object
    if the path doesn't exist - return the default object
    uses cloudpath to be deployment agnostic

    Args:
        read_path (str): where to read from - if None... will return the value "default"
        default (Any):
        return_model (pydantic Models): any Model to read/validate as

    Returns:
        either the object from the JSON file, or None
    """

    if read_path is None:
        get_logger().error('read_json_from_path was passed the path "None"')
        return default

    assert isinstance(read_path, str)
    read_anypath = to_anypath(read_path)

    if not read_anypath.exists():
        get_logger().error(f'{read_path} did not exist')
        return default

    with read_anypath.open() as handle:
        json_data = json.load(handle)
        if return_model:
            return return_model.model_validate(json_data)
        return json_data


def get_non_ref_samples(variant, samples: list[str]) -> tuple[set[str], set[str]]:
    """
    for this variant, find all samples with a call
    cyvcf2 uses 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT

    Args:
        variant (cyvcf2.Variant):
        samples (list[str]):

    Returns:
        2 sets of strings; het and hom
    """
    het_samples = set()
    hom_samples = set()

    # this iteration is based on the cyvcf2 representations
    for sam, genotype_int in zip(samples, variant.gt_types):
        if genotype_int in BAD_GENOTYPES:
            continue
        if genotype_int == HETALT:
            het_samples.add(sam)
        if genotype_int == HOMALT:
            hom_samples.add(sam)

    return het_samples, hom_samples


def extract_csq(csq_contents: str) -> list[dict]:
    """
    handle extraction of the CSQ entries based on string in config

    Args:
        csq_contents ():
    """

    # allow for no CSQ data, i.e. splice variant
    if not csq_contents:
        return []

    # break mono-CSQ-string into components
    csq_categories = config_retrieve(['RunHailFiltering', 'csq_string'])

    # iterate over all consequences, and make each into a dict
    txc_dict = [dict(zip(csq_categories, each_csq.split('|'), strict=True)) for each_csq in csq_contents.split(',')]

    # update this String to be either a float, or missing
    for each_dict in txc_dict:
        am_path = each_dict.get('am_pathogenicity')
        each_dict['am_pathogenicity'] = float(am_path) if am_path else ''

    return txc_dict


# todo: reinstate peds pedigree reading
def find_comp_hets(var_list: list[SmallVariant], pedigree) -> CompHetDict:
    """
    manual implementation to find compound hets
    variants provided in the format

    [var1, var2, ..]

    generate pair content in the form
    {sample: {var_as_string: [partner_variant, ...]}}

    Args:
        var_list (list[SmallVariant]): all variants in this gene
        pedigree (): Pedigree
    """

    # create an empty dictionary
    comp_het_results: CompHetDict = defaultdict(dict)

    # use combinations_with_replacement to find all gene pairs
    for var_1, var_2 in combinations_with_replacement(var_list, 2):
        assert var_1.coordinates.chrom == var_2.coordinates.chrom

        if (var_1.coordinates == var_2.coordinates) or var_1.coordinates.chrom in NON_HOM_CHROM:
            continue

        # iterate over any samples with a het overlap
        for sample in var_1.het_samples.intersection(var_2.het_samples):
            phased = False

            # don't assess male compound hets on sex chromosomes
            if pedigree.by_id[sample].sex == '1' and var_1.coordinates.chrom in X_CHROMOSOME:
                continue

            # check for both variants being in the same phase set
            if sample in var_1.phased and sample in var_2.phased:
                # check for presence of the same phase set
                for phase_set in [ps for ps in var_1.phased[sample] if ps in var_2.phased[sample]]:
                    if var_1.phased[sample][phase_set] == var_2.phased[sample][phase_set]:
                        phased = True
            if not phased:
                comp_het_results[sample].setdefault(var_1.coordinates.string_format, []).append(var_2)
                comp_het_results[sample].setdefault(var_2.coordinates.string_format, []).append(var_1)

    return comp_het_results
