"""
Read, filter, annotate, classify, and write Genetic data
- read MT
- read PanelApp data through GCP client
- hard-filter (FILTERS, AC)
- extract generic fields
- remove all rows and consequences not relevant to GREEN genes
- extract vep data into CSQ string(s)
- annotate with categories 1, 2, 3, 4, 5, 6, and Support
- remove un-categorised variants
- write as VCF
"""

from argparse import ArgumentParser

import hail as hl

from talos_af.config import config_retrieve
from talos_af.logger import get_logger
from talos_af.models import AFSpecification
from talos_af.utils import read_json_from_path

# set some Hail constants
MISSING_INT = hl.int32(0)
MISSING_FLOAT_LO = hl.float64(0.0)
MISSING_STRING = hl.str('missing')
ONE_INT = hl.int32(1)
BENIGN = hl.str('benign')
LOFTEE_HC = hl.str('HC')
PATHOGENIC = hl.str('pathogenic')

# decide whether to repartition the data before processing starts
MAX_PARTITIONS = 10000


def annotate_clinvarbitration(mt: hl.MatrixTable, clinvar: str) -> hl.MatrixTable:
    """
    Don't allow these annotations to be missing
    - Talos has been co-developed with ClinvArbitration, a ClinVar re-summary effort
    - We no longer permit this to be missing (this has slipped in the past, causing odd results)

    See: https://github.com/populationgenomics/ClinvArbitration

    We replace any existing ClinVar annotations with our own version, then identify Pathogenic/Benign variants

    Args:
        mt (): the MatrixTable of all variants
        clinvar (): the table of private ClinVar annotations to use

    Returns:
        The same MatrixTable but with additional annotations
    """

    get_logger().info(f'loading private clinvar annotations from {clinvar}')
    ht = hl.read_table(clinvar)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            clinvar_significance=hl.or_else(ht[mt.row_key].clinical_significance, MISSING_STRING),
            clinvar_stars=hl.or_else(ht[mt.row_key].gold_stars, MISSING_INT),
            clinvar_allele=hl.or_else(ht[mt.row_key].allele_id, MISSING_INT),
        ),
    )

    # remove all confidently benign
    mt = mt.filter_rows(
        (mt.info.clinvar_significance.lower().contains(BENIGN)) & (mt.info.clinvar_stars > 0),
        keep=False,
    )

    # annotate as either strong or regular, return the result
    return mt.annotate_rows(
        info=mt.info.annotate(
            clinvarbitration=hl.if_else(
                mt.info.clinvar_significance.lower().contains(PATHOGENIC),
                ONE_INT,
                MISSING_INT,
            ),
            clinvarbitration_strong=hl.if_else(
                (mt.info.clinvar_significance.lower().contains(PATHOGENIC)) & (mt.info.clinvar_stars > 0),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def extract_annotations(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    pull out select fields which aren't per-consequence
    store these in INFO (required to be included in VCF export)
    replace with placeholder (least consequential) if empty
    e.g. most tools score 0, but for Sift 1 is the least important

    Args:
        mt ():
    Returns:
        Same matrix with re-positioned attributes
    """

    get_logger().info('Pulling VEP annotations into INFO field')

    return mt.annotate_rows(
        info=mt.info.annotate(
            gnomad_ex_af=hl.or_else(mt.gnomad_exomes.AF, MISSING_FLOAT_LO),
            gnomad_ex_an=hl.or_else(mt.gnomad_exomes.AN, MISSING_INT),
            gnomad_ex_ac=hl.or_else(mt.gnomad_exomes.AC, MISSING_INT),
            gnomad_ex_hom=hl.or_else(mt.gnomad_exomes.Hom, MISSING_INT),
            gnomad_ex_hemi=hl.or_else(mt.gnomad_exomes.Hemi, MISSING_INT),
            gnomad_af=hl.or_else(mt.gnomad_genomes.AF, MISSING_FLOAT_LO),
            gnomad_an=hl.or_else(mt.gnomad_genomes.AN, MISSING_INT),
            gnomad_ac=hl.or_else(mt.gnomad_genomes.AC, MISSING_INT),
            gnomad_hom=hl.or_else(mt.gnomad_genomes.Hom, MISSING_INT),
            gnomad_hemi=hl.or_else(mt.gnomad_genomes.Hemi, MISSING_INT),
        ),
    )


def quality_filter_matrix(mt: hl.MatrixTable, ac_threshold: float = 0.01, min_callset_ac: int = 5) -> hl.MatrixTable:
    """
    Remove variants with VCF quality filters

    Remove variants with AC in joint-call over threshold
    Will never remove variants with 5 or fewer instances

    Also overridden by having a Clinvar Pathogenic anno.

    Args:
        mt (hl.MatrixTable):
        ac_threshold (float):
        min_callset_ac (int): only AC filter if there are at least this many calls

    Returns:
        MT with all common-in-this-JC variants removed
        (unless overridden by clinvar path)
    """

    return mt.filter_rows(
        # keep variants with no applied filters
        (hl.is_missing(mt.filters) | (mt.filters.length() == 0))
        # keep variants below the threshold frequency
        | ((min_callset_ac >= mt.info.AC[0]) | (ac_threshold > mt.info.AC[0] / mt.info.AN))
        # permit clinvar pathogenic to slip through these filters
        | (mt.info.clinvarbitration == ONE_INT),
    )


def filter_to_population_rare(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    run the rare filter, using Gnomad Exomes and Genomes
    allow clinvar pathogenic to slip through this filter
    """
    # gnomad exomes and genomes below threshold or missing
    # if missing they were previously replaced with 0.0
    # 'semi-rare' as dominant filters will be more strictly filtered later
    rare_af_threshold = config_retrieve(['RunHailFiltering', 'af_semi_rare'])
    return mt.filter_rows(
        (
            (hl.or_else(mt.gnomad_exomes.AF, MISSING_FLOAT_LO) < rare_af_threshold)
            & (hl.or_else(mt.gnomad_genomes.AF, MISSING_FLOAT_LO) < rare_af_threshold)
        )
        | (mt.info.clinvarbitration == ONE_INT),
    )


def split_rows_by_gene_and_filter_to_relevant(mt: hl.MatrixTable, relevant_genes: hl.SetExpression) -> hl.MatrixTable:
    """
    splits each GeneId onto a new row, then filters any rows not containing annotation on a relevant gene

    - first explode the matrix, separate gene per row
    - throw away all rows without a green gene
    - on all remaining rows, filter transcript consequences to match _this_ gene

    this is the single most powerful filtering row, effectively leaving just the genes we're interested in

    Args:
        mt ():
        relevant_genes (): set of all relevant genes
    Returns:
        exploded MatrixTable
    """

    # split each gene onto a separate row, transforms 'geneIds' field from set to string
    mt = mt.explode_rows(mt.geneIds)

    # filter rows without a green gene (removes empty geneIds)
    mt = mt.filter_rows(relevant_genes.contains(mt.geneIds))

    # limit the per-row transcript CSQ to those relevant to the single
    # gene now present on each row
    return mt.annotate_rows(
        vep=mt.vep.annotate(
            transcript_consequences=mt.vep.transcript_consequences.filter(
                lambda x: (mt.geneIds == x.gene_id)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM'))),
            ),
        ),
    )


def identify_high_impact_variants(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    identifies high impact variants, based on the following criteria:
    - Critical protein consequence on at least one transcript
    - either predicted NMD or
    - any star Pathogenic or Likely_pathogenic in Clinvar

    Args:
        mt (hl.MatrixTable):
    Returns:
        same variants, high_impact set to 1 or 0
    """

    critical_consequences = hl.set(config_retrieve(['RunHailFiltering', 'critical_csq']))

    # First check if we have any HIGH consequences
    # then explicitly link the LOFTEE check with HIGH consequences
    # OR allow for a pathogenic ClinVar, any Stars
    return mt.annotate_rows(
        info=mt.info.annotate(
            high_impact=hl.if_else(
                (
                    hl.len(
                        mt.vep.transcript_consequences.filter(
                            lambda x: (hl.len(critical_consequences.intersection(hl.set(x.consequence_terms))) > 0),
                        ),
                    )
                    > 0
                )
                & (
                    (
                        hl.len(
                            mt.vep.transcript_consequences.filter(
                                lambda x: (hl.len(critical_consequences.intersection(hl.set(x.consequence_terms))) > 0)
                                & ((x.lof == LOFTEE_HC) | (hl.is_missing(x.lof))),
                            ),
                        )
                        > 0
                    )
                    | (mt.info.clinvarbitration == ONE_INT)
                ),
                ONE_INT,
                MISSING_INT,
            ),
        ),
    )


def vep_struct_to_csq(vep_expr: hl.expr.StructExpression) -> hl.expr.ArrayExpression:
    """
    Taken shamelessly from the gnomad library source code
    Given a VEP Struct, returns an array of VEP VCF CSQ strings
    (1 per csq in the struct).
    Fields & order correspond to those in `csq_fields`, corresponding to the
    VCF header that is required to interpret the VCF CSQ INFO field.
    Order is flexible & all fields in the default value are supported.
    These fields are formatted in the same way that their VEP CSQ counterparts are.

    Args:
        vep_expr (hl.Struct):
    Returns:
        generates an array of Strings for each CSQ
    """

    def get_csq_from_struct(element: hl.expr.StructExpression) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update(
            {
                'consequence': hl.delimit(element.consequence_terms, delimiter='&'),
                'feature': element.transcript_id,
                'ensp': element.protein_id,
                'gene': element.gene_id,
                'symbol': element.gene_symbol,
                'mane_select': element.mane_select,
            },
        )

        # pull the required fields and ordering from config
        csq_fields = config_retrieve(['RunHailFiltering', 'csq_string'])

        return hl.delimit([hl.or_else(hl.str(fields.get(f, '')), '') for f in csq_fields], '|')

    csq = hl.empty_array(hl.tstr)
    csq = csq.extend(
        hl.or_else(vep_expr['transcript_consequences'].map(lambda x: get_csq_from_struct(x)), hl.empty_array(hl.tstr)),
    )

    # previous consequence filters may make this caution unnecessary
    return hl.or_missing(hl.len(csq) > 0, csq)


def write_matrix_to_vcf(mt: hl.MatrixTable, vcf_out: str):
    """
    write the remaining MatrixTable content to file as a VCF
    generate a custom header containing the CSQ contents which
    were retained during this run

    Args:
        mt (): the whole MatrixTable
        vcf_out (str): where to write the VCF
    """

    # this temp file needs to be in GCP, not local
    # otherwise the batch that generates the file won't be able to read
    header_path = 'additional_header.txt'

    # generate a CSQ string specific to the config file for decoding later
    csq_contents = '|'.join(config_retrieve(['RunHailFiltering', 'csq_string']))

    # write this custom header locally
    with open(header_path, 'w') as handle:
        handle.write(f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {csq_contents}">')
    get_logger().info(f'Writing categorised variants out to {vcf_out}')
    hl.export_vcf(mt, vcf_out, append_to_header=header_path, tabix=True)


def generate_a_checkpoint(mt: hl.MatrixTable, checkpoint_path: str) -> hl.MatrixTable:
    """
    wrapper around a few repeated lines which write a checkpoint
    Args:
        mt ():
        checkpoint_path (str): where to write the checkpoint

    Returns:

    """
    get_logger().info(f'Checkpointing to {checkpoint_path} after filtering out a ton of variants')
    mt = mt.checkpoint(checkpoint_path)

    # die if there are no variants remaining. Only ever count rows after a checkpoint
    if not (current_rows := mt.count_rows()):
        raise ValueError('No remaining rows to process!')

    get_logger().info(f'Local checkpoint written, {current_rows} rows remain')
    return mt


def cli_main():
    """
    Read MT, filter, and apply category annotation, export as a VCF
    """

    parser = ArgumentParser()
    parser.add_argument('--input', required=True, help='path to input MT')
    parser.add_argument('--output', help='Where to write the VCF', required=True)
    parser.add_argument('--af_spec', help='Path to JSON containing AF specification', required=True)
    parser.add_argument('--clinvar', help='HT containing ClinvArbitration annotations', required=True)
    parser.add_argument('--checkpoint', help='Where/whether to checkpoint, String path', default=None)
    args = parser.parse_args()
    main(
        mt_path=args.input,
        vcf_out=args.output,
        af_spec=args.af_spec,
        clinvar=args.clinvar,
        checkpoint=args.checkpoint,
    )


def main(
    mt_path: str,
    vcf_out: str,
    af_spec: str,
    clinvar: str,
    checkpoint: str | None = None,
):
    """
    Read MT, filter, and apply category annotation, export as a VCF

    Args:
        mt_path (str): Location of the input MT
        vcf_out (str): where to write VCF out
        af_spec (str): location to a JSON containing AF specification
        clinvar (str): location to a ClinVar HT, or unspecified
        checkpoint (str): path to checkpoint data to - serves as checkpoint trigger
    """

    get_logger(__file__).info(
        r"""Welcome To
███████████   █████████   █████          ███████     █████████
█   ███   █  ███     ███   ███         ███     ███  ███     ███
    ███      ███     ███   ███        ███       ███ ███
    ███      ███████████   ███        ███       ███  █████████
    ███      ███     ███   ███        ███       ███         ███
    ███      ███     ███   ███      █  ███     ███  ███     ███
   █████    █████   █████ ███████████    ███████     █████████""",
    )

    # initiate Hail as a local cluster
    number_of_cores = config_retrieve(['RunHailFiltering', 'cores', 'small_variants'], 8)
    get_logger().info(f'Starting Hail with reference genome GRCh38, as a {number_of_cores} core local cluster')

    hl.context.init_spark(master=f'local[{number_of_cores}]', default_reference='GRCh38', quiet=True)

    # read the matrix table from a localised directory
    mt = hl.read_matrix_table(mt_path)
    get_logger().info(f'Loaded annotated MT from {mt_path}, size: {mt.count_rows()}, partitions: {mt.n_partitions()}')

    # repartition if required - local Hail with finite resources has struggled with some really high (~120k) partitions
    # this creates a local duplicate of the input data with far smaller partition counts, for less processing overhead
    if mt.n_partitions() > MAX_PARTITIONS:
        get_logger().info('Shrinking partitions way down with a unshuffled repartition')
        mt = mt.repartition(shuffle=False, n_partitions=number_of_cores * 10)
        if checkpoint:
            get_logger().info('Trying to write the result locally, might need more space on disk...')
            mt = generate_a_checkpoint(mt, f'{checkpoint}_repartitioned')

    # read in the AF specification
    af_spec = read_json_from_path(af_spec, return_model=AFSpecification)

    # pull out the relevant genes from the AFSpecification, cast as a hl.set
    relevant_gene_expression = hl.set(set(af_spec.genes.keys()))

    # remove any rows which have no genes of interest
    mt = split_rows_by_gene_and_filter_to_relevant(mt=mt, relevant_genes=relevant_gene_expression)

    if checkpoint:
        mt = generate_a_checkpoint(mt, f'{checkpoint}_relevant_genes')

    # swap out the default clinvar annotations with private clinvar
    mt = annotate_clinvarbitration(mt=mt, clinvar=clinvar)

    # remove common-in-gnomad variants (also includes ClinVar annotation)
    mt = filter_to_population_rare(mt=mt)

    # filter variants by frequency
    mt = quality_filter_matrix(mt=mt)

    # rearrange the row annotation to make syntax nicer downstream
    mt = extract_annotations(mt=mt)

    # add Labels to the MT based on our selection criteria
    get_logger().info('Identifying variants of interest')
    mt = identify_high_impact_variants(mt=mt)

    # remove all variants without positively assigned labels
    mt = mt.filter_rows((mt.info.clinvarbitration_strong == 1) | (mt.info.high_impact == 1))

    # obtain the massive CSQ string using method stolen from the Broad's Gnomad library
    # also take the single gene_id (from the exploded attribute)
    mt = mt.annotate_rows(info=mt.info.annotate(CSQ=vep_struct_to_csq(mt.vep), gene_id=mt.geneIds))

    write_matrix_to_vcf(mt=mt, vcf_out=vcf_out)


if __name__ == '__main__':
    cli_main()
