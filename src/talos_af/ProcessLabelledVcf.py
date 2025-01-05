"""
This is the TalosAF equivalent of Talos' ValidateMoi (which I hate as a name)
We read the VCF, then parse the variants contained to determine which are a match for the ACMG reportable variants
We then write the results to a file

Due to the nature of this process (based on a very limited number of variants in a limited number of genes) we don't
need to do any tricky fragmented/parallelised uparsing, we just read the whole VCF in one go, then process the variants
in each gene according to its own rules
"""

from argparse import ArgumentParser
from itertools import combinations_with_replacement

from peds import open_ped

from talos_af.config import config_retrieve
from talos_af.logger import get_logger
from talos_af.models import (
    AFCategoryC,
    AFSpecification,
    MinimalVariant,
    ReportableVariant,
    SmallVariant,
)
from talos_af.utils import (
    gather_gene_dict,
    read_json_from_path,
)


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to the input VCF', required=True)
    parser.add_argument('--output', help='Where to write the output', required=True)
    parser.add_argument('--pedigree', help='Path to the pedigree file', required=True)
    parser.add_argument('--af_spec', help='Path to JSON containing AF specification', required=True)
    args = parser.parse_args()

    main(
        input_vcf=args.input,
        output=args.output,
        pedigree=args.pedigree,
        af_spec=args.af_spec,
    )


def type_a_gene(gene: str, variants: list[SmallVariant]) -> list[ReportableVariant]:
    """
    simplest version - we check for any variant (dominant is tolerated)

    Args:
        gene (str): gene ID
        variants (list[SmallVariant]): list of variants to check

    Returns:
        list[ReportableVariant]: a list of all the variants which are type A (monoallelic/dominant)
    """

    outgoing_results: list[ReportableVariant] = []

    for variant in variants:
        # todo embed a more stringent AF/consequence check here?

        for sample in variant.het_samples:
            # create the minimal representation
            minrep_var = MinimalVariant(
                gene=gene,
                coordinates=variant.coordinates,
                info=variant.info,
                genotype='Het',
            )
            outgoing_results.append(
                ReportableVariant(
                    sample=sample,
                    var_1=minrep_var,
                    gene=gene,
                    rule='Type A',
                )
            )

        for sample in variant.hom_samples:
            minrep_var = MinimalVariant(
                gene=gene,
                coordinates=variant.coordinates,
                info=variant.info,
                genotype='Hom',
            )
            outgoing_results.append(
                ReportableVariant(
                    sample=sample,
                    var_1=minrep_var,
                    gene=gene,
                    rule='Type A',
                )
            )

    return outgoing_results


def type_b_gene(gene: str, variants: list[SmallVariant]) -> list[ReportableVariant]:
    """
    fairly simple version - we check for any variant with 2-hits (hom or comp-het)

    Args:
        gene (str): gene ID
        variants (list[SmallVariant]): list of variants to check

    Returns:
        list[ReportableVariant]: a list of all the variants which are type B (recessive/bi-allelic)
    """

    outgoing_results: list[ReportableVariant] = []

    # for homozygotes, we just iterate directly over the variants and look for Homs
    for variant in variants:
        for sample in variant.hom_samples:
            minrep_var = MinimalVariant(
                gene=gene,
                coordinates=variant.coordinates,
                info=variant.info,
                genotype='Hom',
            )
            outgoing_results.append(
                ReportableVariant(
                    sample=sample,
                    var_1=minrep_var,
                    gene=gene,
                    rule='Type B (Hom)',
                )
            )

    # generate all possible comp-het permutations
    for variant_1, variant_2 in combinations_with_replacement(variants, 2):
        if common_samples := variant_1.het_samples & variant_2.het_samples:
            for sample in common_samples:
                minrep_var_1 = MinimalVariant(
                    gene=gene,
                    coordinates=variant_1.coordinates,
                    info=variant_1.info,
                    genotype='Het',
                )
                minrep_var_2 = MinimalVariant(
                    gene=gene,
                    coordinates=variant_2.coordinates,
                    info=variant_2.info,
                    genotype='Het',
                )
                outgoing_results.append(
                    ReportableVariant(
                        sample=sample,
                        var_1=minrep_var_1,
                        var_2=minrep_var_2,
                        gene=gene,
                        rule='Type B (Comp-Het)',
                    ),
                )

    return outgoing_results


def type_c_gene(gene: str, af_details: AFCategoryC, variants: list[SmallVariant]) -> list[ReportableVariant]:
    """
    For genes assigned a category C, we check for a range of specific conditions

    Args:
        gene ():
        af_details ():
        variants ():

    Returns:

    """

    result_list: list[ReportableVariant] = []


def main(input_vcf: str, output: str, pedigree: str, af_spec: str):
    """
    Read the VCF, and write the results to a file
    """

    # read the pedigree file - why?
    pedigree = open_ped(pedigree)

    # read the AF specification
    af_spec = read_json_from_path(af_spec, return_model=AFSpecification)

    # get a dict of all variants indexed by gene
    gene_dict = gather_gene_dict(input_vcf)

    results = []
    for gene, variants in gene_dict.items():
        gene_af_spec = af_spec.genes.get(gene)

        # get results from the variants, spec, and pedigree
        if gene_af_spec.a:
            results.extend(type_a_gene(gene, variants))
        if gene_af_spec.b:
            results.extend(type_b_gene(gene, variants))
        if gene_af_spec.c:
            results.extend(type_c_gene(gene, gene_af_spec.c, variants))


if __name__ == '__main__':
    cli_main()