import argparse
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator import install as genInstall

# Usage
# python sigprofilerassgnment.py --data-directory /path/to/vcf/files --output-directory /path/to/output --maximum-signatures 25 --cosmic-version 3.2


def main(args):
    sample = args.input_vcf
    output = args.output_directory or './'
    genome = args.genome or 'GRCh37'
    cosmic_version = args.cosmic_version or 3.4
    print('running...')
    Analyze.cosmic_fit(
        sample,
        output,
        input_type="vcf",
        context_type="96",
        collapse_to_SBS96=True,
        cosmic_version=cosmic_version,
        exome=False,
        genome_build=genome,
        signature_database=None,
        exclude_signature_subgroups=None,
        export_probabilities=False,
        export_probabilities_per_mutation=False,
        make_plots=False,
        sample_reconstruction_plots=False,
        verbose=True
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run SigProfilerAssigner with a CLI interface.')
    parser.add_argument('--genome', type=str, default='GRCh37',
                        help='Genome reference to download. Default is GRCh37.')
    parser.add_argument('--input-vcf', type=str, required=True,
                        help='Path to the input VCF file.')
    parser.add_argument('--output-directory', type=str, required=True,
                        help='Directory to save the output results.')
    parser.add_argument('--cosmic-version', type=float, default=3.4,
                        help='COSMIC database version to use. Default is 3.4.')

    args = parser.parse_args()
    main(args)
