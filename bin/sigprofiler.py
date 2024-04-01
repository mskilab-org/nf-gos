import argparse
from SigProfilerExtractor import sigpro as sig
from SigProfilerMatrixGenerator import install as genInstall

# Usage
# python sigprofiler.py --data-directory /path/to/vcf/files --output-directory /path/to/output --maximum-signatures 25 --cosmic-version 3.2


def main(args):
    genome = args.genome
    print(f'downloading {genome} reference')
    genInstall.install(genome, rsync=False, bash=True)
    print('running...')
    data = args.data_directory
    sig.sigProfilerExtractor(
        "vcf",
        args.output_directory,
        data,
        maximum_signatures=args.maximum_signatures,
        cosmic_version=args.cosmic_version,
        gpu=args.gpu
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run SigProfilerExtractor with a CLI interface.')
    parser.add_argument('--genome', type=str, default='GRCh37',
                        help='Genome reference to download. Default is GRCh37.')
    parser.add_argument('--data-directory', type=str, required=True,
                        help='Directory containing the VCF files.')
    parser.add_argument('--output-directory', type=str, required=True,
                        help='Directory to save the output results.')
    parser.add_argument('--maximum-signatures', type=int, default=30,
                        help='Maximum number of signatures to extract. Default is 30.')
    parser.add_argument('--cosmic-version', type=float, default=3.2,
                        help='COSMIC database version to use. Default is 3.2.')
    parser.add_argument('--gpu', action='store_true',
                        help='Use GPU acceleration if available.')

    args = parser.parse_args()
    main(args)
