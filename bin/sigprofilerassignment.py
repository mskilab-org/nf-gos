import argparse
import os
import gzip
import shutil
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# Usage
# python sigprofilerassgnment.py --data-directory /path/to/vcf/files --output-directory /path/to/output --maximum-signatures 25 --cosmic-version 3.2

def none_or_str(value):
    return None if value == "None" else value

def unzip_vcf(input_vcf, input_dir):
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)

    # Extract the base name of the file (without .gz)
    base_name = os.path.basename(input_vcf).rsplit('.', 1)[0]
    unzipped_vcf_path = os.path.join(input_dir, base_name)

    with gzip.open(input_vcf, 'rb') as f_in:
        with open(unzipped_vcf_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return unzipped_vcf_path


def main(args):
    input_dir = args.input_directory
    output_sbs = './sbs_results'
    output_indel = './indel_results'
    sigmat_project = 'sigmat_results'
    genome = args.genome or 'GRCh37'
    cosmic_version = args.cosmic_version or 3.4

    # run twice for SBS and indel
    print('running sbs...')
    Analyze.cosmic_fit(
        samples=input_dir,
        output=output_sbs,
        input_type="vcf",
        context_type="96",
        collapse_to_SBS96=True,
        cosmic_version=cosmic_version,
        exome=False,
        genome_build=genome,
        signature_database=args.sbs_signatures,
        exclude_signature_subgroups=None,
        export_probabilities=True,
        export_probabilities_per_mutation=True,
        make_plots=False,
        sample_reconstruction_plots=False,
        verbose=True
    )

    print('running indel...')
    Analyze.cosmic_fit(
        samples=input_dir,
        output=output_indel,
        input_type="vcf",
        context_type="ID",
        collapse_to_SBS96=False,
        cosmic_version=cosmic_version,
        exome=False,
        genome_build=genome,
        signature_database=args.id_signatures,
        export_probabilities=True,
        export_probabilities_per_mutation=True,
        make_plots=False,
        sample_reconstruction_plots=False,
        verbose=True
    )

    print('running matrix generator...')
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        sigmat_project,
        genome,
        input_dir,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run SigProfilerAssigner with a CLI interface.')
    parser.add_argument('--genome', type=str, default='GRCh37',
                        help='Genome reference to download. Default is GRCh37.')
    parser.add_argument('--input-vcf', type=str, required=True,
                        help='Path to the input gzipped VCF file.')
    parser.add_argument('--cosmic-version', type=float, default=3.4,
                        help='COSMIC database version to use. Default is 3.4.')
    parser.add_argument('--id-signatures', type=none_or_str, default=None,
                        help='Custom id signature matrix.')
    parser.add_argument('--sbs-signatures', type=none_or_str, default=None,
                        help='Custom sbs signature matrix.')

    args = parser.parse_args()

    # Create input directory and unzip the VCF file
    input_dir = os.path.join(os.getcwd(), 'sig_inputs')
    unzipped_vcf_path = unzip_vcf(args.input_vcf, input_dir)

    # Update args to use the input directory
    args.input_directory = input_dir

    main(args)
