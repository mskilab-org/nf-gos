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

def recursive_touch(path):
    """
    Creates a file at the given path recursively, including any necessary parent directories.
    If the file exists, it updates the modification time.
    """
    # Create parent directories if they don't exist
    os.makedirs(os.path.dirname(path), exist_ok=True)
    
    # Create or update the file
    with open(path, 'a'):
        os.utime(path, None)


def main(args):
    input_dir = args.input_directory
    output_sbs = './sbs_results'
    output_indel = './indel_results'
    sigmat_project = 'sigmat_results'
    genome = args.genome or 'GRCh37'
    cosmic_version = args.cosmic_version or 3.4
    
    sbs_touch = [
        "sbs_results/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
        "sbs_results/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed_Mutation_Probabilities.txt",
        "sig_inputs/output/SBS/Input_vcffiles.SBS1536.all",
        "sig_inputs/output/SBS/Input_vcffiles.SBS18.all",
        "sig_inputs/output/SBS/Input_vcffiles.SBS24.all",
        "sig_inputs/output/SBS/Input_vcffiles.SBS288.all",
        "sig_inputs/output/SBS/Input_vcffiles.SBS384.all",
        "sig_inputs/output/SBS/Input_vcffiles.SBS4608.all",
        "sig_inputs/output/SBS/Input_vcffiles.SBS6.all",
        "sig_inputs/output/SBS/Input_vcffiles.SBS6144.all",
        "sig_inputs/output/SBS/Input_vcffiles.SBS96.all",
        "sig_inputs/output/SBS/sigmat_results.SBS1536.all",
        "sig_inputs/output/SBS/sigmat_results.SBS18.all",
        "sig_inputs/output/SBS/sigmat_results.SBS24.all",
        "sig_inputs/output/SBS/sigmat_results.SBS288.all",
        "sig_inputs/output/SBS/sigmat_results.SBS384.all",
        "sig_inputs/output/SBS/sigmat_results.SBS4608.all",
        "sig_inputs/output/SBS/sigmat_results.SBS6.all",
        "sig_inputs/output/SBS/sigmat_results.SBS6144.all",
        "sig_inputs/output/SBS/sigmat_results.SBS96.all"
    ]
    
    indel_touch = [
        "indel_results/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
        "indel_results/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed_Mutation_Probabilities_WG-25-08___vs___WG-25-80_somatic.txt",
        "sig_inputs/output/ID/Input_vcffiles.ID28.all",
        "sig_inputs/output/ID/Input_vcffiles.ID332.all",
        "sig_inputs/output/ID/Input_vcffiles.ID415.all",
        "sig_inputs/output/ID/Input_vcffiles.ID83.all",
        "sig_inputs/output/ID/Input_vcffiles.ID8628.all",
        "sig_inputs/output/ID/Input_vcffiles.ID96.all",
        "sig_inputs/output/ID/sigmat_results.ID28.all",
        "sig_inputs/output/ID/sigmat_results.ID332.all",
        "sig_inputs/output/ID/sigmat_results.ID415.all",
        "sig_inputs/output/ID/sigmat_results.ID83.all",
        "sig_inputs/output/ID/sigmat_results.ID8628.all",
        "sig_inputs/output/ID/sigmat_results.ID96.all",
        "sig_inputs/output/vcf_files/ID/10_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/11_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/12_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/13_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/14_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/15_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/16_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/17_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/18_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/19_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/1_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/20_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/21_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/22_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/2_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/3_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/4_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/5_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/6_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/7_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/8_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/9_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/MT_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/X_seqinfo.txt",
        "sig_inputs/output/vcf_files/ID/Y_seqinfo.txt"
    ]


    # run twice for SBS and indel
    print('running sbs...')
    try:
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
    except:
        print(f"SNVs might be empty")
        # Touch files to ensure they exist    
        for f in sbs_touch:
            recursive_touch(f)

    print('running indel...')
    try:
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
    except:
        for f in indel_touch:
            recursive_touch(f)


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
