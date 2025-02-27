import os
import sys
import json

def create_params_file():
    default_input = "./samplesheet.csv"
    default_outdir = "./results/"
    default_tools = "all"
    default_genome = "hg19"
    default_email = "example_email@gmail.com"

    # Check if default input exists
    input_exists = os.path.exists(default_input)
    input_status = "found!" if input_exists else "not found!"

    # Prompt for input
    input_prompt = f"Enter input file path [{default_input} ({input_status})]: "
    input_path = input(input_prompt).strip()
    if not input_path:
        if input_exists:
            input_path = default_input
        else:
            print("Error: You must provide a path to the samplesheet.")
            sys.exit(1)
    else:
        if not os.path.exists(input_path):
            print(f"Error: The provided input file '{input_path}' does not exist.")
            sys.exit(1)

    # Check if default outdir exists
    outdir_exists = os.path.exists(default_outdir)
    outdir_status = "found!" if outdir_exists else "not found!"

    # Prompt for outdir
    outdir_prompt = f"Enter output directory [{default_outdir} ({outdir_status})]: "
    outdir = input(outdir_prompt).strip() or default_outdir

    # Prompt for tools
    tools_list = "[ aligner, bamqc, gridss, amber, fragcounter, dryclean, cbs, sage, purple, jabba, non_integer_balance, lp_phased_balance, events, fusions, snpeff, snv_multiplicity, signatures, hrdetect ]"
    tools_prompt = (
        f"Enter tools to use (comma-separated) or leave blank for default [{default_tools}]:\n"
        f"Available tools: {tools_list}\n"
        "Hint: You can pick any combination of tools, or leave it as 'all' which is recommended.\n"
    )
    tools = input(tools_prompt).strip() or default_tools

    # Prompt for genome
    genome_prompt = (
        f"Enter genome [{default_genome}]: (options: hg19, hg38)\n"
        "Hint: Choose 'hg19' or 'hg38'.\n"
    )
    genome_input = input(genome_prompt).strip() or default_genome
    genome_map = {'hg19': 'GATK.GRCh37', 'hg38': 'GATK.GRCh38'}
    genome = genome_map.get(genome_input.lower())
    if not genome:
        print(f"Warning: Invalid genome '{genome_input}'. Using default '{default_genome}'.")
        genome = genome_map[default_genome]

    # Prompt for email
    email_prompt = f"Enter email address [{default_email}]: "
    email = input(email_prompt).strip() or default_email

    # Create params dictionary
    params = {
        "input": input_path,
        "outdir": outdir,
        "tools": tools,
        "genome": genome,
        "email": email
    }

    # Write to params.json
    with open("params.json", "w") as f:
        json.dump(params, f, indent=4)
    print("Created 'params.json' with the provided parameters.")
