import socket
import shutil
import os

def get_environment_defaults():
    nyu_defaults = {
        'pipeline-dir': "/gpfs/data/imielinskilab/git/mskilab/nf-gos",
        'profile': "nyu",
        'nextflow_module': "nextflow/23.04.4",
        'JAVA_HOME': '/gpfs/share/apps/jdk/17u028',
        'aws_module': 'aws-cli'
    }
    nygc_defaults = {
        'pipeline-dir': "/gpfs/commons/groups/imielinski_lab/git/nf-gos",
        'profile': "nygc",
        'nextflow_module': "nextflow/23.10.0",
        'JAVA_HOME': '/nfs/sw/easybuild/software/Java/17.0.6',
        'aws_module': 'awscli'
    }
    mapping = {
        'fn-': nyu_defaults,
        'a100-': nyu_defaults,
        'a40-': nyu_defaults,
        'cn-': nyu_defaults,
        'gn-': nyu_defaults,
        'gpu-': nyu_defaults,
        'bigpurple': nyu_defaults,

        'mskilab0': nygc_defaults,
        'pe2': nygc_defaults,
        'ne1': nygc_defaults,
    }
    hostname = socket.gethostname()
    for prefix, defaults in mapping.items():
        if hostname.startswith(prefix):
            print(f"Detected environment: {hostname}")
            print(f"Using defaults for {hostname}")
            return defaults
    return {}

def load_required_modules(env_defaults):
    required_commands = ['nextflow', 'aws', 'singularity']
    modules_to_load = []
    load_modules_command = ""

    # Check for 'nextflow' command
    if shutil.which('nextflow') is None:
        nextflow_module = env_defaults.get('nextflow_module', 'nextflow')
        if not nextflow_module:
            nextflow_module = 'nextflow'
        modules_to_load.append(nextflow_module)
        print(f"'nextflow' command not found. Loading module '{nextflow_module}'.")
    else:
        print("'nextflow' command is already available.")

    # Check for 'aws' command
    if shutil.which('aws') is None:
        aws_module = env_defaults.get('aws_module', 'aws-cli')
        modules_to_load.append(aws_module)
        print(f"'aws' command not found. Loading module '{aws_module}'.")
    else:
        print("'aws' command is already available.")

    # Check for 'singularity' command
    if shutil.which('singularity') is None:
        modules_to_load.append('singularity')
        print("'singularity' command not found. Loading module 'singularity'.")
    else:
        print("'singularity' command is already available.")

    # Check if JAVA_HOME is set to correct path
    if 'JAVA_HOME' not in os.environ or os.environ['JAVA_HOME'] != env_defaults['JAVA_HOME']:
        print(f"Setting JAVA_HOME to {env_defaults['JAVA_HOME']}")
        os.environ['JAVA_HOME'] = env_defaults['JAVA_HOME']
    else:
        print(f"JAVA_HOME is already set to {env_defaults['JAVA_HOME']}")

    # Build the load modules command string
    for module in modules_to_load:
        load_modules_command += f"module load {module} && "

    return load_modules_command
