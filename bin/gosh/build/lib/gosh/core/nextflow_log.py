import logging

logging.basicConfig(level=logging.WARNING)

# FIELD_NAMES = [
#     'attempt', 'complete', 'container', 'cpus', 'disk', 'duration', 'env', 'error_action',
#     'exit', 'hash', 'hostname', 'inv_ctxt', 'log', 'memory', 'module', 'name', 'native_id',
#     'pcpu', 'peak_rss', 'peak_vmem', 'pmem', 'process', 'queue', 'rchar', 'read_bytes',
#     'realtime', 'rss', 'scratch', 'script', 'start', 'status', 'stderr', 'stdout', 'submit',
#     'syscr', 'syscw', 'tag', 'task_id', 'time', 'vmem', 'vol_ctxt', 'wchar', 'workdir',
#     'write_bytes'
# ]
FIELD_NAMES = [ 'name', 'process', 'workdir', 'status', 'cpus', 'pcpu', 'memory', 'pmem' ]

def run_nextflow_log(args):
    import subprocess
    command = ['nextflow', 'log'] + args

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        # print(result)
        raise RuntimeError(f"Error running nextflow log: {result.stderr}")
    return result.stdout

def get_all_run_names():
    import shutil
    from ..core.module_loader import get_environment_defaults
    env_defaults = get_environment_defaults()
    if shutil.which('nextflow') is None:
        nextflow_module = env_defaults.get('nextflow_module', 'nextflow')
        # Error out if nextflow is not found
        # in the error message use the nextflow_module variable if it is set
        if nextflow_module:
            raise RuntimeError(f"Nextflow command not found. Load the module '{nextflow_module}' to use Nextflow.")
        else:
            raise RuntimeError("Nextflow command not found. Load the Nextflow module to use Nextflow.")

    output = run_nextflow_log(['-q'])
    run_names = output.strip().split('\n')
    return run_names

def get_all_entries():
    run_names = get_all_run_names()
    results = []
    for run_name in run_names:
        args = ['-f', ','.join(FIELD_NAMES), run_name]
        try:
            output = run_nextflow_log(args)
            entries = parse_log_output(output)
            results.extend(entries)
        except RuntimeError as e:
            logging.warning(f"Skipping run '{run_name}' due to error: {e}")
            continue
    return results

def get_entries_with_string_in_name(search_string):
    run_names = get_all_run_names()
    results = []
    for run_name in run_names:
        args = [
            '-f', ','.join(FIELD_NAMES),
            '-F', f"name =~ '.*{search_string}.*'",
            run_name
        ]
        try:
            output = run_nextflow_log(args)
            entries = parse_log_output(output)
            results.extend(entries)
        except RuntimeError as e:
            logging.warning(f"Skipping run '{run_name}' due to error: {e}")
            continue
    return results

def get_entries_with_sample_names(sample_names):
    results = []
    for sample_name in sample_names:
        entries = get_entries_with_string_in_name(sample_name)
        results.extend(entries)
    return results

def get_entries_with_process_names(process_names):
    results = []
    for process_name in process_names:
        entries = get_entries_with_string_in_name(process_name)
        results.extend(entries)
    return results

def parse_log_output(output):
    lines = output.strip().split('\n')
    entries = []
    for line in lines:
        values = line.split('\t')
        if len(values) == 0 or all([v == '' for v in values]):
            continue
        entry = dict(zip(FIELD_NAMES, values))
        entries.append(entry)
    return entries
