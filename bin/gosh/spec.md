# gOSh

A CLI for the nf-gOS pipeline

## Features

### Run
- Sensible mode-specific defaults for parameters (heme, tumor-only, paired)
  - these should be relegated to profiles (likewise for the hpc) in nf-gOS
    which are then passed to the --profile parameter via run command
    - Set up common parameters in profile
       - CPLEX
       - HPC config
       - Dependent files (PON, reference, etc.)
- Resume/overwrite from a process (by deleting the working directories)
- Resource sequestering
- [ ] Auto-resume on OOM
    - [ ] detect if all expected results are present, auto-resume if not
    - [ ] detect if the run was interrupted by user, don't auto-resume in this case
- More informative email notifications
- Help with samplesheet construction
- Stream samplesheet from remote
- refactor:
  - move the module load environment spec to its own script
  - save the debug log on first command then re-use instead of recalling nextflow log

### Debug

- Run information (resource usage, run parameters)
- Working directory by sample ID(s) x process(es)
- AI parsing of nextflow.log files containing errors
- Results directory by sample ID(s) x process(es)
- Export output to .csv for loading in R
- Samplesheet validation

### Help

- Help for each command
- Information on common pipeline parameters

### Skilift

- Skilift current results into gOS

### Architecture

 1 Approach: Mixed OOP & Functional

 • Use OOP for core components that maintain state or represent complex entities (e.g., a
   NextflowRunner class, SampleSheet class)
 • Use functional approach for utilities and pure operations (e.g., log parsing, validation functions)

 2 Program Structure


 gosh/
 ├── cli/
 │   ├── run.py      # Run command group
 │   ├── debug.py    # Debug command group
 │   ├── help.py     # Help command group
 │   ├── skilift.py  # Skilift command group
 ├── core/
 │   ├── nextflow.py    # NextflowRunner class
 │   ├── samplesheet.py # SampleSheet handling
 │   ├── config.py      # Configuration management
 │   ├── resources.py   # Resource management
 ├── utils/
 │   ├── log_parser.py  # Nextflow log parsing
 │   ├── validators.py  # Input validation
 │   ├── ai_helper.py   # AI log parsing
 │   ├── email.py       # Email notifications
 ├── main.py        # CLI entry point


### Prompts
- [ ] Working directory by sample ID(s) x process(es)

feature: create a script `nextflow_log.py` that will serve as a wrapper on the `nextflow log` command

- put this script in `gosh/core/nextflow_log.py`
- there should be a general method to run `nextflow log` command (with passed arguments)
- that general method should be used in other methods which do specific things:
  1. get a list of all run names (`nextflow log -q`)
  2. get a list of all entries across all runs that contain a certain string in the name
    - to do this you will have to first get all run names and then iterate over them using the command: `nextflow log -f {list_of_all_fields} -F 'name ~= /.*{string}.*/' {run_name}`
    - then append the results to a list
  3. a method that uses #2 to get a list of entries across all runs that contain a any of the sample names passed to it (in a list) in the name
  4. a method that uses #2 to get a list of entries across all runs that contain a any of the process names passed to it (in a list) in the name

list of all fields:
```
  attempt complete container cpus disk duration env error_action exit hash
hostname inv_ctxt log memory module name native_id pcpu peak_rss peak_vmem pmem
process queue rchar read_bytes realtime rss scratch script start status stderr
stdout submit syscr syscw tag task_id time vmem vol_ctxt wchar workdir
write_bytes
```

- [o] ai ask command for asking questions
feature: add some new methods and prompts to the ai_helper.py for handling help related questions

- add a system prompt for answering questions about the nf-gOS pipeline, the gOSh CLI, and nextflow bioinformatics pipelines
  - it will include the contents of a `help_context.txt` file that includes some information about the pipeline and the CLI
  - it will include a section detailing how to use the `params.json` file if it has been appended to the question and is relevant to the user's question
- add a method that uses this prompt to answer user questions
  - the method will take a user inputted question as string input
  - it will look for a `params.json` file in the current directory, if it exists, it will append it to the question


### Tasks

- [x] Run command
  - [x] set up profiles/configs
    - [x] hpc
      - nyu
      - nygc
      - cplex
    - [x] alignment
      - gpu vs. cpu
  - [x] predict correct set of profiles input params using env information and mode
      - [x] echo $HOSTNAME to get the host -> -profile nyu/nygc if recognized
      - [x] check if normals are in the samplesheet to determine mode (tumor-only vs. paired)
  - [x] Resume/overwrite from a process (by deleting the working directories)
- [x] Debug command
  - [x] Working directory by sample ID(s) x process(es)
  - [x] Run information (resource usage, run parameters)
  - [x] ai ask command for asking questions
  - [x] ai parsing of nextflow (gosh debug eye [.nextflow.log])
    - [x] copy some real error .nextflow.log files for testing (create a test dir)
    - [x] write the ai_helper.py (use personal api key for now)
      - [x] parser for extracting error messages from nextflow.log
      - [x] base method for sending the query/prompt to the ai
      - [x] prompt for interpreting the error message
      - [x] prompt for advising the user on how to fix the error
      - [x] method to chain the prompts into a single response
      - [x] send the ai response to debug
- [x] Skilift command

