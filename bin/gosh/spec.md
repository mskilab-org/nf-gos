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

