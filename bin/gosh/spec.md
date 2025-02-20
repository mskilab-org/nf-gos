# gOSh

A CLI for the nf-gOS pipeline

## Features

### Run
- Sensible mode-specific defaults for parameters (heme, tumor-only, paired)
  - these should be relegated to profiles (likewise for the hpc) in nf-gOS
    which are then passed to the --profile parameter
- Resume/overwrite from a process (by deleting the working directories)
- Resource sequestering
- Auto-resume on OOM
- More informative email notifications
- Exposure of common parameters
   - CPLEX
   - HPC config
   - Dependent files (PON, reference, etc.)
- Stream samplesheet from remote
- Help with samplesheet construction

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


 3 Key Abstractions

 • NextflowRunner: Core class to handle nextflow execution
    • Manages process resumption
    • Handles resource allocation
    • Monitors execution
    • Auto-resume on OOM
 • SampleSheet: Class for samplesheet operations
    • Validation
    • Remote streaming
    • Construction helpers
 • ConfigManager: Handle different execution modes
    • Mode-specific defaults (heme, tumor-only, paired)
    • Resource configurations
    • Parameter management
 • LogAnalyzer: Parse and analyze nextflow logs
    • Error detection
    • Resource usage analysis
    • AI-powered error interpretation

 4 Testing Strategy

 • Unit Tests
    • Mock nextflow execution
    • Test configuration generation
    • Validate samplesheet parsing
 • Integration Tests
    • Test with minimal nextflow pipelines
    • Verify command execution
    • Test resource management
 • Fixtures
    • Sample nextflow logs
    • Example samplesheets
    • Mock configurations
 • Test Categories
    • Command-line interface tests
    • Configuration management
    • Log parsing
    • Resource management
    • Error handling

 5 Additional Considerations

 • Use dependency injection for better testability
 • Implement proper logging throughout
 • Use type hints for better code maintainability
 • Consider using pydantic for config/data validation
 • Use click.testing for CLI tests

### Tasks

- ai parsing of nextflow (gosh debug eye [.nextflow.log])
  - [x] copy some real error .nextflow.log files for testing (create a test dir)
  - write the ai_helper.py (use personal api key for now)
    - [x] parser for extracting error messages from nextflow.log
    - [x] base method for sending the query/prompt to the ai
    - [x] prompt for interpreting the error message
    - [x] prompt for advising the user on how to fix the error
    - method to chain the prompts into a single response
    - send the ai response to debug
  - pretty printing the response (via click?) in debug
