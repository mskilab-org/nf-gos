# mskilab-org/nf-jabba: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of mskilab-org/nf-jabba, created with the [nf-core](https://nf-co.re/) template.

### `Added`
- [JaBbA container optimization handoff](docs/jabba-container-optimization-handoff.md) covering package migration, runtime-patch removal, solver parameter precedence, and acceptance checks.

### `Changed`
- Tighten JaBbA's CPLEX formulation by deriving big-M from observed copy number and explicitly bounding binary indicator and loose-end variables; tighten the default relative MIP gap from 0.1% to 0.01%. `--epgap_jabba` remains configurable for runtime-sensitive runs.
- Run non-integer balance through a staged runtime wrapper with adaptive big-M, explicit binary/loose-end upper bounds, preserved fractional-CN domains, allocation-matched CPLEX threads, and configurable `--mipemphasis_non_integer_balance` (default 0). No container rebuild required.

### `Fixed`
- Preserve numeric zero in scalar parameter channels so the default `mipemphasis_non_integer_balance = 0` schedules non-integer balance after JaBbA instead of silently closing its input. Select one tumor per patient for new and supplied non-integer balance outputs, excluding matched-normal rows from patient-keyed joins.
- Add a scalar-channel regression workflow: `nextflow run tests/value_channels.nf -lib lib` (covers numeric zero, booleans, and absent values).
- Route Esvee somatic and unfiltered VCF chimera tagging through native VF/site-QUAL support fields and VF matched-normal evidence; preserve advisory tags, existing filters, output channels, and supplied-output bypass. Keep chimera publication/resources separate from the Esvee caller.
- Preserve paired tumor and reference Cobalt filenames when staging PURPLE inputs.
- Allow mixed JaBbA inputs with missing `j_supp` values to survive upstream channel joins and pass through sequence-name coercion.
- Match the non-integer balance stub output filename to the declared `non_integer.balanced.gg.rds` output.

### `Dependencies`

### `Deprecated`
