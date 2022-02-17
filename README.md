# RNA2D-data-schema

RNA2D-data-schema defines a standard JSON Schema to facilitate export of RNA secondary structures (2D) in [specific layouts](https://www.nature.com/articles/s41467-021-23555-5). 

For example, work is underway to use the schema to export RNA 2D structures from [R2DT](https://github.com/rnacentral/r2dt) as JSON files that can be imported into another software for manual editing. 

:warning: This is work in progress and the schema is subject to change.

### Features

The following information can currently be stored:

- nucleotide sequence
- (x, y) coordinates of each nucleotides
- labels and classes can be optionally attached to each nucleotide

A script can be used to validate a JSON file and ensure that the data is formatted correctly.

### Get involved

We invite anyone interested in using the schema in your work to get involved by creating [an issue](https://github.com/LDWLab/RNA2D-data-schema/issues).
