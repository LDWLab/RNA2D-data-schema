import os
import json
import logging
import click
from collections import Counter
import jsonschema as js

HERE = os.path.abspath(os.path.dirname(__file__))
SECTIONS = os.path.join(HERE, 'sections')
SCHEMA_NAME = 'xrna-schema.json'
LOGGER = logging.getLogger(__name__)

# 1. Validates that there are no repeated nucleotide indices per RNA molecule.
# 2. Validates that there are no gaps betwwen nucleotide indices per RNA molecule.
# 3. Validates that no basePairs refer to an invalid nucleotide index (Nucleotide indices outside the sequence resId range).
# 4. Validates that no labels or annotations refer to an invalid nucleotide index (nucleotide indices outside the sequence resId range).
def validate_indices(data):
    rna_complexes = data.get("rna complexes", [])
    for rna_complex in rna_complexes:
        rna_molecules = rna_complex.get("rna molecules", [])
        for rna_molecule in rna_molecules:
            sequence = rna_molecule.get("sequence", [])
            seen_residue_ids = set()
            for nucleotide in sequence:
                if "resId" in nucleotide:
                    res_id = nucleotide["resId"]
                    if res_id in seen_residue_ids:
                        yield js.ValidationError("Within a single RNA molecule, resId values must not be repeated.")
                    else:
                        seen_residue_ids.add(res_id)
            basePairs = rna_molecule.get("basePairs", [])
            for basePair in basePairs:
                if "resId1" in basePair and "resId2" in basePair:
                    res_id_1 = basePair["resId1"]
                    res_id_2 = basePair["resId2"]
                    absolute_distance = abs(res_id_1 - res_id_2)
                    if absolute_distance == 0:
                        yield js.ValidationError("A nucleotide may not base pair with itself.")
                    if absolute_distance == 1:
                        yield js.ValidationError("A nucleotide may not base pair with adjacent nucleotides.")
                    if not res_id_1 in seen_residue_ids:
                        yield js.ValidationError("resId1 values within basePair objects must reference a valid nucleotide resId.")
                    if not res_id_2 in seen_residue_ids:
                        yield js.ValidationError("resId2 values within basePair objects must reference a valid nucleotide resId.")
            minimum_residue_id = min(seen_residue_ids, key=lambda x: int(x))
            maximum_residue_id = max(seen_residue_ids, key=lambda x: int(x))
            labels_and_annotations = rna_molecule.get("labelsAndAnnotations", [])
            if len(sequence) > 0 and (maximum_residue_id - minimum_residue_id) != len(sequence) - 1:
                yield js.ValidationError("RNA molecule sequences must not leave gaps between residue ids")
            for label_or_annotation in labels_and_annotations:
                if "resId" in label_or_annotation:
                    res_id = label_or_annotation["resId"]
                    if not res_id in seen_residue_ids:
                        yield js.ValidationError("ResId values within labels and annotations must reference a valid nucleotide resId.")
            

def validate(data, schema_path, sections_path):

    with open(schema_path, 'r') as raw:
        schema = json.load(raw)

    base = 'file://%s/' % sections_path
    validator = js.Draft4Validator(
        schema,
        format_checker=js.FormatChecker(),
        resolver=js.RefResolver(base, None),
    )

    found = False
    counts = Counter()
    for error in validator.iter_errors(data):
        counts[error.validator] += 1
        found = True
        LOGGER.error(error.message)

    for error in validate_indices(data):
        counts[error.validator] += 1
        found = True
        LOGGER.error(error.message)

    if found:
        summary = ', '.join('%s: %s' % (k, v) for k, v in counts.items())
        raise click.ClickException("Validation failed: %s" % summary)
    else:
        print("Validation succeeded.")

@click.command()
@click.argument('filename')
@click.option('--schema', default=SCHEMA_NAME,
              help='Filename of the schema to use')
@click.option('--sections', default=SECTIONS,
              help='Directory where schema parts are kept')
def main(filename, schema=None, sections=None):
    with open(filename, 'r') as raw:
        data = json.load(raw)

    validate(data, schema, os.path.abspath(sections))

if __name__ == '__main__':

    logging.basicConfig(
        format='%(levelname)s: %(message)s',
        level=logging.WARNING,
    )

    main()