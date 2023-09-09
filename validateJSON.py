import os
import json
import logging
import click
from collections import Counter
import jsonschema as js 
from jsonschema import Draft202012Validator
from pathlib import Path

from referencing import Registry, Resource
from referencing.exceptions import NoSuchResource

HERE = os.path.abspath(os.path.dirname(__file__))
SCHEMA_NAME = 'xrna-schema.json'
SCHEMAS = Path("./sections")
LOGGER = logging.getLogger(__name__)

def retrieve_from_filesystem(uri):
    expected_prefix = "http://localhost/"
    if not uri.startswith(expected_prefix):
        raise NoSuchResource(ref=uri)
    path = SCHEMAS / Path(uri.removeprefix(expected_prefix))
    contents = json.loads(path.read_text())
    return Resource.from_contents(contents)

# 1. Validates that there are no repeated nucleotide (residue) indices per RNA molecule.
# 2. Validates that there are no gaps betwwen nucleotide indices per RNA molecule.
# 3. Validates that no basePairs refer to an invalid nucleotide index (Nucleotide indices outside the sequence residueIndex range).
# 4. Validates that no labels or annotations refer to an invalid nucleotide index (nucleotide indices outside the sequence residueIndex range).
def validate_indices(data):
    rna_complexes = data.get("rnaComplexes", [])
    for rna_complex in rna_complexes:
        rna_molecules = rna_complex.get("rnaMolecules", [])
        for rna_molecule in rna_molecules:
            sequence = rna_molecule.get("sequence", [])
            seen_residue_indices = set()
            for nucleotide in sequence:
                if "residueIndex" in nucleotide:
                    residue_index = nucleotide["residueIndex"]
                    if residue_index in seen_residue_indices:
                        yield js.ValidationError("Within a single RNA molecule, residueIndex values must not be repeated.")
                    else:
                        seen_residue_indices.add(residue_index)
            basePairs = rna_molecule.get("basePairs", [])
            for basePair in basePairs:
                if "residueIndex1" in basePair and "residueIndex2" in basePair:
                    residue_index_1 = basePair["residueIndex1"]
                    residue_index_2 = basePair["residueIndex2"]
                    absolute_distance = abs(residue_index_1 - residue_index_2)
                    if absolute_distance == 0:
                        yield js.ValidationError("A nucleotide may not base pair with itself.")
                    if absolute_distance == 1:
                        yield js.ValidationError("A nucleotide may not base pair with adjacent nucleotides.")
                    if not residue_index_1 in seen_residue_indices:
                        yield js.ValidationError("residueIndex1 values within basePair objects must reference a valid nucleotide residueIndex.")
                    if not residue_index_2 in seen_residue_indices:
                        yield js.ValidationError("residueIndex2 values within basePair objects must reference a valid nucleotide residueIndex.")
            if len(sequence) > 0:
                minimum_residue_id = min(seen_residue_indices, key=lambda residue_index: int(residue_index))
                maximum_residue_id = max(seen_residue_indices, key=lambda residue_index: int(residue_index))
                if maximum_residue_id - minimum_residue_id != len(sequence) - 1:
                    yield js.ValidationError("RNA molecule sequences must not leave gaps between residue indices")
            labels = rna_molecule.get("labels", [])
            for label in labels:
                if "residueIndex" in label:
                    residue_index = label["residueIndex"]
                    if not residue_index in seen_residue_indices:
                        yield js.ValidationError("residueIndex values within labels must reference a valid nucleotide residueIndex.")

# Validates that all referenced classes are present within the top-level "classes" object
def validate_classes(data):
    seen_classes = set()
    rna_complexes = data.get("rnaComplexes", [])
    for rna_complex in rna_complexes:
        rna_molecules = rna_complex.get("rnaMolecules", [])
        for rna_molecule in rna_molecules:
            sequence = rna_molecule.get("sequence", [])
            for residue in sequence:
                classes = residue.get("classes", [])
                for class_ in classes:
                    seen_classes.add(class_)
            labels = rna_molecule.get("labels", [])
            for class_ in labels:
                if "labelLine" in class_:
                    label_line = class_["labelLine"]
                    classes = label_line.get("classes", [])
                    for class_ in classes:
                        seen_classes.add(class_)
                if "labelContent" in class_:
                    label_content = class_["labelContent"]
                    classes = label_content.get("classes", [])
                    for class_ in classes:
                        seen_classes.add(class_)
            basePairs = rna_molecule.get("basePairs", [])
            for basePair in basePairs:
                classes = basePair.get("classes", [])
                for class_ in classes:
                    seen_classes.add(class_)
            classes_for_sequence = rna_molecule.get("classesForSequence", [])
            for class_ in classes_for_sequence:
                seen_classes.add(class_)
            classes_for_labels = rna_molecule.get("classesForLabels", [])
            for class_ in classes_for_labels:
                seen_classes.add(class_)
            classes_for_base_pairs = rna_molecule.get("classesForBasePairs", [])
            for class_ in classes_for_base_pairs:
                seen_classes.add(class_)
    classes = set()
    for class_ in data.get("classes", []):
        if "name" in class_:
            classes.add(class_["name"])
    for seen_class in seen_classes:
        if not seen_class in classes:
            yield js.ValidationError("All referenced classes must be defined within the top-level \"classes\" array.")

def validate(data, schema_path):
    registry = Registry(retrieve=retrieve_from_filesystem)

    with open(schema_path, "r") as raw:
        schema = json.load(raw)
    validator = Draft202012Validator(
        schema,
        registry=registry
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

    for error in validate_classes(data):
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
def main(filename, schema=None):
    with open(filename, 'r') as raw:
        data = json.load(raw)

    validate(data, schema)

if __name__ == '__main__':

    logging.basicConfig(
        format='%(levelname)s: %(message)s',
        level=logging.WARNING,
    )

    main()