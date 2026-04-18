#Possible implementations:
# Perform cross-validation for file inputs
# Validate exon ordering per transcript
# Validate read alignment plausibility

import sys
import re
from pathlib import Path

# file response
def error(msg):
    raise ValueError(f"[ERROR] {msg}")

def warn(msg):
    print(f"[WARNING] {msg}")

def info(msg):
    print(f"[INFO] {msg}")


# detect file type by extension
def detect_file_type(filepath: Path):
    ext = filepath.suffix.lower()

    if ext == ".gtf":
        return "gtf"
    if ext in {".fa", ".fasta"}:
        return "fasta"
    if ext in {".fq", ".fastq"}:
        return "fastq"

    error(f"Unsupported file extension: {ext}")


# validation for fasta
def validate_fasta(filepath: Path):
    info(f"Validating FASTA: {filepath}")

    headers = []
    sequences = []

    with open(filepath) as f:
        current_seq = []
        current_header = None

        for line_no, line in enumerate(f, 1):
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_header:
                    sequences.append("".join(current_seq))
                    current_seq = []

                current_header = line[1:].strip()
                if not current_header:
                    # non-fatal error
                    error(f"Empty FASTA header at line {line_no}")

                headers.append(current_header)
            else:
                # possible errors
                if current_header is None:
                    error(f"Sequence before first header at line {line_no}")

                if not re.fullmatch(r"[ACGTNacgtn]+", line):
                    error(f"Invalid nucleotide characters at line {line_no}: {line}")

                current_seq.append(line)

        if current_header:
            sequences.append("".join(current_seq))

    if not headers:
        error("No FASTA headers found")

    if len(headers) != len(sequences):
        error("Mismatch between FASTA headers and sequences")

    info(f"FASTA OK: {len(headers)} sequences")
    return set(headers)


# file validation for fastq
def validate_fastq(filepath: Path):
    info(f"Validating FASTQ: {filepath}")

    headers = set()

    with open(filepath) as f:
        line_no = 0

        while True:
            h = f.readline()
            if not h:
                break

            s = f.readline()
            p = f.readline()
            q = f.readline()

            line_no += 4

            if not (s and p and q):
                error(f"Incomplete FASTQ record ending at line {line_no}")

            h = h.strip()
            s = s.strip()
            p = p.strip()
            q = q.strip()

            if not h.startswith("@"):
                error(f"FASTQ header does not start with @ at line {line_no - 3}")

            if not p.startswith("+"):
                error(f"FASTQ separator line does not start with + at line {line_no - 1}")

            header_id = h[1:].split()[0]
            headers.add(header_id)

            if not re.fullmatch(r"[ACGTNacgtn]+", s):
                error(f"Invalid FASTQ sequence characters at line {line_no - 2}")

            if len(s) != len(q):
                error(f"Sequence and quality length mismatch at line {line_no - 2}")

    if not headers:
        error("No FASTQ records found")

    info(f"FASTQ OK: {len(headers)} reads")
    return headers


# file validation for gtf
def parse_gtf_attributes(attr_field: str):
    attrs = {}
    for part in attr_field.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            error(f"Malformed GTF attribute: {part}")
        key, value = part.split(" ", 1)
        attrs[key] = value.strip('"')
    return attrs


def validate_gtf(filepath: Path):
    info(f"Validating GTF: {filepath}")

    transcript_ids = set()
    gene_ids = set()

    with open(filepath) as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                error(f"GTF line {line_no} does not have 9 fields")

            seqname, source, feature, start, end, score, strand, frame, attrs = fields

            if not start.isdigit() or not end.isdigit():
                error(f"Non-integer coordinates at line {line_no}")

            start = int(start)
            end = int(end)

            if start < 1 or end < start:
                error(f"Invalid coordinate range at line {line_no}: {start}-{end}")

            if strand not in {"+", "-", "."}:
                error(f"Invalid strand at line {line_no}: {strand}")

            attr_dict = parse_gtf_attributes(attrs)

            if "gene_id" not in attr_dict:
                error(f"Missing gene_id at line {line_no}")

            gene_ids.add(attr_dict["gene_id"])

            if "transcript_id" in attr_dict:
                transcript_ids.add(attr_dict["transcript_id"])

    if not gene_ids:
        error("No gene_id values found in GTF")

    info(f"GTF OK: {len(gene_ids)} genes, {len(transcript_ids)} transcripts")
    return gene_ids, transcript_ids


# perform cross validation for input files
def cross_validate(gtf_genes, gtf_transcripts, fasta_headers, fastq_headers):
    info("Cross-validating identifiers")

    missing_in_fasta = gtf_transcripts - fasta_headers
    if missing_in_fasta:
        warn(f"{len(missing_in_fasta)} GTF transcripts missing in FASTA")

    extra_in_fasta = fasta_headers - gtf_transcripts
    if extra_in_fasta:
        warn(f"{len(extra_in_fasta)} FASTA sequences not present in GTF")

    if fastq_headers:
        overlap = fasta_headers & fastq_headers
        if not overlap:
            warn("No overlap between FASTA headers and FASTQ read IDs")

    info("Cross-validation complete")


# MAIN FUNCTION
def validate_files(filepaths):
    gtf_genes = set()
    gtf_transcripts = set()
    fasta_headers = set()
    fastq_headers = set()

    for fp in filepaths:
        fp = Path(fp)
        ftype = detect_file_type(fp)

        if ftype == "gtf":
            gtf_genes, gtf_transcripts = validate_gtf(fp)

        elif ftype == "fasta":
            fasta_headers = validate_fasta(fp)

        elif ftype == "fastq":
            fastq_headers = validate_fastq(fp)

    if gtf_transcripts and fasta_headers:
        cross_validate(gtf_genes, gtf_transcripts, fasta_headers, fastq_headers)

    info("All validations finished successfully")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python validate_omics.py <file1> <file2> <file3>")
        sys.exit(1)

    try:
        validate_files(sys.argv[1:])
    except ValueError as e:
        print(e)
        sys.exit(2)


