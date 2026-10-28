#/usr/bin/python3

from __future__ import annotations
import argparse
import random
from pathlib import Path

FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas"}
HEADER_INSERT= ["SNP","INS","DEL","ALL"]
ROOT = Path(__file__).resolve().parent
if ROOT.name == "scripts":
    ROOT = ROOT.parent



def project_path(path: str) -> Path:
    """Interpret relative paths from the project root, not the current cwd."""
    p = Path(path).expanduser()
    return (p if p.is_absolute() else ROOT / p).resolve()


def fasta_records(path: Path):
    """Read multi-record FASTA; yield (header, uppercased sequence)."""
    header = None
    parts = []
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None:
                    if not parts:
                        raise ValueError(f"Empty FASTA record: {path} >{header}")
                    yield header, "".join(parts)
                header, parts = line[1:].strip(), []
                if not header:
                    raise ValueError(f"Empty FASTA header: {path}")
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                parts.append("".join(line.split()).upper())
    if header is not None:
        if not parts:
            raise ValueError(f"Empty FASTA record: {path} >{header}")
        yield header, "".join(parts)


def pick_reference(ref_path: Path, rng: random.Random):
    """Uniformly sample a FASTA record across files and records."""
    if ref_path.is_file():
        files = [ref_path]
    elif ref_path.is_dir():
        files = sorted(p for p in ref_path.rglob("*")
                       if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES)
    else:
        raise ValueError(f"Reference not found: {ref_path}")
    if not files:
        raise ValueError(f"No FASTA files found: {ref_path}")

    selected = None
    seen = 0
    for file in files:
        for header, sequence in fasta_records(file):
            seen += 1
            if rng.randrange(seen) == 0:
                selected = (file, header, sequence)
    if selected is None:
        raise ValueError("No FASTA records found")
    return selected

def header_partition(header):
    serial_num,name=header.split(" ",1)
    split_header=[serial_num,name]
    return split_header

def write_fasta(path, header, sequence):
    with path.open("w", encoding="utf-8") as out:
        for s in HEADER_INSERT:
            serial_num,name=header_partition(header)
            out.write(f">{serial_num}|{s} {name}\n")
            for i in range(0, len(sequence), 80):
                out.write(sequence[i:i + 80] + "\n")




def main():
    parser = argparse.ArgumentParser(
        description="From one random reference genome, output SNP/DEL/INS/ALL mutants."
    )
    parser.add_argument("--ref", "--input", dest="ref", default="data/ref/ref.fasta",
                        help="FASTA file (default: data/ref/ref.fasta)")

    args = parser.parse_args()

    ref_path = project_path(args.ref)
    selected_output = project_path("data/selected/selected.fasta")
    if not ref_path.exists(): 
        parser.error("input must be an existing file")
    if ref_path.resolve() == selected_output.resolve():
        parser.error("reference FASTA file cannot be overwirte")
    if not ref_path.suffix.lower() in FASTA_SUFFIXES :
        parser.error("cannot open the file")
    seed= random.randrange(30000)
    print(seed)
    rng = random.Random(seed)
    source_file,header,sequence=pick_reference(ref_path,rng)
    write_fasta(selected_output,header,sequence)


if __name__ == "__main__":
    main()
