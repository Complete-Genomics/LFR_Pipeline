#!/usr/bin/env python3
"""Split one primary-contig FASTA into four named, non-overlapping quarters."""
import argparse


def records(path):
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.out, "w") as handle:
        for name, sequence in records(args.input):
            for quarter in range(4):
                start = quarter * len(sequence) // 4
                end = (quarter + 1) * len(sequence) // 4
                handle.write(">{0}||Q{1}\n{2}\n".format(name, quarter, sequence[start:end]))


if __name__ == "__main__":
    main()
