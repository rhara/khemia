#!/usr/bin/env python

from rdkit import Chem
from khemia.io import MolReader, MolWriter
from khemia.timer import Timer
import argparse

def main():
    parser = argparse.ArgumentParser(description='Substructure search molecule')
    parser.add_argument('--in', '-i', type=str, required=True, dest='iname', help='Input file')
    parser.add_argument('--out', '-o', type=str, required=True, dest='oname', help='Output file')
    parser.add_argument('--pattern', '-p', type=str, required=True, help='Pattern')
    args = parser.parse_args()

    pat = Chem.MolFromSmarts(args.pattern)

    reader = MolReader(args.iname)
    writer = MolWriter(args.oname)

    count = 0
    match_count = 0
    write_count = 0
    timer = Timer(dt=2.0)
    for mol in reader:
        count += 1
        match = mol.HasSubstructMatch(pat)
        if match:
            match_count += 1
            success = writer.write(mol)
            if success:
                write_count += 1
        ret = timer.check()
        if ret:
            lap = timer.elapsed()
            print(f'{lap:7.2f}s {count}/{match_count}/{write_count}')
    lap = timer.elapsed()
    print(f'{lap:7.2f}s {count}/{match_count}/{write_count}')
    reader.close()
    writer.close()

if __name__ == '__main__':
    main()
