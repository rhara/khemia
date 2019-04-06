#!/usr/bin/env python

from khemia import MolReader, MolWriter
import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert molecule file')
    parser.add_argument('--in', '-i', type=str, required=True, dest='iname', help='Input file')
    parser.add_argument('--out', '-o', type=str, required=True, dest='oname', help='Output file')
    parser.add_argument('--tag', '-t', type=str, help='Tag to set title')
    parser.add_argument('-n', type=int, default=None, help='Max molecules in output')
    args = parser.parse_args()

    reader = MolReader(args.iname)
    writer = MolWriter(args.oname)

    count = 0
    for mol in reader:
        if args.tag is not None:
            name = mol.GetProp(args.tag)
            mol.SetProp('_Name', name)
        success = writer.write(mol)
        if success:
            count += 1
        if args.n is not None and args.n <= count:
            break
    print(count)
    reader.close()
    writer.close()

if __name__ == '__main__':
    main()
