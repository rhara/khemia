#!/usr/bin/env python

from khemia import MolReader
import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert molecule file')
    parser.add_argument('--in', '-i', type=str, required=True, dest='iname', help='Input file')
    args = parser.parse_args()

    reader = MolReader(args.iname)

    count = 0
    error = 0
    for mol in reader:
        if mol is None:
            error += 1
        else:
            count += 1
    reader.close()

    print('mol: %d, error: %d' % (count, error))

if __name__ == '__main__':
    main()
