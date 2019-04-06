#!/usr/bin/env python

from khemia.io import MolReader, MolWriter
import re, argparse

class GE:
    def __init__(self, a):
        self.a = a
    def __call__(self, v):
        return self.a <= v

class LE:
    def __init__(self, a):
        self.a = a
    def __call__(self, v):
        return v <= self.a

class EQ:
    def __init__(self, a):
        self.a = a
    def __call__(self, v):
        return v == self.a

class BTWN:
    def __init__(self, a, b):
        self.a = a
        self.b = b
    def __call__(self, v):
        return self.a <= v <= self.b

def main():
    parser = argparse.ArgumentParser(description='Convert molecule file')
    parser.add_argument('--in', '-i', type=str, dest='infname', help='Input file')
    parser.add_argument('--out', '-o', type=str, help='Output file')
    parser.add_argument('--numbers', '-n', type=str, help='List format')
    args = parser.parse_args()

    pat_ab = re.compile('^([0-9]+)-([0-9]+)$')
    pat_a_ = re.compile('^([0-9]*)-$')
    pat__b = re.compile('^-([0-9]*)$')
    pat_a  = re.compile('^([0-9]+)$')
    upper = 0
    lower = float('inf')
    evals = []
    for it in args.numbers.split(','):
        m = pat_ab.match(it)
        if m:
            a = int(m.group(1))
            b = int(m.group(2))
            lower = min(lower, a)
            upper = max(upper, b)
            evals.append(BTWN(a, b))
        m = pat_a_.match(it)
        if m:
            a = int(m.group(1))
            lower = min(lower, a)
            evals.append(GE(a))
            continue
        m = pat__b.match(it)
        if m:
            b = int(m.group(1))
            upper = max(upper, b)
            evals.append(LE(b))
            continue
        m = pat_a.match(it)
        if m:
            a = int(m.group(1))
            lower = min(lower, a)
            upper = max(upper, a)
            evals.append(EQ(a))
            continue

    g_le = LE(lower-1)
    g_ge = GE(upper+1)

    reader = MolReader(args.infname)
    writer = MolWriter(args.out)

    count = 0
    written = 0
    for mol in reader:
        if mol is None:
            continue
        else:
            count += 1
        if lower != float('inf') and g_le(count):
            continue
        if g_ge(count) and upper != 0:
            break
        for ev in evals:
            if ev(count):
                writer.write(mol)
                written += 1

    reader.close()
    writer.close()

    print('written: %d' % written)

if __name__ == '__main__':
    main()
