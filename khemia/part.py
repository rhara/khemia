import numpy as np
from rdkit import Chem


def submolFromPart(mol, part, n):
    submol = Chem.RWMol()
    mapping = {}
    for idx in np.where(part == n)[0]:
        idx = idx.item()
        _idx = submol.AddAtom(mol.GetAtomWithIdx(idx))
        mapping[idx] = _idx
    for i in mapping:
        for j in mapping:
            if j <= i:
                continue
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond:
                _type = bond.GetBondType()
                submol.AddBond(mapping[i], mapping[j], _type)
    return submol.GetMol()


def getRingLinkerPart(mol):
    def trav_infrag(atom, ring=False):
        def rec(atom):
            if atom.GetIdx() in idxs:
                return
            idxs.add(atom.GetIdx())
            for neighbor in atom.GetNeighbors():
                if ring:
                    if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).IsInRing():
                        rec(neighbor)
                else:
                    if not neighbor.IsInRing():
                        rec(neighbor)
        idxs = set()
        rec(atom)
        return idxs

    part = np.zeros(mol.GetNumAtoms(), dtype=np.int16)

    i, n, N = 0, 1, mol.GetNumAtoms()
    while True:
        if N <= i:
            break
        if 0 < part[i]:
            i += 1
            continue
        atom = mol.GetAtomWithIdx(i)
        for j in trav_infrag(atom, ring=atom.IsInRing()):
            part[j] = n
        n += 1
        i += 1

    return part

