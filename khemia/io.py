from rdkit import Chem
import os, io, gzip, bz2, time


def MolReader(iname):
    if iname.endswith('.smi.gz') or iname.endswith('.smi'):
        return SmilesReader(iname)
    elif iname.endswith('.sdf.gz') or iname.endswith('.sdf'):
        return SDFReader(iname)
    else:
        raise Exception('Unknown reader format')


def MolWriter(oname):
    if oname.endswith('.smi.gz') or oname.endswith('.smi'):
        return SmilesWriter(oname)
    elif oname.endswith('.sdf.gz') or oname.endswith('.sdf'):
        return SDFWriter(oname)
    else:
        raise Exception('Unknown reader format')


class BaseReader:
    def __init__(self, iname):
        self.iname = iname

    def __iter__(self):
        pass

    def close(self):
        try:
            self.ofs.close()
        except:
            pass


class SmilesReader(BaseReader):
    def __init__(self, iname):
        super(SmilesReader, self).__init__(iname)
        fopen = gzip.open if iname.endswith('.gz') else open
        self.ifs = fopen(iname, 'rt')

    def __iter__(self):
        for line in self.ifs.readlines():
            line = line.rstrip()
            it = line.split(' ', maxsplit=1)
            smiles = it[0]
            supp = it[1] if len(it) == 2 else ''
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol.SetProp('_Name', supp)
                yield mol


class SDFReader(BaseReader):
    def __init__(self, iname):
        super(SDFReader, self).__init__(iname)
        if iname.endswith('.gz'):
            self.ifs = Chem.ForwardSDMolSupplier(gzip.open(iname, 'rb'))
        else:
            self.ifs = Chem.SDMolSupplier(iname)

    def __iter__(self):
        for mol in self.ifs:
            yield mol


class BaseWriter:
    def __init__(self, oname):
        self.oname = oname
        self.ofs = None

    def write(self, mol):
        try:
            self.ofs.write(mol)
            return True
        except:
            return False

    def close(self):
        try:
            self.ofs.flush()
        except:
            pass
        try:
            self.ofs.close()
        except:
            pass


class SmilesWriter(BaseWriter):
    def __init__(self, oname):
        super(SmilesWriter, self).__init__(oname)
        if oname.endswith('.gz'):
            self.ofs = Chem.SmilesWriter(gzip.open(oname, 'wt'))
        else:
            self.ofs = Chem.SmilesWriter(oname)


class SDFWriter(BaseWriter):
    def __init__(self, oname):
        super(SDFWriter, self).__init__(oname)
        if oname.endswith('.gz'):
            self.ofs = Chem.SDWriter(gzip.open(oname, 'wt'))
        else:
            self.ofs = Chem.SDWriter(oname)

    def write(self, mol, confId=None):
        try:
            if confId is None:
                self.ofs.write(mol)
            else:
                self.ofs.write(mol, confId=confId)
            return True
        except:
            return False


