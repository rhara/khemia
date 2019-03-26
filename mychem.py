from rdkit import Chem
import os, io, gzip, time

def getExtension(fname):
    dirname, bname = os.path.split(fname)
    i = bname.rindex('.')
    ext = bname[i+1:]
    if ext in ['bz2', 'gz', 'xz', 'zip', 'zst']:
        try:
            i = bname[:i].rindex('.')
            ext = bname[i+1:]
        except:
            pass
    return ext


class MolReader:
    def __init__(self, fname):
        ext = getExtension(fname)
        self.mode = 'mol'
        if ext == 'sdf':
            self.ifs = Chem.SDMolSupplier(fname)
        elif ext == 'sdf.gz':
            self.istrm = gzip.open(fname, 'rb')
            self.ifs = Chem.ForwardSDMolSupplier(self.istrm)
        elif ext in ['ism', 'smi', 'smiles']:
            self.mode = 'line'
            self.ifs = open(fname, 'rt')
        elif ext in ['ism.gz', 'smi.gz', 'smiles.gz']:
            self.mode = 'line'
            self.ifs = gzip.open(fname, 'rt')
        else:
            raise Exception(f'Unknown input format "{ext}"')

    def __iter__(self):
        for obj in self.ifs:
            if self.mode == 'line':
                line = obj.strip()
                it = line.split(maxsplit=1)
                mol = Chem.MolFromSmiles(it[0])
                if mol is None:
                    continue
                mol.SetProp('_Name', it[1])
            else:
                mol = obj
            yield mol

    def close(self):
        try:
            self.ifs.close()
        except:
            pass

class MolWriter:
    def __init__(self, fname):
        ext = getExtension(fname)
        self.ostrm = None
        if ext in ['ism', 'smi', 'smiles']:
            self.ofs = Chem.rdmolfiles.SmilesWriter(fname)
        elif ext in ['ism.gz', 'smi.gz', 'smiles.gz']:
            self.ostrm = gzip.open(fname, 'wt')
            self.ofs = Chem.rdmolfiles.SmilesWriter(self.ostrm)
        elif ext == 'sdf':
            self.ofs = Chem.rdmolfiles.SDWriter(fname)
        elif ext == 'sdf.gz':
            self.ostrm = gzip.open(fname, 'wt')
            self.ofs = Chem.rdmolfiles.SDWriter(self.ostrm)
        else:
            raise Exception(f'Unknown output format "{ext}"')

    def write(self, mol):
        try:
            self.ofs.write(mol)
            return True
        except:
            return False

    def close(self):
        self.ofs.close()
        try:
            self.ostrm.close()
        except:
            pass

class Timer:
    def __init__(self, dt=10.0):
        self.dt = dt
        self.T0 = time.time()
        self.T = self.T0

    def check(self):
        ret = False
        t = time.time()
        if self.dt < t - self.T:
            self.T += self.dt
            ret = True
        return ret

    def elapsed(self):
        return time.time() - self.T0
