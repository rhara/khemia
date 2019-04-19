import os


def showAttr(v):
    print(type(v), ' '.join([d for d in dir(v) if d[0] != '_']))


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


def rename(mol, tagname):
    try:
        title = mol.GetProp(tagname)
        mol.SetProp('_Name', title)
    except:
        pass
