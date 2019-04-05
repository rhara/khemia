def showattr(v):
    print(type(v), ' '.join([d for d in dir(v) if d[0] != '_']))
