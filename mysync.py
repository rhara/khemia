from ftplib import FTP
from urllib import request
import os


class MyFTP:
    def __init__(self, host, user='anonymous', password='jvalao@auvia.com'):
        self.host = host
        self.user = user
        self.password = password

        self.ftp_instance = FTP(host, user, password)
        self.remote_dir = None
        self.ls = None
        self.local_dir = '.'

    def set_remote_list(self, remote_dir):
        lines = []
        self.remote_dir = remote_dir
        self.ftp_instance.dir(self.remote_dir, lines.append)
        self.ls = {}
        for line in lines:
            it = line.rstrip().split()
            name = it[-1]
            size = int(it[4])
            self.ls[name] = size

    def set_local_dir(self, local_dir):
        self.local_dir = local_dir

    def get_size(self, name):
        return self.ls[name]

    def download(self, name, mode='binary', check_size=True):
        dest_name = os.path.join(self.local_dir, name)
        local_size = None
        try:
            st = os.stat(dest_name)
            local_size = st.st_size
        except:
            pass
        try:
            remote_size = self.get_size(name)
            if check_size and (local_size is None or local_size != remote_size):
                self.ftp_instance.cwd(self.remote_dir)
                print(f'{dest_name} <- {self.host}:{self.remote_dir}/{name}')
                if mode[0] == 'b':
                    with open(dest_name, 'wb') as f:
                        self.ftp_instance.retrbinary(f'RETR {name}', f.write)
                else:
                    with open(dest_name, 'wt') as f:
                        self.ftp_instance.retrlines(f'RETR {name}', f.write)
            else:
                print(f'{dest_name} up-to-date')
        except KeyboardInterrupt as e:
            raise e
        except:
            print(f'{dest_name} download error')


class MyHTTP:
    def __init__(self, base):
        self.base = base
        self.local_dir = '.'

    def set_local_dir(self, local_dir):
        self.local_dir = local_dir

    def download(self, name, check_size=True):
        dest_name = os.path.join(self.local_dir, name)
        src = f'{self.base}/{name}'
        local_size = None
        try:
            st = os.stat(dest_name)
            local_size = st.st_size
        except:
            pass
        try:
            res = request.urlopen(src)
            remote_size = int(res.getheader('Content-Length'))
            if check_size and (local_size is None or local_size != remote_size):
                print(f'{dest_name} <- {src}')
                request.urlretrieve(src, dest_name)
            else:
                print(f'{dest_name} up-to-date')
        except KeyboardInterrupt as e:
            raise e
        except:
            print(f'{dest_name} download error')

