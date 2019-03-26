from ftplib import FTP
from urllib import request
import os, re


class MyFTP:
    def __init__(self, base, user='anonymous', password='jvalao@auvia.com'):
        pat = re.compile('^ftp://([^/]+)(/.*)$')
        self.base = base
        m = pat.match(base)
        host = m.group(1)
        self.remote_dir = m.group(2)

        self.ftp_instance = FTP(host, user, password)
        self.ftp_instance.cwd(self.remote_dir)
        self.local_dir = '.'

    def get_file_list(self):
        lines = []
        ls = []
        self.ftp_instance.retrlines('LIST', lines.append)
        for line in lines:
            it = line.strip().split()
            name = it[8]
            ls.append(name)
        return ls

    def set_local_dir(self, local_dir):
        self.local_dir = local_dir

    def download(self, name, mode='binary', check_size=True):
        dest_name = os.path.join(self.local_dir, name)
        local_size = None
        try:
            st = os.stat(dest_name)
            local_size = st.st_size
        except:
            pass
        try:
            lines = []
            self.ftp_instance.dir(name, lines.append)
            remote_size = int(lines[0].strip().split()[4])
            megabytes = remote_size/1024/1024
            if check_size and (local_size is None or local_size != remote_size):
                dir_name = os.path.dirname(dest_name)
                os.makedirs(dir_name, exist_ok=True)
                self.ftp_instance.cwd(self.remote_dir)
                print(f'{dest_name} {megabytes:.2f}MB <- {self.base}/{name}')
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
        except Exception as e:
            print(e)
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
            megabytes = remote_size/1024/1024
            if check_size and (local_size is None or local_size != remote_size):
                dir_name = os.path.dirname(dest_name)
                os.makedirs(dir_name, exist_ok=True)
                print(f'{dest_name} {megabytes:.2f}MB <- {src}')
                request.urlretrieve(src, dest_name)
            else:
                print(f'{dest_name} up-to-date')
        except KeyboardInterrupt as e:
            raise e
        except:
            print(f'{dest_name} download error')

