from setuptools import setup
from khemia import version

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='khemia',
    version=version.KHEMIA_VERSION,
    author='Ryuichiro Hara',
    author_email='r@harara.com',
    description='Khemia Chemistry higher library',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['khemia'],
    scripts=['tools/khconvert.py', 'tools/khcount.py', 'tools/khcut.py', 'tools/khgrep.py'],
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
