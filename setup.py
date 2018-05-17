"""Installation instructions for cisVar."""
import os
from setuptools import setup

import cisVar  # For version

VERSION=cisVar.__version__
GITHUB='https://github.com/TheFraserLab/cisVar'

REQUIREMENTS = ['pandas', 'numpy', 'scipy', 'psutil', 'snakemake', 'wget']


def read(fname):
    """Read the contents of a file in this dir."""
    with open(os.path.join(os.path.dirname(__file__), fname)) as fin:
        return fin.read()


# Actual setup instructions
setup(
    name         = 'cisVar',
    version      = VERSION,
    author       = 'Ashley Tehranchi',
    author_email = 'mike.dacre@gmail.com',
    description  = (
        "Calculate both pre and post frequencies for ChIP or ATAC style data"
    ),
    keywords = (
        "ATACseq ChIPseq regression bioinformatics"
    ),
    long_description = read('README.md'),
    license = 'MIT',

    # URLs
    url = GITHUB,
    download_url='{0}/archive/v{1}.tar.gz'.format(GITHUB, VERSION),

    py_modules=['cisVar'],

    # Entry points and scripts
    # Entry points are functions that can use sys.argv (e.g. main())
    # Scripts are independent pieces of code intended to be executed as is
    entry_points = {
        'console_scripts': [
            'cisVar = cisVar:main',
        ],
    },
    scripts = ['cisVar.py', 'regression_qtls.R',
               'scripts/combine.py', 'scripts/plot_fun.R'],


    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Utilities',
    ],

    # Requirements
    requires=REQUIREMENTS,
    install_requires=REQUIREMENTS,
)
