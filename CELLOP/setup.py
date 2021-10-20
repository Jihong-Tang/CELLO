import sys
import os
import platform
from distutils.core import setup

"""
setup script for CELLO -- untangling Cancer EvoLution from LOngitudinal sequencing
"""

#if sys.version_info[0] != 2 or sys.version_info[1] < 7:
    #print >> sys.stderr, "ERROR: CELLO requires Python 2.7"


def main():
    setup(
        name='CELLO',
        version='1.0',
        description='CELLO (untangling Cancer EvoLution from LOngitudinal sequencing)',
        packages=['mutmodule'],
        package_dir={'foo': 'lib'},
        scripts=['bin/mutProcess.py', 'bin/mutLandscape.py',
                 'bin/mutCorrelation.py', 'bin/mutCorrelation.py',
                 'bin/mutSignature.py'],
        author="Jihong Tang",
        platforms=['Linux', 'MacOS'],
        long_description="CELLO is is a toolbox for comprehensive analysis of longitudinal genomic sequencing data in cancer.",
    )


if __name__ == "__main__":
    main()
