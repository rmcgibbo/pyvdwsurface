"""pyvdwsurface

"""

DOCLINES = __doc__.split("\n")

import sys
from distutils.core import setup
from distutils.extension import Extension

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)
Programming Language :: C
Programming Language :: Python
Operating System :: OS Independent
"""

from Cython.Distutils import build_ext
src = ['pyvdwsurface.pyx', 'src/dotsphere.cc', 'src/vdwsurface.cc']
libraries = []
if sys.platform == 'darwin':
    libraries.append('stdc++')

ext = Extension(
    "pyvdwsurface", src, language='c++', include_dirs=['include/'],
    libraries=libraries)

setup(name='pyvdwsurface',
  author='Robert McGibbon',
  author_email='rmcgibbo@gmail.com',
  url = "http://github.com/rmcgibbo/pyvdwsurface",
  description=DOCLINES[0],
  long_description="\n".join(DOCLINES[2:]),
  version='0.1',
  license='LGPLv2+',
  cmdclass={'build_ext': build_ext},
  ext_modules = [ext]
)
