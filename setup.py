#Compilation file
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("c_antisense",sources=["c_antisense.pyx","antisense.c","utils.c"],extra_compile_args=["-std=c99","-Qunused-arguments","-Wno-error=unused-command-line-argument-hard-error-in-future"])]

setup(
  name = 'Antisense C Module',
    cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules
      )
