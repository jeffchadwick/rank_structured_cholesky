from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension("rchol", ["rchol.pyx", "pcg_pysolver.cpp"],
              include_dirs = ["/home/dsb253/local/rsc/include/rsc"],
              libraries = ["rsclib"],
              library_dirs = ["/home/dsb253/local/rsc/lib"],
              language="c++",
              extra_compile_args=["-std=c++11"],
              extra_linker_args=["-std=c++11"])
]

setup(
    name = "pcg_solver python binding",
    ext_modules = cythonize(extensions)
)
