cdef extern from "python/pcg_pysolver.h":
  int pcg_pysolver(const char* name)

def run(name):
  pcg_pysolver(name)
