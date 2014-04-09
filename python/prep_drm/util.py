import bz2
import collections
import functools
import itertools
import gzip
import os.path
import sys

# From http://wiki.python.org/moin/PythonDecoratorLibrary
class memoized(object):
   """Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   """
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)


def isolate_region(aligned_pairs, start, end):
    """
    Prune a collection of aligned_pairs to a region of interest.
    returns the region from [start,end) on the reference sequence.
    """
    e = enumerate(aligned_pairs)

    istart = next((i for i, (q, r) in e if r is not None and r >= start and r <= end), None)
    iend = next((i for i, (q, r) in e if r is not None and r >= end), None)

    if istart is None and iend is None:
        return []

    return aligned_pairs[istart:iend]


def opener(mode, *args, **kwargs):
    """
    Open a file, with optional compression based on extension
    """
    exts = {'.bz2': bz2.BZ2File,
            '.gz': gzip.open}
    def open_file(path):
        if path == '-':
            if mode.startswith('r'):
                return sys.stdin
            else:
                return sys.stdout

        open_fn = exts.get(os.path.splitext(path)[1], open)
        return open_fn(path, mode, *args, **kwargs)
    return open_file

# From http://docs.python.org/release/2.3.5/lib/itertools-example.html
def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
