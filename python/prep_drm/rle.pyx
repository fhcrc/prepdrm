# distutils: language = c
import re
from libc.stdlib cimport free, malloc
from libc.string cimport strdup

cdef class RLEItem:
    cdef char c
    cdef size_t length

    def __cinit__(self, char c, size_t length):
        self.c = c
        self.length = length

    property c:
        def __get__(self):
            return chr(self.c)
        def __set__(self, char c):
            self.c = c

    property length:
        def __get__(self):
            return self.length
        def __set__(self, size_t l):
            self.length = l

    def __len__(self):
        return self.length

    def __repr__(self):
        return '<RLEItem [{0},{1}]>'.format(chr(self.c), self.length)

def encode(bytes b not None):
    cdef size_t blen = len(b), i = 0
    if blen == 0:
        return []

    cdef list result = []
    cdef char cur = b[0], c
    cdef size_t l = 1

    for i in range(1, blen):
        c = b[i]
        if c != cur:
            result.append(RLEItem.__new__(RLEItem, cur, l))
            cur = c
            l = 1
        else:
            l += 1
    result.append(RLEItem.__new__(RLEItem, cur, l))

    return result

def decode(list items not None):
    cdef size_t i = 0, j, l = 0
    cdef bytes res_str
    cdef char* res
    cdef RLEItem item
    for item in items:
        l += item.length
    if l == 0:
        return ''
    res = <char*>malloc(sizeof(char) * (l + 1))
    res[l] = '\0'
    try:
        for item in items:
            for j in range(item.length):
                res[i] = item.c
                i = i + 1
        assert i == l
        res_str = res
        return res_str
    finally:
        free(res)
