#!/usr/bin/env python
"""Python implementation of Philips Mira decoding logic.
"""

import numpy as np
import struct
from io import IOBase

# little endian single element read functions
#   -reads from either file objects or buffers
def getInt(fil):
    fmt = '<i'  # little endian: <, int: i
    bytesize = 4
    if isinstance(fil, IOBase):
        return(struct.unpack(fmt, fil.read(bytesize))[0])
    else:
        return(struct.unpack(fmt, fil[:bytesize])[0])


def getShort(fil):
    fmt = '<h'  # little endian: <, short: h
    bytesize = 2
    if isinstance(fil, IOBase):
        return(struct.unpack(fmt, fil.read(bytesize))[0])
    else:
        return(struct.unpack(fmt, fil[:bytesize])[0])

# from BitStream.java
# a bit shift get function


def getBits(bitbuf, b, signed, val, p):

    t = int(0)
    nb = 32  # bits per integer
    read_cnt = 0

    if p == 0:
        val = np.int32(getInt(bitbuf))
        read_cnt += 4
        t = val
        p = nb - b
    else:
        t = np.int32(val) << (nb - p)
        if p >= b:
            p -= b
        else:
            val = np.int32(getInt(bitbuf))
            read_cnt += 4
            t |= np.uint32(val) >> p
            p += nb - b

    if signed:
        tmp = np.int32(t) >> (nb - b)
        return read_cnt, tmp, val, p
    else:
        return read_cnt, np.uint32(t) >> (nb - b), val, p

# decode coil-AtoD transmission


def decode(bitbuf, out_data_arr, offset, nr_out_elements):

    read_cnt = 0  # for keeping track of buffer reads
    assn_cnt = 0  # number of actual output elems
    val = int(0)
    p = int(0)

    end_pos = offset + nr_out_elements
    i = int(offset)
    while i < end_pos:

        rcnt, b, val, p = getBits(bitbuf, 5, False, val, p)
        bitbuf = bitbuf[rcnt:]
        read_cnt += 1
        rcnt, s, val, p = getBits(bitbuf, 5, False, val, p)
        bitbuf = bitbuf[rcnt:]
        r = np.uint32(1 << s) >> 1
        read_cnt += 1

        # for n in xrange(min(16,end_pos-i),0,-1):
        n = min(16, end_pos-i)
        while n > 0:
            try:
                rcnt, tmp, val, p = getBits(bitbuf, b, True, val, p)
                out_data_arr[i] = (tmp << s) + r
                i += 1
                assn_cnt += 1
            except:
                print("i:"+str(i)+"; n:"+str(n)+"; len(bitbuf): "+str(len(bitbuf)))
                raise
            bitbuf = bitbuf[rcnt:]
            read_cnt += 1
            n -= 1

    return(read_cnt, assn_cnt)


def decodeMira(bytebuff, temp_data, total_data_size):

    temp_data_size = 0

    while temp_data_size < total_data_size:

        data_size = getShort(bytebuff)
        bytebuff = bytebuff[2:]  # consume 2 bytes

        tmp_coded_size = getShort(bytebuff)
        bytebuff = bytebuff[2:]  # consume 2 bytes

        offset1 = getInt(bytebuff)
        bytebuff = bytebuff[4:]  # consume 4 bytes

        temp_data_size += data_size

        read_cnt, assn_cnt = decode(
            bytebuff, temp_data, offset1/4, data_size/4)

        # decode does not actually modify bytebuff.
        # need to remove the appropriate number of samples.
        bytebuff = bytebuff[(
            tmp_coded_size):]

        return temp_data
