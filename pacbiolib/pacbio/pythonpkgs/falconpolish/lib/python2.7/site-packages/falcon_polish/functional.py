"""Purely functional, somewhat generic code.

(Of course, these might exhaust input iterators.)
"""
import re

def total_length(pairs):
    return sum(length*count for (length, count) in pairs)

def calc_cutoff(target, pairs):
    """Return read_length such that sum(lens for len >= rl) >= target.
    Raise on empty pairs, which are (length, count) tuples.
    """
    accum = 0
    for length, count in reversed(sorted(pairs)):
        accum += length*count
        if accum >= target:
            break
    else:
        raise Exception('Total=%d < target=%d' %(accum, target))
    return length

def fns_from_fofn(fofn):
    for line in fofn:
        fn = line.strip()
        if not fn:
            continue
        yield fn

def joined_strs(pieces, olen):
    """Reduce len to olen by joining some strings.

    >>> list(joined_strs(['a', 'b', 'c'], 2))
    ['ab', 'c']
    """
    ilen = len(pieces)
    rem = ilen
    while rem:
        n = ((rem-1)//olen) + 1
        yield ''.join(pieces[ilen-rem : n+ilen-rem])
        rem -= n
        olen -= 1

def stricter_json(json_str):
    """Remove trailing commas.
    This works except in highly unlikely cases.
    http://stackoverflow.com/a/23705538
    >>> stricter_json('{"foo": "bar",}')
    '{"foo": "bar"}'
    """
    json_str = re.sub(r',\s*}', '}', json_str)
    json_str = re.sub(r',\s*]', ']', json_str)
    return json_str
