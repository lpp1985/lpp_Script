import sys
import os
import logging
import numpy as np

log = logging.getLogger(__name__)


class CsvReader(object):

    """
    Load a csv file into a numpy array.
    Typical usage would be:

    ReadId,#reads
    'myRead', 4

    #provide a dtype dict for indexing into the array
    >dtype = {  "ReadId" : ( "|S64", str ),
                "#Bases"  : ( int,    int ) }
    >c =CsvReader(dtype, 'my.csv')


    #print the ReadId value for each row
    >c.load()
    >data = c.data_as_numpy_array()
    >for row in data:
    >    print( row('ReadId')

    #interate over the 1st dimension, print the value
    >for i in b:
    >    print( i['#Bases')
    """

    def __init__(self, path, column_dtypes, strict):
        if not os.path.exists(path):
            raise IOError('File does not exist: {f}'.format(f=path))

        self._column_dtypes = column_dtypes
        self._path = path
        self._np_array = None
        self._num_records = 0
        self._strict = strict

    def load(self):
        '''
        Load data into a numpy array.
        Should be called after reader construction.
        '''
        titles = None
        f = None
        try:
            data = []
            f = open(self._path, "r")
            for line in f:
                line_parts = line.strip().split(",")
                if titles is None:
                    titles = self._validate_titles(line_parts)
                else:
                    line_values = [line_parts[v[1]] for v in titles]
                    data.append(tuple(line_values))
                    self._num_records += 1

            dTypes = [(name, self._column_dtypes[name][0])
                      for name, idx in titles]
            self._np_array = np.array(data, dtype=dTypes)
        finally:
            if f is not None:
                f.close()

    def _validate_titles(self, titles):
        """
        Return a dict of title to position in csv. Each title map to the keys in _column_dtypes
        Raise ValueError if a key in the column_dtypes dict is not in the titles list.
        :param titles: (list) list of the column header strings
        """
        title_tuples = []
        for t in self._column_dtypes.keys():

            if t not in titles:
                if self._strict:
                    raise ValueError('Column title {t} in map does not exist in csv titles {c}'.format(t=t,
                                                                                                       c=str(titles)))
                else:
                    continue
            title_tuples.append((t, titles.index(t)))

        return sorted(title_tuples, key=lambda x: x[1])

    def data_as_numpy_array(self):
        '''
        Get csv row data as a np array
        '''
        return self._np_array

    @property
    def num_records(self):
        '''
        Return the number of data records in the csv file
        '''
        return self._num_records
