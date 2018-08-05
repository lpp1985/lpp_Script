"""Define Read Annotation class."""

__all__ = ["ReadAnnotation"]

def xorNA(x):
    """Return x if x is not None, or return 'NA'."""
    return str(x) if x is not None else 'NA'


def intorNone(x):
    """Return int(x) if x is not None."""
    return int(x) if x is not None else x


def hasNonPropertyAttr(obj, attr):
    """Return true if obj has a non-property attribute 'attr'."""
    if hasattr(obj, attr) and \
       (not hasattr(type(obj), attr) or
            not isinstance(getattr(type(obj), attr), property)):
        return True
    else:
        return False


class ReadAnnotation(object):

    """Read annotation class, including the following fields
       "ID", "strand",
       "fiveseen", "polyAseen", "threeseen",
       "fiveend", "polyAend", "threeend",
       "primer"
    """

    def __init__(self, ID, strand=None, fiveend=None, polyAend=None,
                 threeend=None, primer=None, chimera=None, ignore_polyA=False):
        self.ID = ID
        self.strand = strand
        self.fiveend = intorNone(fiveend)
        self.polyAend = intorNone(polyAend)
        self.threeend = intorNone(threeend)
        self.primer = intorNone(primer)
        self.chimera = intorNone(chimera)
        self.ignore_polyA = ignore_polyA

    @classmethod
    def fromString(cls, line, ignore_polyA=False):
        """Construct and return a ReadAnnotation object."""
        try:
            sid, desc = line.strip().split(' ')[0:2]
            ret = ReadAnnotation(sid, ignore_polyA=ignore_polyA)
            for d in desc.split(';'):
                attr, val = d.split('=')
                if hasNonPropertyAttr(ret, attr):
                    if attr == "strand" or attr == "ID":
                        setattr(ret, attr, val if val not in ["None", 'NA']
                                else None)
                    else:
                        setattr(ret, attr, int(val) if val not in ["None", 'NA']
                                else None)
            return ret
        except Exception as e:
            errMsg = "String not recognized as a valid read annotation:" + \
                str(e)
            raise ValueError(errMsg)

    @property
    def fiveseen(self):
        """Return whether 5' primer has been seen."""
        return 1 if self.fiveend is not None else 0

    @property
    def polyAseen(self):
        """Return whether polyA tail has been seen."""
        return 1 if (self.polyAend is not None and
                     self.polyAend >= 0) else 0

    @property
    def threeseen(self):
        """Return whether 3' primer has been seen."""
        return 1 if self.threeend is not None else 0

    @property
    def isFullLength(self):
        """Return whether it is a full length read or not."""
        if self.threeseen == 1 and \
           self.fiveseen == 1 and \
           (self.polyAseen == 1 or self.ignore_polyA):
            return True
        else:
            return False

    @staticmethod
    def fieldsNames():
        """Return names of annotation fields."""
        return ("ID", "strand", "fiveseen", "polyAseen", "threeseen",
                "fiveend", "polyAend", "threeend", "primer", "chimera")

    @staticmethod
    def header(delimiter=","):
        """Return header string of annotation fields."""
        fds = list(ReadAnnotation.fieldsNames())
        assert(fds[0] == "ID")
        fds[0] = "id"
        # bug 25728, because of a known word xls bug which mistakenly
        # recognizing a csv file starting with 'ID' as a SYLK file, we have to
        # print field name as "id" rather than "ID".
        return delimiter.join(fds)

    def fields(self):
        """Return all fields defined within class."""
        return (getattr(self, attr, None)
                for attr in ReadAnnotation.fieldsNames())
        #(self.ID, self.strand, self.fiveseen, self.polyAseen,
        #        self.threeseen, self.fiveend, self.polyAend,
        #        self.threeend, self.primer)

    def __repr__(self):
        return "\t".join([xorNA(x) for x in self.fields()])

    def toReportRecord(self, delimitor=","):
        """Return a string to represent an annotation in a report."""
        return delimitor.join([xorNA(x) for x in self.fields()])

    def toAnnotation(self):
        """Return a string to represent an annotation."""
        return str(self.ID) + " " + ";".join(
            [x + "=" + xorNA(y) for x, y in
             zip(self.fieldsNames(), self.fields())[1:]])
