#! python3

import gzip, codecs
from contextlib import contextmanager

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        f.close()
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            f.close()
            return "utf-16"
        except UnicodeDecodeError:
            print(f"'{fileName}' is neither utf-8 nor utf-16 encoded; please convert to one of these formats.")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

class Emptyfile(object):
    def write(self, data):
        pass # ignore the data
    def __enter__(self): return self
    def __exit__(*x): pass

@contextmanager
def write_conditionally(fileName):
    if fileName == None:
        empty = Emptyfile()
        yield empty
    else:
        if fileName.endswith(".gz"):
            with gzip.open(fileName, "wt") as f:
                yield f
        else:
            with open(fileName, "w") as f:
                yield f

def parse_annotation_table(fileName, delimiter="\t"):
    '''
    Returns:
        dataDict -- a dictionary where keys are column headers and
                    values are the text contents of this row's value.
                    An additional "leftcolumn" key is added with the
                    first column's value which should correspond to
                    the ID.
    '''
    alreadyWarned = False
    header = None
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            sl = line.strip().split(delimiter)
            if header == None:
                header = sl
                headerIndex = { h:i for i,h in enumerate(header) }
                continue
            
            dataDict = {"leftcolumn": sl[0]}
            try:
                dataDict.update({ h:sl[i] for h, i in headerIndex.items() })
            except:
                dataDict.update({ h:"." for h, i in headerIndex.items() })
                if not alreadyWarned:
                    print(f"WARNING: line '{sl}' does not have '{len(headerIndex)}' columns as expected;" + 
                          " this and similar gene lines will have '.' values imputed for all columns")
                    alreadyWarned = True
            yield dataDict
