from optparse import OptionParser
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-d", "--Database", action="store",
                  dest="DB_FILE",

                  help="Database File")


parser.add_option("-n", "--n", action="store",
                  dest="db_num",
                  type='int',
                  help="database number")


parser.add_option("-o", "--OUTPUT", action="store",
                  dest="output",

                  help="TSV output")