from falcon_polish.pypeflow import hgap
import argparse
import sys

def main(argv=sys.argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--logging',
            help='.ini or .json config file for Python logging module')
    parser.add_argument('config',
            help='.ini or .json of HGAP config. Available sections: "general", "hgap", "falcon", "pbsmrtpipe", "blasr", "quiver", ...')
    args = parser.parse_args(argv[1:])
    return hgap.run(args.config, args.logging)

if __name__ == "__main__":
    main(sys.argv)
