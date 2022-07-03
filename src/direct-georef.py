
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', default = False, action='store_true')
parser.add_argument('-c', '--config', type = str, default='test.json')
args = parser.parse_args()
if args.verbose == True:
    print(args)


