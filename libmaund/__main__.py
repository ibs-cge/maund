import argparse
import logging
import uuid

from . import run_maund

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
ch.setFormatter(formatter)
logger = logging.getLogger('maund.'+str(uuid.uuid1())[2:8])
logger.setLevel(logging.DEBUG)
logger.addHandler(ch)

parser = argparse.ArgumentParser(description='MAUND')
#parser.add_argument("-v", "--verbosity", action="count", default=0,
#                    help="increase output verbosity")
parser.add_argument('aseq')
parser.add_argument('rgen')
parser.add_argument('files', nargs='*')
parser.add_argument('-c','--comparison_range', type=int, default=60)
parser.add_argument('-b','--window_beg', type=int, default=4)
parser.add_argument('-e','--window_end', type=int, default=7)
parser.add_argument('-ib','--idxseq_beg', type=int, default=13)
parser.add_argument('-ie','--idxseq_end', type=int, default=22)
parser.add_argument('-t','--target_nt',  default='A', choices=['A','C','G','T'])
parser.add_argument('-mcut','--mismatch_cutoff', type=int, default=4)
parser.add_argument('-name','--target_name', type=str, default=None)
parser.add_argument('-otag','--output_nametag', type=str, default='out')
args = parser.parse_args()

if __name__=='__main__':
    run_maund(args, logger)

