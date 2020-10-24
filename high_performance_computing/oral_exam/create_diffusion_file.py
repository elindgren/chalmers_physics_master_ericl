#! /usr/bin/python3

import argparse
from random import randint, uniform


parser = argparse.ArgumentParser()
parser.add_argument(
  '-r', type=int, default=3, help='number of rows')
parser.add_argument(
  '-c', type=int, default=3, help='number of columns')
parser.add_argument(
  '-n', type=int, default=3, help='number of fields')
parser.add_argument(
  '-v', type=float, default=3, help='maximal value')
parser.add_argument(
  'output_name', metavar='output', type=str,
  help='output file name')
args = parser.parse_args()


def create_point_string(r,c,v):
  return "{} {} {}\n".format(randint(0,r-1), randint(0,c-1), uniform(0.,v))


with open(args.output_name, 'w') as output:
  output.write("{} {}\n".format(args.r,args.c))
  output.writelines( create_point_string(args.r, args.c, args.v)
                     for _ in range(args.n) )
