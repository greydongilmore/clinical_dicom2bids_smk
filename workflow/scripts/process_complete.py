#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 20:14:39 2020

@author: greydon
"""
import datetime
import os

def get_parser():
	"""
	Argument Parser
	"""
	from argparse import ArgumentParser
	
	parser = ArgumentParser(description=('Group bids sessions into clinically meaningful'))
	
	# Required arguments
	g_req = parser.add_argument_group('required arguments')
	g_req.add_argument('--out_file', action='store', required=True, help='the directory with dicom data')
	
	return parser

def main(args):
	if not os.path.exists(os.path.dirname(args.out_file)):
		os.makedirs(os.path.dirname(args.out_file))
		
	with open(args.out_file, "w") as text_file:
		text_file.write("Completed at: {0}".format(datetime.datetime.now().strftime('%b_%d_%Y %H:%M:%S')))
	
if __name__ == "__main__":
	
	args = get_parser().parse_args()
	
	main(args)