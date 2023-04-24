#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #############################################
import glob
import os
import shutil
###############
# File Select
def file_listing():
	print('ttttt')
	target_list = ['**/all.dat', '**/ave.dat', '**/averaged.png', '**/plot_all.plt', '**/series.png']
	for target in target_list:
		src_list = glob.glob(target, recursive=True)
		print(src_list)
		for src in src_list:
			dst = os.path.join('./data', src)
			dirname = os.path.dirname(dst)
			os.makedirs(dirname, exist_ok=True)
			shutil.copy(src, dst)
	return 

##########################
# Main
if __name__=='__main__':
	file_listing()

