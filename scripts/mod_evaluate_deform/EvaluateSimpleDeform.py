#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #############################################
from UDFManager import *
import argparse
from pathlib import Path
import glob
# from operator import itemgetter
# import os
import platform
import subprocess
import sys

import mod_evaluate_deform.variables as var
###########################################################
# print("This is module!")
###########################################################
def simple_deform():
	read_arg()
	if not var.f_average:
		read_stress()
	elif var.f_average:
		average()
	#
	plot()

	return


##############
# Read argument 
def read_arg():
	parser = argparse.ArgumentParser(description='Evaluate deformed simulations !')
	parser.add_argument('-f','--func', type=int, help="Functionality of junction point (int).")
	parser.add_argument('-n', '--nu', type=float, help="Strand density of network (float).")
	parser.add_argument('-m', '--mode', help="Mode of deformation; shear or stretch")
	parser.add_argument('-a', '--average', help="Average multi data of different deformation", action='store_true')
	parser.add_argument('-s', '--series', help="Average multi data of different deformation", action='store_true')
	args = parser.parse_args()
	if args.func and args.nu:
		var.func = args.func
		var.nu = args.nu
	else:
		print('\n#####\nfunctionality and/or nu is not specified')
		print('Default value will be used!')
	if args.mode:
		print('deform mode is set to ', args.mode.lower())
		var.simple_def_mode = args.mode.lower()
	else:
		print('\n#####\ndeformation mode is not set!')
		print('according to file name(shear or stretch), evaluation mode  will be set!')
	if args.average:
		var.f_average = True
	elif args.series:
		var.f_series = True
	return

# Read Stress data
def read_stress():
	outudf = '*_out.udf'
	udf_list = glob.glob(outudf)
	if udf_list:
		var.target_udf = udf_list[0]
		if var.simple_def_mode == '':
			sys.exit('\n#####\nfile name is not start from either shear or stretch.\ndefault mode of stretch will be used!')
	else:
		sys.exit('\n#####\nNo effective *_out.udf file in this directory !!\nSomething wrong !!\n')
	#
	var.target_rate = Path.cwd().parent.name
	var.target_rotate = Path.cwd().name
	#
	var.ss_data = read_and_calc()
	save_data()
	return

# Read Data
def read_and_calc():
	print("Readin file = ", var.target_udf)
	uobj = UDFManager(var.target_udf)
	data = []
	#
	uobj.jump(0)
	cell = uobj.get("Structure.Unit_Cell.Cell_Size")

	for i in range(1, uobj.totalRecord()):
		print("Reading Rec.=", i)
		uobj.jump(i)
		if var.simple_def_mode == 'shear':
			stress = uobj.get('Statistics_Data.Stress.Total.Batch_Average.xy')
			strain = uobj.get('Structure.Unit_Cell.Shear_Strain')
		elif var.simple_def_mode == 'stretch':
			cell = uobj.get("Structure.Unit_Cell.Cell_Size")
			stress_list = uobj.get("Statistics_Data.Stress.Total.Batch_Average")
			stress = stress_list[2]-(stress_list[0] + stress_list[1])/2.
			strain = uobj.get("Structure.Unit_Cell.Cell_Size.c")/(cell[0]*cell[1]*cell[2])**(1/3)
		data.append([strain, stress])
	return data

# 計算結果をターゲットファイル名で保存
def save_data():
	var.ss_data_file = "SS_" + var.target_rotate + '_' + var.target_rate + '.dat'
	with open(var.ss_data_file,'w') as f:
		for line in var.ss_data:
			f.write(str(line[0]) + '\t' + str(line[1]) + '\n')
	return

#####
# average
def average():
	var.ss_data_file = "SS_" + var.target_rate + '_all.dat'
	dat_list = glob.glob('./*/SS_rotate*.dat')
	gathered = []
	for file in dat_list:
		with open(file, 'r') as f:
			for line in f.readlines():
				if line[0] not in ['#', '', '\n']:
					gathered.append(line)
	with open(var.ss_data_file,'w') as f:
		for line in gathered:
			f.write(line)
	return

############################
# 結果をプロット
def plot():
	plt_file = 'plot.plt'
	SS_script()
	if var.simple_def_mode == 'stretch':
		MR_script()
	#
	with open(plt_file, 'w') as f:
		f.write(var.script)
	#
	if platform.system() == "Windows":
		subprocess.call([plt_file], shell=True)
	elif platform.system() == "Linux":
		subprocess.call(['gnuplot ' + plt_file], shell=True)
	return

# スクリプトの中身
def SS_script():
	var.script = 'set term pngcairo font "Arial,14"\n\n'
	var.script += '#set mono\nset colorsequence classic\n\n'
	var.script += 'data = "' + var.ss_data_file + '"\n'
	if var.f_average:
		var.script += 'set output "SS_' + var.target_rate + '_averaged.png"\n\n'
	else:
		var.script += 'set output "SS_' + var.target_rate + '.png"\n\n'
	var.script += 'set key left\nset size square\n'
	var.script += 'set xlabel "{/Symbol l}^2-1/{/Symbol l}"\nset ylabel "{/Symbol s}_{true}"\n\n'
	var.script += 'G=' + str(var.nu) + '\nfunc=' + str(var.func) + '\nphi=(1.-2./func)\n'
	if var.simple_def_mode == 'stretch':
		var.script += 'set xrange [0.:]\nset yrange [0.:]\n#set xtics 0.5\n#set ytics 0.01\n\n'
		var.script += 'gamma=1.\np(x)=gamma*G*phi*(x**2-1./x)\nfit [1.1:2] p(x) data via gamma\n\n' 
		var.script += 'plot	'
		if var.f_average:
			var.script += 'data usi ($1**2-1./$1):($2) lc rgb "gray" ti "raw data", \\\n'
			var.script += 'data usi ($1**2-1./$1):($2) smooth unique w l lw 2 lt 1 ti "averaged", \\\n'
		else:
			var.script += 'data usi ($1**2-1./$1):($2) w l lw 2 lt 1 ti "rate: ' + var.target_rate.split('_')[1] + '", \\\n'
		var.script += 'G*x w l lw 2 lt 2 ti "Affine", \\\nG*phi*x w l lw 2 lt 3 ti "Phantom", \\\ngamma*G*phi*x w l  lw 2 lt 4 ti "Fitted"'
		var.script += '\n\nreset\n\n'
	elif var.simple_def_mode == 'shear':
		var.script += 'set xrange [0:]\nset yrange [0.:]\n#set xtics 0.5\n#set ytics 0.01\n'
		var.script += 'gamma=1.\np(x)=gamma*G*phi*x\nfit [0.1:1] p(x) data via gamma\n\n' 
		var.script += 'plot	'
		var.script += 'data w l lw 2 lt 1 ti "rate: ' + var.target_rate.split('_')[1] + '", \\\n'
		var.script += 'G*x w l lw 2 lt 2 ti "Affine", \\\nG*phi*x w l lw 2 lt 3 ti "Phantom", \\\ngamma*G*phi*x w l  lw 2 lt 4 ti "Fitted"'
		var.script += '\n\nreset\n\n'
	
	return

# スクリプトの中身
def MR_script():
	var.script += 'set term pngcairo font "Arial,14"\n\n'
	var.script += '#set mono\nset colorsequence classic\n\n'
	var.script += 'data = "' + var.ss_data_file + '"\n'
	if var.f_average:
		var.script += 'set output "MR_' + var.target_rate + '_averaged.png"\n\n'
	else:
		var.script += 'set output "MR_' + var.target_rate + '.png"\n\n'
	var.script += 'set key left\nset size square\n'
	var.script += 'set xrange [0:1]\nset yrange [0.:]\n#set xtics 0.5\n#set ytics 0.02\n'
	var.script += 'set xlabel "1/{/Symbol l}"\nset ylabel "{/Symbol s}_{true}/({/Symbol l}^2-1/{/Symbol l})"\n\n'
	var.script += '## Fit Range\n\nlow = 0.5\nhigh = 0.7\n\n'
	var.script += 'fit [low:high] a*x+b data usi ( 1/$1 ):( $2/( $1**2 - 1/$1 ) ) via a,b\n\n'
	var.script += 'set label 1 sprintf("C1 = %.3f", b/2) left at graph 0.2,0.7\n'
	var.script += 'set label 2 sprintf("C2 = %.3f", a/2) left at graph 0.2,0.6\n'
	var.script += 'set label 3 sprintf("fit range = %.3f to %.3f", low, high) left at graph 0.2,0.5\n\n#\n'
	var.script += 'plot	'
	if var.f_average:
		var.script += 'data usi ( 1/$1 ):( $2/( $1**2 - 1/$1 ) ) lc rgb "gray" ti "raw data", \\\n'
		var.script += 'data usi ( 1/$1 ):( $2/( $1**2 - 1/$1 ) ) smooth unique w lp pt 7 lt 1 ti "averaged", \\\n'
	else:
		var.script += 'data usi ( 1/$1 ):( $2/( $1**2 - 1/$1 ) ) w lp pt 7 lt 1 ti "Original Data", \\\n'
	var.script += '[low:high] a*x + b w l lw 5 lt 2 ti "Fitted Line"'
	var.script += '\n\nreset'

	return