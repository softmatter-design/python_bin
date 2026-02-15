#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #############################################
# from unittest import result
from UDFManager import *
import argparse
import numpy as np
import os
import glob
import math
import cmath
import platform
import subprocess
import sys

from scipy.signal import savgol_filter

import mod_evaluate_deform.variables as var
###########################################################
def step_deform():
	read_arg()	
	if not var.ave_flag:
		calc_gt()
	else:
		average_gt()
		gt2gw()
	return

# Read argument 
def read_arg():
	parser = argparse.ArgumentParser(description='Evaluate deformed simulations !')
	parser.add_argument('-f','--func', type=int, help="Functionality of junction point (int).")
	parser.add_argument('-n', '--nu', type=float, help="Strand density of network (float).")
	parser.add_argument('-m', '--mode', help="Mode of deformation; shear or stretch")
	parser.add_argument('-a', '--average', help="Flag for averaging subdir data", action='store_true')
	args = parser.parse_args()
	if args.func and args.nu:
		var.func = args.func
		var.nu = args.nu
	else:
		print('\n#####\nfunctionality and/or nu is not specified')
		print('Default value will be used!')
	if args.mode:
		var.step_def_mode = args.mode.lower()
	else:
		print('\n#####\ndeformation mode is not set!')
		sys.exit('either mode of shear or stretch should be set!')
	var.ave_flag = args.average
	return

############################
# Calculate stress either for shear or stretch deformation
def calc_gt():
	calc_deform()
	calc_relax()
	return

def calc_deform():
	t_udf = 'deform_out.udf'
	read_data(t_udf)
	return

def calc_relax():
	files = glob.glob("relaxation*out.udf")
	files.sort()
	for t_udf in files:
		read_data(t_udf)
	return

def read_data(t_udf):
	print("Readin file = ", t_udf)
	uobj = UDFManager(t_udf)
	time = []
	strain = []
	g = []
	stress = []
	temp = []
	z_init = read_z(uobj)

	# データ読み込み
	prev_g = 0
	for rec in range(1, uobj.totalRecord()):
		print("Reading Rec.=", rec)
		uobj.jump(rec)
		temp.append(uobj.get('Statistics_Data.Temperature.Batch_Average'))
		if t_udf == 'deform_out.udf':
			time.append(round(uobj.get("Time"), 5))
		else:
			time.append(round(uobj.get("Time"), 5) + var.deform_time)
		if var.step_def_mode == 'stretch':
			stress_list = uobj.get("Statistics_Data.Stress.Total.Batch_Average")
			tmp_stress = stress_list[2]-(stress_list[0] + stress_list[1])/2.
			if t_udf == 'deform_out.udf':
				var.step_strain = round(uobj.get("Structure.Unit_Cell.Cell_Size.c")/z_init, 5)
			tmp_g = tmp_stress/(var.step_strain**2 - 1/var.step_strain)
		elif var.step_def_mode == 'shear':
			tmp_stress = uobj.get('Statistics_Data.Stress.Total.Batch_Average.xy')
			if t_udf == 'deform_out.udf':
				var.step_strain = round(uobj.get('Structure.Unit_Cell.Shear_Strain'), 5)
			tmp_g = tmp_stress/var.step_strain
		#
		if tmp_g < 0:
			tmp_g = prev_g*.8
		strain.append(var.step_strain)
		stress.append(tmp_stress)
		g.append(tmp_g)
		#
		prev_g = tmp_g
	#
	gt = np.stack([time, g, temp], 1)
	ss = np.stack([strain, stress, temp], 1)
	#
	if t_udf == 'deform_out.udf':
		var.step_strain = strain[-1]
		var.deform_time = time[-1]
		var.deform_gt = gt
		var.accumd_gt = gt
		save_data(ss, 'ss_deform.dat')
		save_data(gt, 'gt_deform.dat')
	else:
		var.accumd_gt = np.concatenate([var.accumd_gt, gt])
		file_base = t_udf.rsplit("_", 1)[0]
		save_data(np.concatenate([var.deform_gt, gt]), 'gt_' + str(file_base) + '.dat')
		save_data(var.accumd_gt, "gt_all.dat")
	return

def read_z(uobj):
	uobj.jump(0)
	z_init = uobj.get("Structure.Unit_Cell.Cell_Size")[2]
	return z_init

###############################
#----- 計算結果をターゲットファイル名で保存
def save_data(target, f_data):
	with open(f_data,'w') as f:
		for line in target:
			for data in line:
				f.write(str(data) + '\t')
			f.write('\n')
	plot(f_data)
	return

#----- 結果をプロット
def plot(f_data):
	plt = make_script(f_data)
	#
	if platform.system() == "Windows":
		subprocess.call([plt], shell=True)
	elif platform.system() == "Linux":
		subprocess.call(['gnuplot ' + plt], shell=True)
	return
	
# 必要なスクリプトを作成
def make_script(f_data):
	script = script_content(f_data)
	plt = f_data.replace('dat', 'plt')
	with open(plt, 'w') as f:
		f.write(script)
	return plt

# スクリプトの中身
def script_content(f_data):
	out_png = f_data.replace('dat', 'png')
	script = 'set term pngcairo font "Arial, 14"\n'
	script += 'set colorsequence classic\n\n'
	script += 'data = "' + f_data + '"\n'
	script += 'set output "' + out_png + '"\n\n'
	script += 'set key left\nset size square\n'
	if f_data == 'ss_deform.dat':
		script += '#set xrange [1:4]\n#set yrange [0:0.2]\n#set xtics 1\n#set ytics 0.1\n'
		script += 'set y2tics\nset xlabel "Strain"\nset ylabel "Stress"\nset y2label "Temp."\n\n'
		script += 'plot data u 1:2 axis x1y1 w l lw 2 lt 1 ti "stress", \\\n'
		script += 'data u 1:3 axis x1y2 w l lw 2 lt 2 ti "Temp."'
	elif ('gt' in f_data) and ("position" not in f_data):
		script += '#set xrange [1e-2:1e6]\n#set yrange [1e-2:1e1]\n#set xtics 1\n#set ytics 0.1\n'
		script += 'set xlabel "Time"\nset ylabel "G(t)"\n'	
		script += 'set logscale xy\n\n'	
		script += 'set format x "10^{%L}"\nset format y "10^{%L}"\n'
		script += 'plot	data u 1:2 axis x1y1 w l lw 2 lt 1 ti "G(t)"'
	else:
		script += '#set xrange [1e-2:1e6]\n#set yrange [1e-2:1e1]\n#set xtics 1\n#set ytics 0.1\n'
		script += f'G={var.nu:}\nfunc ={var.func:}\n'
		script += 'f1 = (func - 1.)/(func + 1.)\nphi = 1. - 2./func\n\n'
		script += 'tau = 1000\nst = G\neq = G\ns = 1000\ne = 100000\n\n'
		script += 'g(x) = eq -(eq - st)*exp(-x/tau)\n'
		script += 'fit [s:e] g(x) data via st, eq, tau\n\n'
		script += 'set label 1 sprintf("{/Symbol n}k_BT = %.2e", G) at graph 0.45, 0.7\n\n'
		script += 'set label 2 sprintf("fitting region: from %.d", s) at graph 0.45, 0.6\n'
		script += 'set label 3 sprintf("G_{eq} = %.1e", eq) at graph 0.45, 0.5\n'
		script += 'set label 4 sprintf("{/Symbol t} = %.1e", tau) at graph 0.45, 0.4\n'
		script += 'set xlabel "Time"\nset ylabel "G(t)"\n'	
		script += 'set logscale xy\n\n'	
		script += 'set format x "10^{%L}"\nset format y "10^{%L}"\n'
		script += 'plot	data u 1:2 smooth unique w l lw 2 lt 1 ti "G(t)", \\\n'
		script += '[s:e] g(x) w l lw 2 lt 2 ti "fit", \\\n'
		script += '[1000:] G w l lw 3 dt (10, 5) lt 8 ti "Affin", \\\n'
		script += '[1000:] G*phi w l lw 3 dt (10, 5) lt 7 ti "Phantom"\n\n'
		script += 'set table "gt_averaged.dat"\n'
		script += 'set format "%.5e"\n'
		script += 'plot	data u 1:2 smooth unique \n'
		script += 'unset table\n'
	script += '\n\nreset\n\n'
	return script

###############################
def average_gt():
	dat_list = glob.glob('./*/gt_all.dat', recursive = True)
	res = []
	for dat in dat_list:
		with open(dat, 'r') as f:
			for line in f.readlines():
				res.append([float(line.split('\t')[0]), float(line.split('\t')[1])])
	save_data(res,  'gt_all_position.dat')

	return

##################################
def gt2gw():
	gt_averaged = []
	while not os.path.isfile('gt_averaged.dat'):
		print("waiting")
	with open('gt_averaged.dat', 'r') as f:
		for line in f.readlines():
			if len(line.split()) > 0 and line.split()[0] != '#':
				gt_averaged.append([float(line.split()[0]), float(line.split()[1])])
	irheo(gt_averaged)
	return


##################################
# 
def irheo(data_list):
	minmax = [1e-6, 1e2]
	div = 8
	#
	gw = calcgw(data_list, minmax, div)
	save_gw_data(gw, 'gw.dat')
	#
	plotgtgw('gw.dat')
	#
	return

def calcgw(gt, minmax, div):
	gw = []
	mag = math.log10(minmax[0])
	while mag < math.log10(minmax[1]):
		for i in range(div):
			omega = 10**(mag+i/div)
			gstar = gs(gt, omega)
			gw.append([omega, gstar.real, abs(gstar.imag)])
		mag += 1
	#
	return gw

def gs(gt, omega):
	gstar = gt[0][1] + (1 - cmath.exp(-1j*omega*gt[1][0]))*(gt[1][1] - gt[0][1])/gt[1][0]/(1j*omega)
	for k in range(len(gt) - 2):
		gstar += (gt[k+2][1] - gt[k+1][1])*(cmath.exp(-1j*omega*gt[k+1][0]) - cmath.exp(-1j*omega*gt[k+2][0]))/(gt[k+2][0] - gt[k+1][0])/(1j*omega)
	#
	return gstar 

#----- 計算結果をターゲットファイル名で保存
def save_gw_data(target, f_data):
	with open(f_data,'w') as f:
		for line in target:
			for data in line:
				f.write(str(data) + '\t')
			f.write('\n')
	return

#----- 結果をプロット
def plotgtgw(f_data):
	plt = make_gtgw(f_data)
	#
	if platform.system() == "Windows":
		subprocess.call([plt], shell=True)
	elif platform.system() == "Linux":
		subprocess.call(['gnuplot ' + plt], shell=True)
	return

# 必要なスクリプトを作成
def make_gtgw(f_data):
	script = gtgw_content(f_data)
	plt = f_data.replace('dat', 'plt')
	with open(plt, 'w') as f:
		f.write(script)
	return plt

# スクリプトの中身
def gtgw_content(f_data):
	out_png = f_data.replace('dat', 'png')
	script = 'set term pngcairo font "Arial,14"\n\n'
	script += 'set colorsequence classic\n\n'
	script += 'data = "' + f_data + '"\n\n'
	script += 'set output "' + out_png + '"\n\n'
	script += 'set key left\nset size square\n'
	script += 'set xrange [:1e2]\nset yrange [1e-4:]\nset y2range [1e-1:1e1]\nset y2tics\n'
	script += 'set logscale xyy2\n'
	script += '# 斜辺の傾きが -2 の三角形の準備\n'
	script += 'a = 10; # グラフの中に入るように三角形の高さを調整\n'
	script += 'x1=5e-5; x2=1e-3;\n'
	script += 'y1=a*x1**(1);y2=a*x2**(1);\n'
	script += 'set object 1 polygon from x1,y1 to x2,y1 to x2,y2 to x1,y1 fs empty border\n\n'
	script += 'set format x "10^{%L}" \nset format y "10^{%L}"\nset format y2 "10^{%L}"\n'
	# if var.step_def_mode == 'stretch':
	# 	script += 'set label 1 sprintf("{/Symbol l} = %.2f", ' + str(var.step_strain) + ') at graph 0.5, 0.9\n\n'
	# else:
	# 	script += 'set label 1 sprintf("{/Symbol g} = %.2f", ' + str(var.step_strain) + ') at graph 0.5, 0.9\n\n'
	script += 'set xlabel "Frequency"\nset ylabel "G' + "', G''" + '"\nset y2label "tan{/Symbol d}"\n\n'
	script += 'plot	data u 1:2 w lp lt 1 ti "G' + "'" + '", \\\n'
	script += 'data u 1:3 w lp lt 2 ti "G' + "''" + '", \\\n'
	script += 'data u 1:($3/$2) axis x1y2 w lp lt 3 ti "tan{/Symbol d}"'
	script += '\n\nreset'

	return script