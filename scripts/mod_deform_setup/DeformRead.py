#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *

import argparse
import codecs
import numpy as np
import os
import platform
import sys
import uuid

import mod_deform_setup.variables as var
import mod_global.glob_var as globvar
################
################
# 各種条件を読み取り
def read_all():
	read_arg()
	read_nw_cond()
	read_sim_cond()
	return

def read_arg():
	parser = argparse.ArgumentParser(description='Select udf file to read !')
	parser.add_argument('udf', help="udf file name to read previous simulation")
	args = parser.parse_args()
	if args.udf:
		if len(args.udf.split('.')) != 2 or args.udf.split('.')[1] != 'udf':
			print('\nthe file name you selected is not udf file !')
			sys.exit('select proper udf file to read.')
		elif not os.access(args.udf, os.R_OK):
			sys.exit('\nSelected udf of ', args.udf, ' seems not exist !\nbye now!!')
		else:
			var.read_udf = args.udf
			var.base_name = args.udf.split('.')[0]
	else:
		print('no udf file is selected')
		sys.exit('select proper udf file to read.')
	return

# 計算対象の条件を読み取る
def read_nw_cond():
	if not os.access(os.path.join(var.base_name, 'target_condition.udf'), os.R_OK):
		sys.exit("\n'target_condition.udf' is not exists in the '" + var.base_name + "' directory!!")
		return
	else:
		cond_u = UDFManager(os.path.join(var.base_name, 'target_condition.udf'))
		var.func = cond_u.get('TargetCond.NetWork.Function')
		var.nu = cond_u.get('TargetCond.System.Nu')
		var.system_size = cond_u.get('TargetCond.System.SystemSize')
		var.title_base = os.path.basename(os.path.dirname(os.path.abspath(var.read_udf)))
	return

# シミュレーション条件を設定する。
def read_sim_cond():
	if os.path.isfile('deform_condition.udf'):
		var.deform_udf = 'deform_condition.udf'
		read_and_set()
	else:
		while not os.path.isfile('../deform_condition.udf'):
			print('\nIn the parent directory, no "deform_condition.udf" is found !')
			print('New one will be generated.')
			print('Please, modify and save it !\n')
			make_newudf()
			input('Press ENTER to continue...')
		var.deform_udf = '../deform_condition.udf'
		read_and_set()
	return

# make new udf when not found.
def make_newudf():
	contents = '''
	\\begin{def}
	CalcConditions:{
		Cognac_ver:select{"cognac112"} "使用する Cognac のバージョン",
		Cores: int "計算に使用するコア数を指定"
		} "Cognac による計算の条件を設定"
	SimpleDeformation:{
		DeformMode:select{"none", "Stretch", "Shear", "both"} "変形モードを選択",
			Stretch:{
				DeformRate[]:float "これらは変形レートのリスト",
				MaxDeformation:float "最大ひずみ",
				Resolution:float "これは１ステップ計算での伸長度 Res = lambda/1_step"
				}
			Shear:{
				DeformRate[]:float "これらは変形レートのリスト",
				MaxDeformation:float "最大ひずみ",
				Resolution:float "これは１ステップ計算での伸長度 Res = lambda/1_step"
				}
			both:{
				DeformRate[]:float "これらは変形レートのリスト",
				MaxDeformation:float "最大ひずみ",
				Resolution:float "これは１ステップ計算での伸長度 Res = lambda/1_step"
				}
		} "計算ターゲットの条件を設定"		
	CycleDeformation:{
		CyclicDeform:select{"none", "CyclicStretch", "CyclicShear"} "変形モードを選択",
		CyclicStretch:{
			StretchConditions[]:{
				MaxDeformation:float "最大ひずみ",
				Repeat:int "サイクルの繰り返し数",
				DeformRate[]:float "これらは変形レートのリスト",
				Resolution:float "これは１ステップ計算での伸長度 Res = lambda/1_step"
				}
			}
		CyclicShear:{
			ShearConditions[]:{
				MaxDeformation:float "最大ひずみ",
				Repeat:int "サイクルの繰り返し数",
				DeformRate[]:float "これらは変形レートのリスト",
				Resolution:float "これは１ステップ計算での伸長度 Res = lambda/1_step"
				}
			}
		} "計算ターゲットの条件を設定"
	StepDeformation:{
		StepDeform:select{"none", "StepStretch", "StepShear"} "変形モードを選択",
		StepStretch:{
			StretchConditions:{
				Deformation:{
					MaxDeformation:float "最大ひずみ",
					DeformRate:float "変形レート",
					DeformSteps:int "シミュレーションのステップ数"
					}
				Relaxation[]:{
					RelaxationTime:int "緩和を観測する時間",
					CalcSteps:int "緩和時間の分割数",
					}
				}
			}
		StepShear:{
			ShearConditions:{
				Deformation:{
					MaxDeformation:float "最大ひずみ",
					DeformRate:float "変形レート",
					DeformSteps:int "シミュレーションのステップ数"
					}
				Relaxation[]:{
					RelaxationTime:int "緩和を観測する時間",
					CalcSteps:int "緩和時間の分割数",
					}
				}
			}
		} "計算ターゲットの条件を設定"
	\end{def}	

	\\begin{data}
	CalcConditions:{"cognac112",1}
	SimpleDeformation:{
		"both",
			{
			[1.0e-03,5.0e-4,1.0e-04,5.0e-05]
			3.00,
			1.0e-02
			}
			{
			[1.0e-03,5.0e-4,1.0e-04,5.0e-05]
			2.0,
			1.0e-02
			}
			{
			[1.0e-03,5.0e-4,1.0e-04,5.0e-05]
			3.00,
			1.0e-02
			}
		}
	CycleDeformation:{
		"CyclicShear",
			{
				[
					{1.0,
					3,
					[1.0e-03,5.0e-04,1.0e-04,5.0e-05,1.0e-05]
					1.0e-02
					}
					{2.0,
					3,
					[1.0e-03,1.00e-04,1.00e-05]
					1.00e-02
					}
				]
			}
			{
				[
					{2.00,
					3,
					[1.0e-03,1.00e-04]
					1.00e-02
					}
					{3.00,
					3,
					[1.00e-03,1.00e-04,1.00e-05]
					1.00e-02
					}
				]
			}
		}
	StepDeformation:{
		"StepStretch",
			{
				{
				{1.50,5.0e-02,200}
				[{10,200}{100,200}{1000,200}{10000,200}{100000,200}]
				}
			}
			{
				{
				{1.00,0.10,200}
				[{10,200}{100,200}{1000,200}{10000,200}{100000,200}]
				}
			}
		}
	\end{data}
	'''
	###
	with codecs.open('../deform_condition.udf', 'w', 'utf_8') as f:
		f.write(contents)
	return

# Read udf and setup initial conditions
def read_and_set():
	dic={'y':True,'yes':True,'q':False,'quit':False}
	while True:
		# read udf
		read_condition()
		# select
		init_calc()
		print('Change UDF: type [r]eload')
		print('Quit input process: type [q]uit')
		inp = input('Condition is OK ==> [y]es >> ').lower()
		if inp in dic:
			inp = dic[inp]
			break
		print('##### \nRead Condition UDF again \n#####\n\n')
	if inp:
		return
	else:
		sys.exit("##### \nQuit !!")

# Read condition udf
def read_condition():
	u = UDFManager(var.deform_udf)
	u.jump(-1)
	# 使用するCognacのバージョン
	globvar.ver_Cognac = u.get('CalcConditions.Cognac_ver')
	# 計算に使用するコア数
	var.core = u.get('CalcConditions.Cores')
	# Simple Deformation
	var.simple_def_mode  = u.get('SimpleDeformation.DeformMode')
	if var.simple_def_mode == 'Stretch':
		var.sim_rate_list = u.get('SimpleDeformation.Stretch.DeformRate[]')
		var.sim_deform_max = u.get('SimpleDeformation.Stretch.MaxDeformation')
		var.sim_resolution = u.get('SimpleDeformation.Stretch.Resolution')
		var.sim_deform = var.simple_def_mode
	elif var.simple_def_mode == 'Shear':
		var.sim_rate_list = u.get('SimpleDeformation.Shear.DeformRate[]')
		var.sim_deform_max = u.get('SimpleDeformation.Shear.MaxDeformation')
		var.sim_resolution = u.get('SimpleDeformation.Shear.Resolution')
		var.sim_deform = var.simple_def_mode
	elif var.simple_def_mode == 'both':
		var.sim_rate_list = u.get('SimpleDeformation.both.DeformRate[]')
		var.sim_deform_max = u.get('SimpleDeformation.both.MaxDeformation')
		var.sim_resolution = u.get('SimpleDeformation.both.Resolution')
	# Cyclic Deformation
	tmp = []
	var.cyclic_deform = u.get('CycleDeformation.CyclicDeform')
	if var.cyclic_deform == 'CyclicStretch':
		tmp = u.get('CycleDeformation.CyclicStretch.StretchConditions[]')
	elif var.cyclic_deform == 'CyclicShear':
		tmp = u.get('CycleDeformation.CyclicShear.ShearConditions[]')
	for data in tmp:
		max_strain, repeat, ratelist, resolution = data
		var.cyc_deform_cond_dic[max_strain] = [repeat, ratelist, resolution]
	# Step Deformation
	var.step_deform = u.get('StepDeformation.StepDeform')
	if var.step_deform == 'StepShear':
		[var.step_deform_max, var.step_rate, var.step_steps] = u.get('StepDeformation.StepShear.ShearConditions.Deformation')
		deform_time = var.step_deform_max/var.step_rate
		#
		var.step_relaxation = u.get('StepDeformation.StepShear.ShearConditions.Relaxation[]')
	elif var.step_deform == 'StepStretch':
		[var.step_deform_max, var.step_rate, var.step_steps] = u.get('StepDeformation.StepStretch.StretchConditions.Deformation')
		if var.step_deform_max == 1.0:
			sys.exit('\nStep Stretch Condition is not proper !!\nMax Deformation should be greater than 1.0 !')
		else:
			deform_time = abs(var.step_deform_max - 1)/var.step_rate
		#
		var.step_relaxation = u.get('StepDeformation.StepStretch.StretchConditions.Relaxation[]')
	#
	if var.step_deform != 'none':
		dt = min(var.sim_time_div, deform_time/var.step_steps)	# dt の暫定値を決定
		total_steps = round(deform_time/dt)
		interval = max(1, round(total_steps/var.step_steps))	# 整数値の interval を決定
		dt = round(deform_time/var.step_steps/interval, 4)		# 小数点４桁で丸めたdtを決定
		var.step_deform_time = [dt, total_steps, interval]
	#
	if var.simple_def_mode == 'none' and var.cyclic_deform == 'none' and var.step_deform == 'none':
		sys.exit('No proper condition is selected.\nBye!')
	return
# 
def init_calc():
	text = "################################################" + "\n"
	text += "Cores used for simulation\t\t" + str(var.core ) + "\n"
	text += "################################################" + "\n"
	if var.simple_def_mode != 'none':
		text += "Deform mode:\t\t\t\t" + str(var.simple_def_mode) + "\n"
		text += "Deform Rate:\t\t" + ', '.join([f"{x:.1e}" for x in var.sim_rate_list]) + "\n"
		text += "Maximum Strain:\t\t\t\t" + str(var.sim_deform_max) + "\n"
		text += "Resolution:\t\t\t\t" + str(round(var.sim_resolution,4)) + "\n"
		text += "################################################" + "\n"
	if var.cyclic_deform != 'none':
		text += "Cyclic Deform mode:\t\t\t" + str(var.cyclic_deform) + "\n"
		count = 0
		for key in var.cyc_deform_cond_dic:
			text += f'Cyclic condition #{count}\n'
			text += f"\tMaximum Strain:\t\t\t{key:.1f}\n"
			text += f"\tRepeat:\t\t\t\t{var.cyc_deform_cond_dic[key][0]}\n"
			text += "\tCyclic Deform Rate:\t" + ', '.join([f"{x:.1e}" for x in var.cyc_deform_cond_dic[key][1]]) + "\n"
			text += "\tResolution:\t\t\t" + str(round(var.cyc_deform_cond_dic[key][2], 4)) + "\n"
		text += "################################################" + "\n"
	if var.step_deform != 'none':
		text += f"Step Deform mode:\t\t\t{var.step_deform:}\n"
		text += f"Step Strain:\t\t\t\t{var.step_deform_max:.1f}\n"
		text += f"Deformation rate:\t\t\t{var.step_rate:.1e}\n"
		text += f"Deformation steps:\t\t\t{var.step_steps:}\n"
		text += f"Simulation time:\t\t{var.step_deform_time:}\n"
		text += "#\n"
		for i, data in enumerate(var.step_relaxation):
			text += f"Relaxation-{i:}\n"
			text += f"\tRelaxation Time:\t\t{data[0]:.1e}\n"
			text += f"\tCalc. steps:\t\t\t{data[1]:}\n"
		text += "################################################" + "\n"
	print(text)
	return