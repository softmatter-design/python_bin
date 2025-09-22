#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *

import os
import platform

import mod_deform_setup.variables as var
import mod_global.glob_var as globvar
import mod_deform_setup.DeformGeneral as gen

#######
# 繰り返し変形の設定
def setup_cyclic_deform():
	set_cyclic_basedir()
	set_each_cycle()
	return

def set_cyclic_basedir():
	var.cyc_dir = var.cyclic_deform + '_read_' + var.read_udf.split('.')[0]
	if os.path.exists(var.cyc_dir):
		print("Use existing dir of ", var.cyc_dir)
	else:
		print("Make new dir of ", var.cyc_dir)
		os.makedirs(var.cyc_dir)
	return
#
def set_each_cycle():
	for cyc_def_max in var.cyc_deform_cond_dic:
		for cyc_rate in var.cyc_deform_cond_dic[cyc_def_max][1]:
			cond_dir = 'Deform_until_' + str(cyc_def_max).replace('.', '_') + "_rate_" + f"{cyc_rate:.1e}".replace('.', '_')
			var.cyc_dirlist.append(cond_dir)
			middle_dir = os.path.join(var.cyc_dir, cond_dir)
			if os.path.exists(middle_dir):
				print("Use existing dir of ", middle_dir)
			else:
				print("Make new dir of ", middle_dir)
				os.makedirs(middle_dir)
			set_cyclic_rotation(middle_dir, cyc_def_max, cyc_rate)
	# プラットフォームに応じて命令を変更
	if platform.system() == "Windows":
		task = 'call calc_series.bat\n'
		filename = 'calc_all.bat'
	elif platform.system() == "Linux":
		task = 'sh calc_series.sh\n'
		filename = 'calc_all.sh'
	gen.make_batch_series(var.cyc_dirlist, var.cyc_dir, task, filename, '')
	return

def set_cyclic_rotation(middle_dir, cyc_def_max, cyc_rate):
	if var.cyclic_deform == 'CyclicShear':
		var.cyc_rotate = ['base', 'x', 'y', 'z', 'yx', 'zx']
		var.sim_deform = 'shear'
	elif var.cyclic_deform == 'CyclicStretch':
		var.cyc_rotate = ['base', 'x', 'y']
		var.sim_deform = 'stretch'
	for rotate in var.cyc_rotate:
		set_cyc_rotate_dir(middle_dir, rotate, cyc_def_max, cyc_rate)
	# プラットフォームに応じて命令を変更
	if platform.system() == "Windows":
		task = 'call calc.bat\n'
		filename = 'calc_series.bat'
		path = os.path.join(globvar.bin_path, 'eval_cyc_def.py')
		option = f'python {str(path):} -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} -a \n'
	elif platform.system() == "Linux":
		task = 'sh calc.sh\n'
		filename = 'calc_all.sh'
		option = f'eval_cyc_def.py -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} -a \n'
	gen.make_batch_series(['rotate_' + dir for dir in var.cyc_rotate], middle_dir, task, filename, option)
	return

def set_cyc_rotate_dir(middle_dir, rotate, cyc_def_max, cyc_rate):
	tmp_dir = f'rotate_{rotate}'
	var.calc_dir = os.path.join(middle_dir, tmp_dir)
	if os.path.exists(var.calc_dir):
		print("Use existing dir of ", var.calc_dir)
	else:
		print("Make new dir of ", var.calc_dir)
		os.makedirs(var.calc_dir)
	var.base_udf = os.path.join(var.calc_dir, 'base.udf')
	u = UDFManager(var.read_udf)
	u.jump(1)
	u.eraseRecord(record_pos=0, record_num=u.totalRecord()-1)
	if rotate != 'base':
		gen.rotate_position(u, rotate)
	u.write(var.base_udf)

	make_cycle_batch(cyc_def_max, cyc_rate)

	return

# ファイル名を設定し、バッチファイルを作成
def make_cycle_batch(cyc_def_max, cyc_rate):
	repeatcount = ''
	var.batch = "#!/bin/bash\n"
	for var.cyc_count in range(var.cyc_deform_cond_dic[cyc_def_max][0]):
		var.cyc_resol = var.cyc_deform_cond_dic[cyc_def_max][2]
		make_cycle(cyc_def_max, cyc_rate)
		repeatcount += str(var.cyc_count) + ' '
	# 評価スクリプトを追加
	if platform.system() == "Windows":
		path = os.path.join(globvar.bin_path, 'eval_cyc_def.py')
		var.batch += f'python {str(path):} -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:}\n'
	elif platform.system() == "Linux":
		var.batch += f'eval_cyc_def.py -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:}\n'
	# バッチファイルを作成
	gen.write_batchfile(var.calc_dir, 'calc.bat', var.batch)
	return
#
def make_cycle(cyc_def_max, cyc_rate):
	for var.cyc_direction in ['_Forward', '_Backward']:
		# UDFファイル名を設定
		uin = 'No_' +str(var.cyc_count) + var.cyc_direction + "_uin.udf"
		uout = uin.replace("uin", "out")
		gen.make_title(var.title_name + "_Calculating_Cycle_until_" + str(cyc_def_max).replace('.', '_') + "_rate_" + f"{cyc_rate:.1e}".replace('.','_') + '_No' + str(var.cyc_count) + var.cyc_direction)
		var.batch += globvar.ver_Cognac + ' -I ' + uin + ' -O ' + uout + ' -n ' + str(var.core) +' \n'
		udf_in =  os.path.join(var.calc_dir, uin)
		if var.cyc_count == 0 and var.cyc_direction == '_Forward':
			var.cyc_readudf = 'base.udf'
		mod_cycle_udf(cyc_def_max, cyc_rate, udf_in)
		var.cyc_readudf = uout
	return

#
def mod_cycle_udf(cyc_def_max, cyc_rate, udf_in):
	if var.cyclic_deform == 'CyclicStretch':
		deform_time = (cyc_def_max - 1)/cyc_rate
		speed = cyc_rate*var.system_size
	elif var.cyclic_deform == 'CyclicShear':
		deform_time = cyc_def_max/cyc_rate
	#
	time_total = round(deform_time/var.sim_time_div)
	time_1_step = round(var.cyc_resol/var.sim_time_div/cyc_rate)
	#
	u = UDFManager(var.base_udf)
	# goto global data
	u.jump(-1)
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(100000.,		p + 'Max_Force')
	u.put(var.sim_time_div,	p + 'Time.delta_T')
	u.put(time_total,	p + 'Time.Total_Steps')
	u.put(time_1_step,	p + 'Time.Output_Interval_Steps')
	u.put(1.0,			p + 'Temperature.Temperature')
	u.put(0, 			p + 'Temperature.Interval_of_Scale_Temp')
	u.put(0,			p + 'Pressure_Stress.Pressure')
	# Deformation
	if var.cyclic_deform == 'CyclicStretch':
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('Cell_Deformation', 	p + 'Method')
		u.put('Simple_Elongation', 	p + 'Cell_Deformation.Method')
		u.put('Deformation_Speed', 	p + 'Cell_Deformation.Simple_Elongation.Input_Method')
		if var.cyc_direction == '_Forward':
			u.put(speed, p + 'Cell_Deformation.Simple_Elongation.Deformation_Speed.Speed')
		else:
			u.put(-1.*speed, p + 'Cell_Deformation.Simple_Elongation.Deformation_Speed.Speed')
		u.put(0.5, 						p + 'Cell_Deformation.Simple_Elongation.Poisson_Ratio')
		u.put('z', 						p + 'Cell_Deformation.Simple_Elongation.Axis')
		u.put(1, 						p + 'Cell_Deformation.Interval_of_Deform')
		u.put(0, 						p + 'Cell_Deformation.Deform_Atom')
	elif var.cyclic_deform == 'CyclicShear':
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('Lees_Edwards', 	p + 'Method')
		u.put('Steady', 		p + 'Lees_Edwards.Method')
		if var.cyc_direction == '_Forward':
			u.put(cyc_rate, 		p + 'Lees_Edwards.Steady.Shear_Rate')
		else:
			u.put(-1.*cyc_rate, 	p + 'Lees_Edwards.Steady.Shear_Rate')
	# Output_Flags
	u.put([1, 1, 1], 'Simulation_Conditions.Output_Flags.Structure')
	# Read_Set_of_Molecules
	p = 'Initial_Structure.Read_Set_of_Molecules'
	u.put(['', -1], p)
	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', 		p + 'Method')
	u.put([var.cyc_readudf, -1, 1, 1], 	p + 'Restart')
	# Relaxation
	p = 'Initial_Structure.Relaxation.'
	u.put(0, p + 'Relaxation')
	#--- Write UDF ---
	u.write(udf_in)
	return