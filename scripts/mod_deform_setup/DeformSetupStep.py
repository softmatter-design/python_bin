#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *

# import numpy as np
import os
import platform
# import uuid

import mod_global.glob_var as globvar
import mod_deform_setup.variables as var
import mod_deform_setup.DeformGeneral as gen

# #######################################
# ステップ変形の設定
def setup_step_deform():
	# 計算用のディレクトリーを作成
	set_step_basedir()
	# 
	set_rotation_step()
	
	return
	
# 計算用のディレクトリーを作成
def set_step_basedir():
	var.step_dir = f"Deform_read_{var.base_name:}/{var.step_deform:}_until_" + f'{var.step_deform_max:.1f}'.replace('.','_') + '_rate_' + f'{var.step_rate:.1e}'.replace('.', '_')
	if os.path.exists(var.step_dir):
		print("Use existing dir of ", var.step_dir)
	else:
		print("Make new dir of ", var.step_dir)
		os.makedirs(var.step_dir)
	return

def set_rotation_step():
	if var.step_deform == 'StepShear':
		var.step_rotate = ['base', 'x', 'y', 'z', 'yx', 'zx']
		var.sim_deform = 'shear'
	elif var.step_deform == 'StepStretch':
		var.step_rotate = ['base', 'x', 'y']
		var.sim_deform = 'stretch'
	for rotate in var.step_rotate:
		set_rotate_dir(rotate)
		set_udf_batch(rotate)
	#
	if platform.system() == "Windows":
		task = 'call calc_series.bat\n'
		filename = 'calc_all.bat'
		path = os.path.join(globvar.bin_path, 'eval_step_def.py')
		option = f'python {str(path):} -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} -a \n'
	elif platform.system() == "Linux":
		task = 'sh calc_series.sh\n'
		filename = 'calc_all.sh'
		option = f'eval_step_def.py -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} -a \n'
	gen.make_batch_series(var.step_dirlist, var.step_dir, task, filename, option)
	return

def set_rotate_dir(rotate):
	tmp_dir = f'rotate_{rotate}'
	var.step_dirlist.append(tmp_dir)
	var.calc_dir = os.path.join(var.step_dir, tmp_dir)
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
	return

def set_udf_batch(rotate):
	# UDFファイル名を設定
	base = f'{var.step_deform}_until_' + f'{var.step_deform_max:.1e}'.replace('.', '_') + '_rate_' + f'{var.step_rate:.1e}'.replace('.', '_') + f'_{rotate}'
	uin = 'deform_uin.udf'
	uout = uin.replace("uin", "out")
	gen.make_title(var.title_base + '_' + base + "_deform")
	var.batch = "#!/bin/bash\n"
	var.batch += globvar.ver_Cognac + ' -I ' + uin + ' -O ' + uout + ' -n ' + str(var.core) +' \n'
	udf_in =  os.path.join(var.calc_dir, uin)
	make_stepdeform_udf(udf_in)
	prev_udf = uout

	# 各放置時間における緩和計算を設定
	for i, condition in enumerate(var.step_relaxation):
		uin = f'relaxation_{i}_uin.udf'
		uout = uin.replace("uin", "out")
		gen.make_title(var.title_base + '_' + base + f'_relaxation_{i}')
		var.batch += globvar.ver_Cognac + ' -I ' + uin + ' -O ' + uout + ' -n ' + str(var.core) +' \n'
		udf_in =  os.path.join(var.calc_dir, uin)
		make_steprelax_udf(udf_in, prev_udf, condition)
	
	# 最長の緩和計算のUDFをリスタートにして長時間計算を繰り返す。
	repeat = var.step_repeat[0]
	for i in range(repeat):
		condition = var.step_repeat[1:]
		uin = f'repeat_{i}_uin.udf'
		udf_in =  os.path.join(var.calc_dir, uin)
		uout = uin.replace("uin", "out")
		gen.make_title(var.title_base + '_' + base + f'_repeat_{i}')
		var.batch += globvar.ver_Cognac + ' -I ' + uin + ' -O ' + uout + ' -n ' + str(var.core) +' \n'
		make_steprelax_udf(udf_in, prev_udf, condition)
		prev_udf = uin.replace("uin", "out")
	#
	if platform.system() == "Windows":
		filename = 'calc_series.bat'
		path = os.path.join(globvar.bin_path, 'eval_step_def.py')
		var.batch += f'python {str(path):} -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} \n'
	elif platform.system() == "Linux":
		filename = 'calc_series.sh'
		var.batch += f'eval_step_def.py -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} \n'
	gen.write_batchfile(var.calc_dir, filename, var.batch)
	return

#-----
def make_stepdeform_udf(udf_in):
	u = UDFManager(var.base_udf)
	# goto global data
	u.jump(-1)

	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(100000.,		p + 'Max_Force')
	u.put(var.step_deform_time[0],	p + 'Time.delta_T')
	u.put(var.step_deform_time[1],	p + 'Time.Total_Steps')
	u.put(var.step_deform_time[2],	p + 'Time.Output_Interval_Steps')
	u.put(1.0,			p + 'Temperature.Temperature')
	u.put(0, 			p + 'Temperature.Interval_of_Scale_Temp')
	u.put(0,			p + 'Pressure_Stress.Pressure')

	# Deformation
	if var.step_deform == 'StepStretch':
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('Cell_Deformation', 		p + 'Method')
		u.put('Simple_Elongation', 		p + 'Cell_Deformation.Method')
		u.put('Initial_Strain_Rate', 	p + 'Cell_Deformation.Simple_Elongation.Input_Method')
		u.put(var.step_rate,	 					p + 'Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate')
		if var.step_deform_max < 1:
			u.put(-1*var.step_rate,	 					p + 'Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate')
		u.put(0.5, 						p + 'Cell_Deformation.Simple_Elongation.Poisson_Ratio')
		u.put('z', 						p + 'Cell_Deformation.Simple_Elongation.Axis')
		u.put(1, 						p + 'Cell_Deformation.Interval_of_Deform')
		u.put(0, 						p + 'Cell_Deformation.Deform_Atom')
	elif var.step_deform == 'StepShear':
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('Lees_Edwards', 	p + 'Method')
		u.put('Steady', 		p + 'Lees_Edwards.Method')
		u.put(var.step_rate, 			p + 'Lees_Edwards.Steady.Shear_Rate')
	
	# Output_Flags
	u.put([1, 1, 1], 'Simulation_Conditions.Output_Flags.Structure')

	# Read_Set_of_Molecules
	p = 'Initial_Structure.Read_Set_of_Molecules'
	u.put(['', -1], p)

	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', 		p + 'Method')
	u.put(['', -1, 1, 1], 	p + 'Restart')

	# Relaxation
	p = 'Initial_Structure.Relaxation.'
	u.put(0, p + 'Relaxation')
	#--- Write UDF ---
	u.write(udf_in)
	return

#-----
def make_steprelax_udf(udf_in, prev_udf, condition):
	u = UDFManager(var.base_udf)
	# goto global data
	u.jump(-1)

	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(100000.,		p + 'Max_Force')
	u.put(var.sim_time_div,	p + 'Time.delta_T')
	u.put(round(condition[0]/var.sim_time_div),	p + 'Time.Total_Steps')
	u.put(round(condition[0]/var.sim_time_div/condition[1]),	p + 'Time.Output_Interval_Steps')
	u.put(1.0,			p + 'Temperature.Temperature')
	u.put(0, 			p + 'Temperature.Interval_of_Scale_Temp')
	u.put(0,			p + 'Pressure_Stress.Pressure')

	# Output_Flags
	u.put([1, 1, 1], 'Simulation_Conditions.Output_Flags.Structure')
	
	# Read_Set_of_Molecules
	p = 'Initial_Structure.Read_Set_of_Molecules'
	u.put([prev_udf, -1], p)

	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', 		p + 'Method')
	u.put([prev_udf, -1, 1, 1], 	p + 'Restart')

	# Relaxation
	p = 'Initial_Structure.Relaxation.'
	u.put(0, p + 'Relaxation')
	#--- Write UDF ---
	u.write(udf_in)
	return