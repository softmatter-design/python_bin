#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *

import os
import platform

import mod_deform_setup.variables as var
import mod_global.glob_var as globvar
import mod_deform_setup.DeformGeneral as gen

######################################
# 単純変形の設定
def setup_simple_deform():
	if var.simple_def_mode == 'both':
		for var.sim_deform in ['shear', 'stretch']:
			set_simple_eachrate()
	else:
		set_simple_eachrate()
	return

def set_simple_eachrate():
	var.sim_basedir = f"{var.sim_deform:}_calculation_read_{var.read_udf.split('.')[0]:}_until_{var.sim_deform_max:}"
	if os.path.exists(var.sim_basedir):
		print("Use existing dir of ", var.sim_basedir)
	else:
		print("Make new dir of ", var.sim_basedir)
		os.makedirs(var.sim_basedir)

	c_dir = os.getcwd().split('\\')[-1]
	var.title_base = str(c_dir.split('_', 2)[-1]) + f"_{var.sim_deform:}_calculation_until_{var.sim_deform_max:}_"
	# プラットフォームに応じて命令を変更
	if platform.system() == "Windows":
		task = 'call calc_series.bat\n'
		filename = 'calc_all.bat'
	elif platform.system() == "Linux":
		task = 'sh calc_series.sh\n'
		filename = 'calc_all.sh'
		#
		# task2 = 'sh eval_all.sh\n'
		# filename2 = 'eval_all.sh'
		# option = f'simple_deform.py -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} -s \n'
		# gen.make_batch_series([f'rate_{rate:4.0e}' for rate in var.sim_rate_list], var.sim_basedir, task2, filename2, option)
	gen.make_batch_series([f'rate_{rate:4.0e}' for rate in var.sim_rate_list], var.sim_basedir, task, filename,'')

	for var.sim_rate in var.sim_rate_list:
		var.sim_ratedir = os.path.join(var.sim_basedir, f"rate_{var.sim_rate:4.0e}")
		#
		if os.path.exists(var.sim_ratedir):
			print("Use existing dir of ", var.sim_ratedir)
		else:
			print("Make new dir of ", var.sim_ratedir)
			os.makedirs(var.sim_ratedir)
		#
		set_rotation_simple()
	return

def set_rotation_simple():
	# 変形方法に応じて回転方向を設定
	if var.sim_deform == 'shear':
		var.step_rotate = ['base', 'x', 'y', 'z', 'yx', 'zx']
	elif var.sim_deform == 'stretch':
		var.step_rotate = ['base', 'x', 'y']
	# プラットフォームに応じて命令を変更
	if platform.system() == "Windows":
		task = 'call calc.bat\n'
		filename = 'calc_series.bat'
		path = os.path.join(globvar.bin_path, 'eval_sim_def.py')
		option = f'python {str(path):} -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} -a \n'
	elif platform.system() == "Linux":
		task = 'sh calc.sh\n'
		filename = 'calc_series.sh'
		option = f'eval_sim_def.py -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} -a \n'
	gen.make_batch_series(['rotate_' + dir for dir in var.step_rotate], var.sim_ratedir, task, filename, option)

	for rotate in var.step_rotate:
		set_rotate_dir_sim(rotate)
		set_udf_batch_sim()
	return

def set_rotate_dir_sim(rotate):
	tmp_dir = f'rotate_{rotate}'
	var.title_name = var.title_base + f"rate_{var.sim_rate:4.0e}" + f'_rotate_{rotate}'
	var.sim_dirlist.append(tmp_dir)
	var.calc_dir = os.path.join(var.sim_ratedir, tmp_dir)
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

# ファイル名を設定し、バッチファイルを作成
def set_udf_batch_sim():
	# UDFファイル名を設定
	uin = f'rate_{var.sim_rate:4.0e}_uin.udf'
	# プラットフォームに応じてバッチファイルを設定
	var.batch = "#!/bin/bash\n"
	gen.make_title(var.title_name)
	var.batch += globvar.ver_Cognac + ' -I ' + uin + ' -O ' + uin.replace("uin", "out") + ' -n ' + str(var.core) +' \n'
	if platform.system() == "Windows":
		target_bat = 'calc.bat'
		path = os.path.join(globvar.bin_path, 'eval_sim_def.py')
		var.batch += f'python {str(path):} -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} \n'
	elif platform.system() == "Linux":
		target_bat = 'calc.sh'
		var.batch += f'eval_sim_def.py -f {str(var.func):} -n {str(var.nu):} -m {var.sim_deform:} \n'
	gen.write_batchfile(var.calc_dir, target_bat, var.batch)

	udf_in =  os.path.join(var.calc_dir, uin)
	make_simpledeform_udf(udf_in)
	return

# UDFを作成
def make_simpledeform_udf(udf_in):
	if var.sim_deform == 'stretch':
		deform_time = abs(var.sim_deform_max - 1)/var.sim_rate
	elif var.sim_deform == 'shear':
		deform_time = var.sim_deform_max/var.sim_rate
	#
	time_total = round(deform_time/var.sim_time_div)
	time_1_step = round(var.sim_resolution/var.sim_time_div/var.sim_rate)
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
	if var.sim_deform == 'stretch':
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('Cell_Deformation', 		p + 'Method')
		u.put('Simple_Elongation', 		p + 'Cell_Deformation.Method')
		u.put('Initial_Strain_Rate', 	p + 'Cell_Deformation.Simple_Elongation.Input_Method')
		if var.sim_deform_max < 1.:
			var.sim_rate = -1*var.sim_rate
		u.put(var.sim_rate,	 					p + 'Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate')
		u.put(0.5, 						p + 'Cell_Deformation.Simple_Elongation.Poisson_Ratio')
		u.put('z', 						p + 'Cell_Deformation.Simple_Elongation.Axis')
		u.put(1, 						p + 'Cell_Deformation.Interval_of_Deform')
		u.put(0, 						p + 'Cell_Deformation.Deform_Atom')
	elif var.sim_deform == 'shear':
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('Lees_Edwards', 	p + 'Method')
		u.put('Steady', 		p + 'Lees_Edwards.Method')
		u.put(var.sim_rate, 			p + 'Lees_Edwards.Steady.Shear_Rate')
	
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