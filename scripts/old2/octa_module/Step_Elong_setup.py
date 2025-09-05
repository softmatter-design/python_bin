#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *
import sys
import os
import shutil
import numpy as np
import platform
import pandas as pd
################################################################################
def main():
	print("このスクリプトを直接読んでも、初期状態が入っていません。")
	exit(1)
################################################################################
class Setup:
	def __init__(self, py_mod, ver_cognac, calc_dir, f_data, core, stress_eval, gt_eval, cond_list):
		self.py_mod = py_mod
		self.ver = ver_cognac
		self.calc_dir = calc_dir
		self.f_data = f_data
		self.core = ' -n ' + str(core) +' \n'
		self.stress_eval = stress_eval
		self.gt_eval = gt_eval
		#
		self.step_rate = cond_list[0]
		self.step_time_div = cond_list[1]
		self.lamd_max = cond_list[2]
		self.res = cond_list[3]
		self.quench_list = cond_list[4]
		#
		self.base_udf = "base_uin.udf"

	############################################################################
	##### Main #####
	# 
	def calc_step(self):
		# 平衡計算したUDFファイルとその場所を選択
		t_udf, data_dir = self.file_select()
		#
		func, nu, structure = self.read_data(data_dir)
		#
		base = self.make_base_udf(t_udf)
		#
		self.make_batch(base)
		#
		self.make_script_ss(func, nu, structure)
		#
		self.make_script_gt(func, nu, structure)
		return

	############################################################################
	##### Function #####
	# 平衡計算したUDFファイルとその場所を選択
	def file_select(self):
		param = sys.argv
		if len(param) == 1:
			print("usage: python", param[0], "Honya_out.udf")
			exit(1)
		elif not os.access(param[1],os.R_OK):
			print(param[1], "not exists.")
			exit(1)
		else:
			target = param[1]
			data_dir = os.path.dirname(target)
		return target, data_dir

	# 
	def read_data(self, data_dir):
		t_data = os.path.join(data_dir, self.f_data)
		if not os.path.exists(t_data):
			func = 0
			nu = 0
			structure = 'none'
		else:
			with open(t_data, 'r') as f:
				calc_cond = f.readlines()[-1].strip().split('\t')
			func = calc_cond[3]
			nu = calc_cond[4]
			structure = calc_cond[5]
		return func, nu, structure

	# 
	def make_base_udf(self, t_udf):
		if os.path.exists(self.calc_dir):
			print("Use existing dir of ", self.calc_dir)
		else:
			print("Make new dir of ", self.calc_dir)
			os.makedirs(self.calc_dir)
		#
		base = os.path.join(self.calc_dir, self.base_udf)
		print("Readin file = ", t_udf)
		u = UDFManager(t_udf)
		u.jump(1)
		u.eraseRecord(record_pos=-999,record_num=-999)
		u.write(base)
		return base

	# ファイル名を設定し、バッチファイルを作成
	def make_batch(self, base):
		batch = "#!/bin/bash\n\n"
		# ステップ伸長
		rate_str = "{0:4.0e}".format(self.step_rate)
		lambda_str = str(self.lamd_max)
		uin = 'Step_Elong_lambda_' + lambda_str + '_rate_' + rate_str + "_uin.udf"
		batch = self.make_title(batch, "Calculating Step rate=" + rate_str)
		read_udf, batch = self.make_step(uin, batch)
		self.step_elong(base, uin)
		pre = read_udf
		template = uin
		# ステップ伸長後にSSカーブ及び温度をプロット
		batch += 'python ' + self.stress_eval + '\n'
		# 
		count = 0
		for time in self.quench_list:
			present_udf = 'Quench_from_' + lambda_str + "_step_" + str(time[0]*time[2]).replace(".", "_") + "_uin.udf"
			batch = self.make_title(batch, "Calculating Quench")
			read_udf, batch = self.make_step(present_udf, batch)
			self.eq_udf(template, pre, present_udf, time)
			pre = read_udf
			template = present_udf
			# ステップ伸長後にSSカーブ及び温度をプロット
			batch += 'python ' + self.gt_eval + '\n'
		
		
		# バッチファイルを作成
		f_batch = os.path.join(self.calc_dir, '_step_elong.bat')
		with open(f_batch, 'w') as f:
			f.write(batch)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)
		return

	###########################
	# ターミナルのタイトルを設定
	def make_title(self, batch, title):
		if platform.system() == "Windows":
			batch += "title " + title + "\n"
		elif platform.system() == "Linux":
			batch += r'echo -ne "\033]0; ' + title + ' \007"' + '\n'
		return batch

	# ファイル名の処理
	def make_step(self, present_udf, batch):
		out_udf = present_udf.replace("uin", "out")
		batch += self.ver + ' -I ' + present_udf + ' -O ' + out_udf + self.core
		read_udf = out_udf
		return read_udf, batch
		

	#-----
	def step_elong(self, base, udf_in):
		time_1_step = int(self.res/self.step_time_div/self.step_rate)
		time_total = time_1_step*(self.lamd_max - 1)/self.res
		#
		u = UDFManager(base)
		u.jump(-1)

		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000.,				p + 'Max_Force')
		u.put(self.step_time_div,	p + 'Time.delta_T')
		u.put(time_total,			p + 'Time.Total_Steps')
		u.put(time_1_step,			p + 'Time.Output_Interval_Steps')
		u.put(1.0,					p + 'Temperature.Temperature')
		u.put(0, 					p + 'Temperature.Interval_of_Scale_Temp')
		u.put(0,					p + 'Pressure_Stress.Pressure')

		# Deformation
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('Cell_Deformation', 		p + 'Method')
		u.put('Simple_Elongation', 		p + 'Cell_Deformation.Method')
		u.put('Initial_Strain_Rate', 	p + 'Cell_Deformation.Simple_Elongation.Input_Method')
		u.put(self.step_rate,	 		p + 'Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate')
		u.put(0.5, 						p + 'Cell_Deformation.Simple_Elongation.Poisson_Ratio')
		u.put('z', 						p + 'Cell_Deformation.Simple_Elongation.Axis')
		u.put(1, 						p + 'Cell_Deformation.Interval_of_Deform')
		u.put(0, 						p + 'Cell_Deformation.Deform_Atom')

		# Read_Set_of_Molecules
		p = 'Initial_Structure.Read_Set_of_Molecules'
		u.put([self.base_udf, -1], p)

		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 					p + 'Method')
		u.put([self.base_udf, -1, 1, 0], 	p + 'Restart')

		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		u.put(0, p + 'Relaxation')

		#--- Write UDF ---
		u.write(os.path.join(self.calc_dir, udf_in))
		return
	
	################################################################################
	def eq_udf(self, template, read_udf, present_udf, time):
		u = UDFManager(os.path.join(self.calc_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time[0],  p+'Time.delta_T')
		u.put(time[1],  p+'Time.Total_Steps')
		u.put(time[2],  p+'Time.Output_Interval_Steps')

		# Deformation
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('None',	p + 'Method')

		#--- Initial_Structure ---
		p = 'Initial_Structure.Read_Set_of_Molecules'
		u.put([read_udf, -1], p)
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 			p+'Method')
		u.put([read_udf, -1, 1, 0], p+'Restart')
		p = 'Initial_Structure.Relaxation.'
		u.put(0, p + 'Relaxation')

		#--- Write UDF ---
		u.write(os.path.join(self.calc_dir, present_udf))
		return

	#######################
	# 必要なスクリプトを作成
	def make_script_ss(self, func, nu, structure):
		script = self.script_content_ss(func, nu, structure)
		with open(os.path.join(self.calc_dir, self.stress_eval),'w') as f:
			f.write(script)
		return

	# スクリプトの中身
	def script_content_ss(self, func, nu, structure):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n################################\n'
		script += 'import os \nimport sys \n'
		script += 'py_mod = "' + self.py_mod + '"\n'
		script += 'sys.path.append(py_mod)\n'
		script += 'from network_evaluation import Read_Step_Stress\n################################\n'
		script += 'func = ' + str(func) + '\nnu = ' + str(nu) + '\n################################\n'
		script += 'rs = Read_Step_Stress.MakeSS()\n'
		script += 't_udf_list = rs.file_select()\n'
		script += 'stress_data = rs.calc_stress_all(t_udf_list)\n'
		script += 'rs.save_data(stress_data, nu, t_udf_list)\n'
		script += 'rs.plot(func, nu)\n'
		return script

	#######################
	# 必要なスクリプトを作成
	def make_script_gt(self, func, nu, structure):
		script = self.script_content_gt(func, nu, structure)
		with open(os.path.join(self.calc_dir, self.gt_eval),'w') as f:
			f.write(script)
		return

	# スクリプトの中身
	def script_content_gt(self, func, nu, structure):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n'
		script += '################################\n'
		script += 'import os \nimport sys \n'
		script += 'py_mod = "' + self.py_mod + '"\n'
		script += 'sys.path.append(py_mod)\n'
		script += 'from network_evaluation import Read_Step_Stress\n'
		script += 'import platform\nimport subprocess\n'
		script += '################################\n'
		script += 'func = ' + str(func) + '\nnu = ' + str(nu) + '\nlamd = ' + str(self.lamd_max) + '\nstructure = "' + structure + '"\n'
		script += '################################\n'
		script += 'qc = Read_Step_Stress.QuenchCalc()\n'
		script += 't_udf_list = qc.file_select()\n'
		script += 'series = qc.calc_gt(t_udf_list, lamd)\n'
		script += 'qc.save_data(series)\n'
		script += 'cmd = "corr2gw < gt_series.dat > gw.dat"\n'
		script += 'subprocess.call(cmd, shell=True)\n'
		script += 'qc.plot(func, lamd, nu, structure)\n\n'
		return script

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()