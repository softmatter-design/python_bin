#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
import numpy as np
import platform
from UDFManager import UDFManager
import os
##########################################
# UDF の作成
##########################################
class SetUpUDF:
	def __init__(self, nw_type, files_cond, target_name, target_dir, names):
		# 
		self.calc_type = nw_type
		#
		self.ver_cognac = files_cond[0]
		self.blank_udf = files_cond[1]
		self.base_udf = files_cond[2]
		self.core = ' -n ' + str(files_cond[3])
		self.shebang = files_cond[4]
		#
		self.target_name = target_name
		# 
		self.f_eval_py = 'evaluate_all.py'
		# ファイルのディレクトリ
		self.target_dir = target_dir
		# Cognac用の名称設定
		self.nw_name = names[0]
		self.atom_name = names[1]
		self.bond_name = names[2]
		self.angle_name = names[3]
		self.site_name = names[4]
		self.pair_name = names[5]
		self.site_pair_name = names[6]

	####################################
	# UDFファイルを設定し、バッチ処理を作成
	def setup_udf(self):
		if platform.system() == "Windows":
			batch = ""
		elif platform.system() == "Linux":
			batch = "#!/bin/bash\n"

		###############
		# nw_typeに応じて計算条件を選択
		if self.calc_type == "homo_KG":
			print("making homo KG")
			batch = self.homo_kg(batch)
			# 評価用のパイソンスクリプトを作成
			self.evaluate_setup("chain")
		elif self.calc_type == "homo_sunuke":
			print("making homo sunuke")
			batch = self.homo_sunuke(batch)
			# 評価用のパイソンスクリプトを作成
			self.evaluate_setup("chain")
		elif self.calc_type == "KG_entangled" or self.calc_type == "KG_gel" or self.calc_type == "KG_multi":
			batch = self.kg_calc(batch)
			# 評価用のパイソンスクリプトを作成
			self.evaluate_setup("strand")
		elif self.calc_type == "Sunuke":
			batch = self.sunuke_calc(batch)
			# 評価用のパイソンスクリプトを作成
			self.evaluate_setup("strand")
		elif self.calc_type == "KG_NPT":
			batch = self.npt_calc(batch)
			# 評価用のパイソンスクリプトを作成
			self.evaluate_setup("strand")

		#####################
		# バッチファイルを作成
		f_batch = os.path.join(self.target_dir, '_Calc_all.bat')
		with open(f_batch,'w') as f:
			# f.write(batch_all)
			f.write(batch)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)
		return

	# 評価用のスクリプトを記述したファイルを作成
	def evaluate_setup(self, target):
		script = self.shebang + '\n# -*- coding: utf-8 -*-\n'
		script += '################################\n'
		script += 'import os \nimport sys \n'
		# if target == "chain":
		script += 'from evaluation import EvaluateChain\n'
		script += '################################\n'
		script += 'ec = EvaluateChain.EvaluateAll()\n'
		script += 'ec.evaluate_all()\n\n'
		# elif target == "strand":
		# 	script += 'from evaluation import EvaluateStrand\n'
		# 	script += '################################\n'
		# 	script += 'EvaluateStrand.evaluate_all()\n\n'
		#
		f_script = os.path.join(self.target_dir, self.f_eval_py)
		with open(f_script, 'w') as f:
			f.write(script)
		return

	#####################################################################################
	# ターミナルのタイトルを設定
	def make_title(self, batch, title):
		if platform.system() == "Windows":
			batch += "title " + title + "\n"
#		elif platform.system() == "Linux":
#			batch += r'echo -ne "\033]0; ' + title + ' \007"' + '\n'
		return batch

	# ファイル名の処理
	def make_step(self, time, fn_ext, batch, f_eval):
		present_udf = fn_ext[0] + self.target_name + fn_ext[1]
		out_udf = present_udf.replace("uin", "out")
		batch += self.ver_cognac + ' -I ' + present_udf + ' -O ' + out_udf + self.core + ' \n'
		if f_eval:
			batch += 'python ' + self.f_eval_py + ' ' + out_udf + '\n'
		read_udf = out_udf
		return present_udf, read_udf, batch

######################
# 各種バッチ条件を設定
######################
	######################################################################
	# ホモポリマーのKG鎖の計算
	def homo_kg(self, batch):
		# Force Capped LJ によりステップワイズに初期化
		r = 0.9558*2**(1/6)
		batch = self.make_title(batch, "Calculating-Init")
		fn_ext = ['Init_', '_uin.udf']
		time = [0.001, 1000000, 10000]
		f_eval = 1
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.step_nonbond_setup(self.base_udf, 'random', present_udf, time, r)
		pre = read_udf
		template = present_udf
		#
		for r in [1.0, 0.9, 0.8]:
			# 平衡化
			batch = self.make_title(batch, "Calculating-Pre_" + str(round(r, 3)).replace('.', '_'))
			fn_ext = ['Pre_rfc_' + str(round(r, 3)).replace('.', '_') + "_", "_uin.udf"]
			time = [0.01, 5000000, 50000]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.step_nonbond_setup(template, pre, present_udf, time, r)
			pre = read_udf
			template = present_udf
		# KG 鎖に設定
		time = [0.01, 10000000, 100000]
		batch = self.make_title(batch, "Calculating-KG")
		fn_ext = ['KG_', "_uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.kg_setup(template, pre, present_udf, time)
		pre = read_udf
		template = present_udf
		# 平衡化計算
		repeat = 4
		time = [0.01, 2000000, 5000]
		for i in range(repeat):
			# 平衡化
			batch = self.make_title(batch, "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "_uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.eq_setup(template, pre, present_udf, time)
			pre = read_udf
			template = present_udf
		# グリーン久保
		repeat = 5
		time = [0.01, 2000000, 100000]
		for i in range(repeat):
			# 平衡化
			batch = self.make_title(batch, "Calculating-GK_" + str(i))
			fn_ext = ['GK_' + str(i) + "_", "_uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.greenkubo_setup(template, pre, present_udf, time)
			pre = read_udf
			template = present_udf
		return batch

	######################################################################
	# ホモポリマーのsunuke鎖の計算
	def homo_sunuke(self, batch):
		# Force Capped LJ によりステップワイズに初期化
		r = 0.9558*2**(1/6)
		batch = self.make_title(batch, "Calculating-Init")
		fn_ext = ['Init_', '_uin.udf']
		time = [0.002, 1000000, 10000]
		f_eval = 0
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.step_nonbond_setup(self.base_udf, 'random', present_udf, time, r)
		pre = read_udf
		template = present_udf
		# sunuke 鎖に設定
		time = [0.01, 1000000, 10000]
		batch = self.make_title(batch, "Calculating-KG")
		fn_ext = ['Sunuke_', "_uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.sunuke_setup(template, pre, present_udf, time)
		pre = read_udf
		template = present_udf
		# 平衡化計算
		repeat = 4
		time = [0.01, 1000000, 5000]
		for i in range(repeat):
			# 平衡化
			batch = self.make_title(batch, "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "_uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.eq_setup(template, pre, present_udf, time)
			pre = read_udf
			template = present_udf
		# # グリーン久保
		# repeat = 5
		# time = [0.01, 2000000, 100000]
		# for i in range(repeat):
		# 	# 平衡化
		# 	batch = self.make_title(batch, "Calculating-GK_" + str(i))
		# 	fn_ext = ['GK_' + str(i) + "_", "_uin.udf"]
		#	f_eval = 1
		# 	present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		# 	self.greenkubo_setup(template, pre, present_udf, time)
		# 	pre = read_udf
			template = present_udf
		return batch

	######################################################################
	# KG鎖の計算
	def kg_calc(self, batch):
		# Force Capped LJ によりステップワイズに初期化
		r = 1.1
		batch = self.make_title(batch, "Calculating-Init")
		fn_ext = ['Init_', "_uin.udf"]
		time = [0.001, 500000, 5000]
		f_eval = 1
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.step_nonbond_setup(self.base_udf, '', present_udf, time, r)
		pre = read_udf
		template = present_udf
		#
		for r in [0.9558*2**(1/6), 1.0, 0.9, 0.8]:
			# 平衡化
			batch = self.make_title(batch, "Calculating-Pre_" + str(round(r, 3)).replace('.', '_'))
			fn_ext = ['Pre_rfc_' + str(round(r, 3)).replace('.', '_') + "_", "_uin.udf"]
			time = [0.01, 1000000, 10000]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.step_nonbond_setup(template, pre, present_udf, time, r)
			pre = read_udf
			template = present_udf
		# KG 鎖に設定
		time = [0.01, 1000000, 10000]
		batch = self.make_title(batch, "Calculating-KG")
		fn_ext = ['KG_', "_uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.kg_setup(template, pre, present_udf, time)
		pre = read_udf
		template = present_udf
		# 平衡化計算
		time = [0.01, 500000, 1000]
		for i in [1,2]:
			# 平衡化
			batch = self.make_title(batch, "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "_uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.eq_setup(template, pre, present_udf, time)
			pre = read_udf
			template = present_udf
		# グリーン久保
		repeat = 5
		time = [0.01, 20000000, 1000000]
		for i in range(repeat):
			# 平衡化
			batch = self.make_title(batch, "Calculating-GK_" + str(i))
			fn_ext = ['GK_' + str(i) + "_", "_uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.greenkubo_setup(template, pre, present_udf, time)
			pre = read_udf
			template = present_udf
		return batch

	###############################################################
	# スヌケ鎖の計算
	def sunuke_calc(self, batch):
		# Force Capped LJ により初期化
		r = 1.1
		batch = self.make_title(batch, "Calculating-Init")
		fn_ext = ['Init_', '_uin.udf']
		time = [0.001, 100000, 1000]
		f_eval = 1
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.step_nonbond_setup(self.base_udf, '', present_udf, time, r)
		pre = read_udf
		template = present_udf

		# sunuke 鎖に設定
		time = [0.01, 200000, 2000]
		batch = self.make_title(batch, "Calculating-sunuke")
		fn_ext = ['Set_', "_uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.sunuke_setup(template, pre, present_udf, time)
		pre = read_udf
		template = present_udf
		# 平衡化計算
		repeat = 4
		time = [0.01, 1000000, 5000]
		for i in range(repeat):
			# 平衡化
			batch = self.make_title(batch, "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "_uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.eq_setup(template, pre, present_udf, time)
			pre = read_udf
			template = present_udf
		# # グリーン久保
		# repeat = 5
		# time = [0.01, 2000000, 100000]
		# for i in range(repeat):
		# 	# 平衡化
		# 	batch = self.make_title(batch, "Calculating-GK_" + str(i))
		# 	fn_ext = ['GK_' + str(i) + "_", "_uin.udf"]
		#	f_eval = 1
		# 	present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		# 	self.greenkubo_setup(template, pre, present_udf, time)
		# 	pre = read_udf
			template = present_udf
		return batch

	###########################################
	# NPT 条件で、設定密度まで圧縮
	def npt_calc(self, batch):
		# NPTの設定
		pres = 0.01
		batch = self.make_title(batch, "Calculating-Ini_NPT")
		fn_ext = ['Init_pres_' + str(pres).replace('.', '_') + '_', "_uin.udf"]
		time = [0.001, 1000, 100]
		f_eval = 0
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.npt_setup(self.base_udf, '', present_udf, time, pres)
		pre = read_udf
		template = present_udf
		# ステップワイズに圧力増加
		for pres in [0.02, 0.05, 0.1, 0.5, 0.75]:
			batch = self.make_title(batch, "Calculating-Step_NPT")
			fn_ext = ['Compress_pres_' + str(pres).replace('.', '_') + '_', "_uin.udf"]
			time = [0.01, 500000, 20000]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.npt_setup(template, pre, present_udf, time, pres)
			pre = read_udf
			template = present_udf
		# ステップワイズに圧力増加
		for pres in [1.0, 1.3, 1.6, 2.0]:
			batch = self.make_title(batch, "Calculating-Step_NPT")
			fn_ext = ['Compress_pres_' + str(pres).replace('.', '_') + '_', "_uin.udf"]
			time = [0.01, 2000000, 50000]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.npt_setup(template, pre, present_udf, time, pres)
			pre = read_udf
			template = present_udf
		for pres in [2.5, 3.0, 4.5]:
			batch = self.make_title(batch, "Calculating-Step_NPT")
			fn_ext = ['Compress_pres_' + str(pres).replace('.', '_') + '_', "_uin.udf"]
			time = [0.01, 500000, 20000]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.npt_setup(template, pre, present_udf, time, pres)
			pre = read_udf
			template = present_udf
		# KG 鎖に設定
		time = [0.01, 100000, 1000]
		batch = self.make_title(batch, "Calculating-KG")
		fn_ext = ['Setup_', "_uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		self.kg_setup(template, pre, present_udf, time)
		pre = read_udf
		template = present_udf
		# 平衡化計算
		time = [0.01, 200000, 1000]
		for i in range(5):
			# 平衡化
			batch = self.make_title(batch, "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "_uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
			self.eq_setup(template, pre, present_udf, time)
			pre = read_udf
			template = present_udf
		# # グリーン久保
		# repeat = 3
		# time = [0.01, 2000000, 100000]
		# for i in range(repeat):
		# 	# 平衡化
		# 	batch = self.make_title(batch, "Calculating-GK_" + str(i))
		# 	fn_ext = ['GK_' + str(i) + "_", "_uin.udf"]
		#	f_eval = 1
		# 	present_udf, read_udf, batch = self.make_step(time, fn_ext, batch, f_eval)
		# 	self.greenkubo_setup(template, pre, present_udf, time)
		# 	pre = read_udf
		# 	template = present_udf
		return batch

###############
# UDF 作成条件
###############
	##############################################################################
	# ユーザーポテンシャルにより、Force Capped LJ で、ステップワイズにノンボンドを増加
	def step_nonbond_setup(self, template, read_udf, present_udf, time, r_fc):
		u = UDFManager(os.path.join(self.target_dir, template))
		# goto global data
		u.jump(-1)
		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000000., 	p + 'Max_Force')
		u.put(time[0], 		p + 'Time.delta_T')
		u.put(time[1], 		p + 'Time.Total_Steps')
		u.put(time[2], 		p + 'Time.Output_Interval_Steps')
		u.put(1.0, 			p + 'Temperature.Temperature')
		u.put(0., 			p + 'Pressure_Stress.Pressure')
		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(1, p + 'Bond')
		u.put(1, p + 'Angle')
		if read_udf == 'random':
			u.put(0, p + 'Non_Bonding_Interchain')
			u.put(0, p + 'Non_Bonding_1_3')
			u.put(0, p + 'Non_Bonding_1_4')
			u.put(0, p + 'Non_Bonding_Intrachain')
		else:
			u.put(1, p + 'Non_Bonding_Interchain')
			u.put(1, p + 'Non_Bonding_1_3')
			u.put(1, p + 'Non_Bonding_1_4')
			u.put(1, p + 'Non_Bonding_Intrachain')
		#--- Initial_Structure ---
		# Initial_Unit_Cell
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0.85, p + 'Density')
		u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
		# Generate_Method
		if read_udf == 'random':
			# Set atoms
			p = 'Initial_Structure.Generate_Method.'
			u.put('Random', p + 'Method')
			u.put(1, p + 'Random.Fix_Angle')
			# Relaxation
			p = 'Initial_Structure.Relaxation.'
			u.put(1, p + 'Relaxation')
			u.put('DYNAMICS', p + 'Method')
			u.put(100, p + 'Max_Relax_Force')
			u.put(100000, p + 'Max_Relax_Steps')
		elif read_udf == '':
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', p+'Method')
			u.put(['', -1, 0, 0], p+'Restart')
			# Relaxation
			p = 'Initial_Structure.Relaxation.'
			u.put(1, p + 'Relaxation')
			u.put('DYNAMICS', p + 'Method')
			u.put(100, p + 'Max_Relax_Force')
			u.put(100000, p + 'Max_Relax_Steps')
		else:
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', p+'Method')
			u.put([read_udf, -1, 1, 0], p+'Restart')
			# Relaxation
			p = 'Initial_Structure.Relaxation.'
			u.put(0, p + 'Relaxation')
		#--- Simulation_Conditions ---
		# bond
		for i, b_name in enumerate(self.bond_name):
			p = 'Molecular_Attributes.Bond_Potential[].'
			u.put(b_name, 		p + 'Name', [i])
			if read_udf == '':
				u.put('Harmonic', 	p + 'Potential_Type', [i])
				u.put(0.97,			p + 'R0', [i])
				u.put(1000, 		p + 'Harmonic.K', [i])
			else:
				u.put('User_Bond', 	p + 'Potential_Type', [i])
				u.put(0.9609,		p + 'R0', [i])
				u.put(1,			p + 'User_Bond.Index', [i])
				#
				u.put('sigma',		p + 'User_Bond.Parameters[].Name', [i, 0])
				u.put('epsilon',	p + 'User_Bond.Parameters[].Name', [i, 1])
				u.put('c0',			p + 'User_Bond.Parameters[].Name', [i, 2])
				u.put('c1',			p + 'User_Bond.Parameters[].Name', [i, 3])
				u.put('c2',			p + 'User_Bond.Parameters[].Name', [i, 4])
				u.put('c3',			p + 'User_Bond.Parameters[].Name', [i, 5])
				#
				u.put('1.0',		p + 'User_Bond.Parameters[].Value', [i, 0])
				u.put('1.0',		p + 'User_Bond.Parameters[].Value', [i, 1])
				u.put('20.2026',	p + 'User_Bond.Parameters[].Value', [i, 2])
				u.put('490.628',	p + 'User_Bond.Parameters[].Value', [i, 3])
				u.put('2256.76',	p + 'User_Bond.Parameters[].Value', [i, 4])
				u.put('9685.31',	p + 'User_Bond.Parameters[].Value', [i, 5])
		# Angle
		for i, anglename in enumerate(self.angle_name):
			p = 'Molecular_Attributes.Angle_Potential[].'
			u.put(anglename, 		p + 'Name', [i])
			u.put('User_Angle', 	p + 'Potential_Type', [i])
			u.put(73, 				p + 'theta0', [i])
			u.put(1, 				p + 'User_Angle.Index', [i])
			#
			u.put('sigma', 			p + 'User_Angle.Parameters[].Name', [i, 0])
			u.put('epsilon', 		p + 'User_Angle.Parameters[].Name', [i, 1])
			u.put('rfc', 			p + 'User_Angle.Parameters[].Name', [i, 2])
			#
			u.put('1.0', 			p + 'User_Angle.Parameters[].Value', [i, 0])
			u.put('1.0', 			p + 'User_Angle.Parameters[].Value', [i, 1])
			u.put('0.8', 			p + 'User_Angle.Parameters[].Value', [i, 2])
		#--- Pair_Interaction[] ---
		for i, pairname in enumerate(self.pair_name):
			p = 'Interactions.Pair_Interaction[].'
			u.put(pairname,   					p + 'Name', [i])
			u.put('User_Pair_Interaction', 		p + 'Potential_Type', [i])
			u.put(self.site_pair_name[i][0],	p + 'Site1_Name', [i])
			u.put(self.site_pair_name[i][1],	p + 'Site2_Name', [i])
			u.put(1.12246204830937,				p + 'Cutoff', [i])
			u.put(1.0,							p + 'Scale_1_4_Pair', [i])
			u.put(1,							p + 'User_Pair_Interaction.Index', [i])
			#
			u.put('sigma',						p + 'User_Pair_Interaction.Parameters[].Name', [i, 0])
			u.put('epsilon',					p + 'User_Pair_Interaction.Parameters[].Name', [i, 1])
			u.put('rfc',						p + 'User_Pair_Interaction.Parameters[].Name', [i, 2])
			#
			u.put('1.0',						p + 'User_Pair_Interaction.Parameters[].Value', [i, 0])
			u.put('1.0',						p + 'User_Pair_Interaction.Parameters[].Value', [i, 1])
			u.put(str(r_fc),						p + 'User_Pair_Interaction.Parameters[].Value', [i, 2])
		#--- Write UDF ---
		u.write(os.path.join(self.target_dir, present_udf))
		return

	###############################################
	# ボンドをFENE、ノンボンドをLJとしてKG鎖を設定
	def kg_setup(self, template, read_udf, present_udf, time):
		u = UDFManager(os.path.join(self.target_dir, template))
		u.jump(-1)
		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time[0],  p+'Time.delta_T')
		u.put(time[1],  p+'Time.Total_Steps')
		u.put(time[2],  p+'Time.Output_Interval_Steps')
		u.put(1.0, 		p + 'Temperature.Temperature')
		u.put(0., 		p + 'Pressure_Stress.Pressure')
		# Solver
		p = 'Simulation_Conditions.Solver.'
		u.put('Dynamics', 			p + 'Solver_Type')
		u.put('NVT_Kremer_Grest', 	p + 'Dynamics.Dynamics_Algorithm')
		u.put(0.5, 					p + 'Dynamics.NVT_Kremer_Grest.Friction')
		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(1, p + 'Bond')
		u.put(0, p + 'Angle')
		u.put(1, p + 'Non_Bonding_Interchain')
		u.put(1, p + 'Non_Bonding_1_3')
		u.put(1, p + 'Non_Bonding_1_4')
		u.put(1, p + 'Non_Bonding_Intrachain')
		#--- Initial_Structure ---
		# Initial_Unit_Cell
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0.85, p + 'Density')
		u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', p+'Method')
		u.put([read_udf, -1, 0, 0], p+'Restart')
		p = 'Initial_Structure.Relaxation.'
		u.put(1, p + 'Relaxation')
		#--- Simulation_Conditions ---
		# Bond
		for i, b_name in enumerate(self.bond_name):
			p = 'Molecular_Attributes.Bond_Potential[].'
			u.put(b_name, 		p + 'Name', [i])
			u.put('FENE_LJ', 	p + 'Potential_Type', [i])
			u.put(1.0,			p + 'R0', [i])
			u.put(1.5,			p + 'FENE_LJ.R_max', [i])
			u.put(30,			p + 'FENE_LJ.K', [i])
			u.put(1.0,			p + 'FENE_LJ.sigma', [i])
			u.put(1.0,			p + 'FENE_LJ.epsilon', [i])
		# Site
		for i, sitename in enumerate(self.site_name):
			p = 'Molecular_Attributes.Interaction_Site_Type[].'
			u.put(sitename, 	p + 'Name', [i])
			u.put(1, 			p + 'Num_of_Atoms', [i])
			u.put(1.0, 			p + 'Range', [i])
		#--- Pair_Interaction[] ---
		for i, pairname in enumerate(self.pair_name):
			p = 'Interactions.Pair_Interaction[].'
			u.put(pairname,   					p + 'Name', [i])
			u.put('Lennard_Jones', 				p + 'Potential_Type', [i])
			u.put(self.site_pair_name[i][0],	p + 'Site1_Name', [i])
			u.put(self.site_pair_name[i][1],	p + 'Site2_Name', [i])
			u.put(2**(1/6),						p + 'Cutoff', [i])
			u.put(1.0,							p + 'Scale_1_4_Pair', [i])
			u.put(1.0,							p + 'Lennard_Jones.sigma', [i])
			u.put(1.0,							p + 'Lennard_Jones.epsilon', [i])
		#--- Write UDF ---
		u.write(os.path.join(self.target_dir, present_udf))
		return

	####################################
	# スヌケ鎖をボンドをFENEで設定
	def sunuke_setup(self, template, read_udf, present_udf, time):
		u = UDFManager(os.path.join(self.target_dir, template))
		# goto global data
		u.jump(-1)
		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000000., 	p + 'Max_Force')
		u.put(time[0], 		p + 'Time.delta_T')
		u.put(time[1], 		p + 'Time.Total_Steps')
		u.put(time[2], 		p + 'Time.Output_Interval_Steps')
		u.put(1.0, 			p + 'Temperature.Temperature')
		u.put(0., 			p + 'Pressure_Stress.Pressure')
		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(1, p + 'Bond')
		u.put(1, p + 'Angle')
		u.put(0, p + 'Non_Bonding_Interchain')
		u.put(0, p + 'Non_Bonding_1_3')
		u.put(0, p + 'Non_Bonding_1_4')
		u.put(0, p + 'Non_Bonding_Intrachain')
		#--- Initial_Structure ---
		# Initial_Unit_Cell
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0.85, p + 'Density')
		u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
		# Generate_Method
		if read_udf == '':
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', 		p + 'Method')
			u.put(['', -1, 0, 0], 	p + 'Restart')
		else:
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', p+'Method')
			u.put([read_udf, -1, 1, 0], p+'Restart')
		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		u.put(1, p + 'Relaxation')
		u.put('DYNAMICS', p + 'Method')
		u.put(200, p + 'Max_Relax_Force')
		u.put(10000, p + 'Max_Relax_Steps')
		#--- Simulation_Conditions ---
		# Bond
		p = 'Molecular_Attributes.Bond_Potential[].'		
		for i, bondname in enumerate(self.bond_name):
			u.put(bondname, 	p + 'Name', [i])
			if read_udf == '':
				u.put('Harmonic', 	p + 'Potential_Type', [i])
				u.put(0.97,			p + 'R0', [i])
				u.put(1000, 		p + 'Harmonic.K', [i])
			else:
				u.put('FENE_LJ', 	p + 'Potential_Type', [i])
				u.put(1.0,			p + 'R0', [i])
				u.put(1.5,			p + 'FENE_LJ.R_max', [i])
				u.put(30,			p + 'FENE_LJ.K', [i])
				u.put(1.0,			p + 'FENE_LJ.sigma', [i])
				u.put(1.0,			p + 'FENE_LJ.epsilon', [i])
		#--- Write UDF ---
		u.write(os.path.join(self.target_dir, present_udf))
		return

	###################################
	# NPT 条件をハーモニックボンドで設定
	def npt_setup(self, template, read_udf, present_udf, time, pressure):
		u = UDFManager(os.path.join(self.target_dir, template))
		u.jump(-1)
		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000000.,	p + 'Max_Force')
		u.put(time[0],		p + 'Time.delta_T')
		u.put(time[1], 		p + 'Time.Total_Steps')
		u.put(time[2], 		p + 'Time.Output_Interval_Steps')
		u.put(1.0, 			p + 'Temperature.Temperature')
		# Pressure
		u.put(pressure, 	p + 'Pressure_Stress.Pressure')
		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(1, p + 'Bond')
		u.put(0, p + 'Angle')
		u.put(1, p + 'Non_Bonding_Interchain')
		u.put(1, p + 'Non_Bonding_1_3')
		u.put(1, p + 'Non_Bonding_1_4')
		u.put(1, p + 'Non_Bonding_Intrachain')
		#--- Initial_Structure ---
		if read_udf == '':
			# Generate_Method
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', 		p + 'Method')
			u.put(['', -1, 0, 0], 	p + 'Restart')
		else:
			# Initial_Unit_Cell
			p = 'Initial_Structure.Initial_Unit_Cell.'
			u.put(0, p + 'Density')
			u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
			# Generate_Method
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', 			p + 'Method')
			u.put([read_udf, -1, 1, 0], p + 'Restart')
		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		u.put(1, p + 'Relaxation')
		u.put('DYNAMICS', p + 'Method')
		u.put(200, p + 'Max_Relax_Force')
		u.put(10000, p + 'Max_Relax_Steps')
		#--- Simulation_Conditions ---
		# Bond
		p = 'Molecular_Attributes.Bond_Potential[].'		
		for i, bondname in enumerate(self.bond_name):
			u.put(bondname, 	p + 'Name', [i])
			if read_udf == '':
				u.put('Harmonic', 	p + 'Potential_Type', [i])
				u.put(0.97,			p + 'R0', [i])
				u.put(1000, 		p + 'Harmonic.K', [i])
			else:
				u.put('FENE_LJ', 	p + 'Potential_Type', [i])
				u.put(1.0,			p + 'R0', [i])
				u.put(1.5,			p + 'FENE_LJ.R_max', [i])
				u.put(30,			p + 'FENE_LJ.K', [i])
				u.put(1.0,			p + 'FENE_LJ.sigma', [i])
				u.put(1.0,			p + 'FENE_LJ.epsilon', [i])
		# #--- Write UDF ---
		u.write(os.path.join(self.target_dir, present_udf))
		return

	#####################################################
	# 直前の条件を維持して平衡化
	def eq_setup(self, template, read_udf, present_udf, time):
		u = UDFManager(os.path.join(self.target_dir, template))
		u.jump(-1)
		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time[0],  p+'Time.delta_T')
		u.put(time[1],  p+'Time.Total_Steps')
		u.put(time[2],  p+'Time.Output_Interval_Steps')
		# Moment
		u.put(0, p + "Moment.Interval_of_Calc_Moment")
		u.put(0, p + "Moment.Calc_Moment")
		u.put(0, p + "Moment.Stop_Translation")
		u.put(0, p + "Moment.Stop_Rotation")

		#--- Initial_Structure ---
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', p+'Method')
		u.put([read_udf, -1, 1, 0], p+'Restart')
		p = 'Initial_Structure.Relaxation.'
		u.put(0, p + 'Relaxation')

		#--- Write UDF ---
		u.write(os.path.join(self.target_dir, present_udf))
		return

	###########################################################
	# グリーン久保計算
	def greenkubo_setup(self, template, read_udf, present_udf, time):
		u = UDFManager(os.path.join(self.target_dir, template))
		u.jump(-1)
		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time[0],  p+'Time.delta_T')
		u.put(time[1],  p+'Time.Total_Steps')
		u.put(time[2],  p+'Time.Output_Interval_Steps')
		# Calc Correlation
		u.put(1, 'Simulation_Conditions.Output_Flags.Correlation_Function.Stress')
		#--- Initial_Structure ---
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', p+'Method')
		u.put([read_udf, -1, 1, 0], p+'Restart')
		p = 'Initial_Structure.Relaxation.'
		u.put(0, p + 'Relaxation')
		#--- Write UDF ---
		u.write(os.path.join(self.target_dir, present_udf))
		return
