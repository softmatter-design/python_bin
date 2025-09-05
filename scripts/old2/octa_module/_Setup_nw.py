#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
import numpy as np
import platform
from UDFManager import UDFManager
import os
#
from network_evaluation import Evaluate_Strand as es
################################################################################
def main():
	exit("このスクリプトを直接読んでも、初期状態が入っていません。")
####################################################################################
##########################################
# Initi_UDF の作成
##########################################
class SetUpUDF:
	def __init__(self, calcd_data_dic_list, files_cond, base_data, nw_type, n_strand, py_mod):
		self.calcd_data_dic_list = calcd_data_dic_list
		#
		self.ver_cognac = files_cond[0]
		self.blank_udf = files_cond[1]
		self.base_udf = files_cond[2]
		self.core = ' -n ' + str(files_cond[3])
		#
		self.nw_type = nw_type
		#
		self.target_name = base_data[0]
		self.multi_nw = base_data[1]
		self.system_size = base_data[2]
		self.a_cell = base_data[3]
		self.nu = base_data[4]
		self.structure = base_data[5]
		self.n_strand = n_strand
		# 条件設定
		# [R0, K]
		self.harmonic_low = [0.97, 100]
		self.harmonic_high = [0.97, 1000]
		# [Potential_Type, theta0, K]		
		self.angle_low = ['Theta2', 74, 2.0]
		self.angle_high = ['Theta2', 74, 10.0]
		# [Cutoff, Scale_1_4_Pair, sigma, epsilon, range]
		self.lj_cond = [2**(1/6), 1.0, 1.0, 1.0, 1.0]	
		# Cognac用の名称設定
		self.nw_name = "Network"
		self.atom_name = ["JP_A", "End_A", "Strand_A", "Side_A"]
		self.bond_name = ["bond_JP-Chn", "bond_Strand", "bond_Side"]
		self.angle_name = ["angle_AAA"]
		self.site_name = ["site_JP", "site_End", "site_Strand"]
		self.pair_name = ["site_JP-site_JP", "site_Strand-site_JP", "site_Strand-site_Strand", 
						"site_JP-site_End", "site_Strand-site_End", "site_End-site_End"]
		self.site_pair_name = [ 
						["site_JP", "site_JP"], 
						["site_Strand", "site_JP"], 
						["site_Strand", "site_Strand"],
						["site_JP", "site_End"], 
						["site_Strand", "site_End"], 
						["site_End", "site_End"]
						]
		#
		self.py_mod = py_mod
		self.f_eval_py = 'evaluate_all.py'

	################################################################################
	# UDFファイルを設定し、バッチ処理を作成
	def setup_all(self):
		# 計算用のディレクトリーを作成
		target_dir = self.make_dir()
		# base_udfを作成
		self.make_base_udf(target_dir)
		# 
		self.setup_udf(target_dir)
		return
	################################################################################
	# 計算用のディレクトリーを作成
	def make_dir(self):
		target_dir = str(self.target_name)
		os.makedirs(target_dir, exist_ok = True)
		with open(os.path.join(target_dir, "network.dat"), "w") as f:
			f.write("# func\tnu\tNW_type\n\n" + str(self.n_strand) + "\t" + str(round(self.nu, 5)) + '\t' + self.structure + "\n")
		return target_dir

	############################################
	# base_udfを作成
	def make_base_udf(self, target_dir):
		# 初期udfの内容を作成する
		self.make_base(target_dir)
		# すべてのアトムの位置座標及びボンド情報を設定
		self.setup_atoms(target_dir)
		return

	############################################
	#----- ファイル名を設定し、バッチファイルを作成
	def setup_udf(self, target_dir):
		if platform.system() == "Windows":
			batch = ""
		elif platform.system() == "Linux":
			batch = "#!/bin/bash\n"

		############################################
		# 計算の初期状態を設定
		if self.nw_type == "NoLJ_Harmonic":
			batch = self.make_title(batch, "Calculating-Pre")
			fn_ext = ['Pre_', '_NoLJ_uin.udf']
			time = [0.01, 1000000, 10000]
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
			self.setup_nolj_udf(self.base_udf, present_udf, time, self.harmonic_high, target_dir)
			pre = read_udf
			template = present_udf

		##################################################
		# 充分に絡み合わせるセットアップ、ボンドはハーモニック
		elif self.nw_type == "LJ_Harmonic" or self.nw_type == "KG_reluxed":
			# 斥力ポテンシャルの設定
			batch = self.make_title(batch, "Calculating-Pre")
			fn_ext = ['Pre_', "_LJ_Harm_uin.udf"]
			time = [0.002, 10000, 1000]
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
			self.lj_harm(self.base_udf, '', present_udf, time, self.harmonic_low, target_dir)
			pre = read_udf
			template = present_udf
			#
			for i in range(5):
				batch = self.make_title(batch, "Calculating-Rep_" + str(i))
				# スヌケでアングルポテンシャル
				fn_ext = ['Rep_' + str(i) + '_', "_noLJ_Angle_uin.udf"]
				time = [0.01, 10000, 1000]
				angle = self.angle_high
				present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
				self.ang_no_lj(template, pre, present_udf, time, angle, target_dir)
				pre = read_udf
				template = present_udf

				# 斥力でアングルポテンシャル
				fn_ext = ['Rep_' + str(i) + '_', "_LJ_Angle_uin.udf"]
				time = [0.002, 200000, 20000]
				angle = self.angle_high
				present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
				self.ang_w_lj(template, pre, present_udf, time, angle, target_dir)
				pre = read_udf
				template = present_udf

			# ボンドポテンシャルの増加
			batch = self.make_title(batch, "Calculating-Harm_high")
			fn_ext = ['Harm_high_', "_uin.udf"]
			time = [0.01, 100000, 10000]
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
			self.lj_harm(template, pre, present_udf, time, self.harmonic_high, target_dir)
			pre = read_udf
			template = present_udf

		######################################
		elif self.nw_type == "KG_simple":
			# NPTの設定
			batch = self.make_title(batch, "Calculating-NPT_low")
			fn_ext = ['Pre_', "_KG_NPT_low_uin.udf"]
			time = [0.01, 200000, 2000]
			pressure = 0.1
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
			self.kg_npt(self.base_udf, '', present_udf, time, pressure, target_dir)
			pre = read_udf
			template = present_udf

			batch = self.make_title(batch, "Calculating-NPT_high")
			fn_ext = ['Pre_', "_KG_NPT_high_uin.udf"]
			time = [0.01, 100000, 1000]
			pressure = 5.
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
			self.kg_npt(template, pre, present_udf, time, pressure, target_dir)
			pre = read_udf
			template = present_udf

			# 斥力でアングルポテンシャル
			batch = self.make_title(batch, "Calculating-NPT_angle")
			fn_ext = ['Pre_', "_KG_NPT_angle_uin.udf"]
			time = [0.01, 100000, 1000]
			angle = self.angle_low
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
			self.ang_w_lj(template, pre, present_udf, time, angle, target_dir)
			pre = read_udf
			template = present_udf

			# FENEを計算
			batch = self.make_title(batch, "Calculating-FENE_" )
			fn_ext = ['Pre_', "_KG_FENE_uin.udf"]
			time = [0.01, 100000, 10000]
			rmax = 1.5
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
			self.bond_fene(template, pre, present_udf, rmax, time, target_dir)
			pre = read_udf
			template = present_udf

		###################################
		# KG_reluxed用のセットアップ
		if self.nw_type == "KG_reluxed":
			# FENE
			rmax_list = [1.8, 1.5]
			for rmax in rmax_list:
				# FENEを計算
				batch = self.make_title(batch, "Calculating-FENE_" + str(rmax).replace('.', '_'))
				fn_ext = ['FENE_' + str(rmax).replace('.', '_') + '_', "_uin.udf"]
				time = [0.01, 100000, 10000]
				present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
				self.bond_fene(template, pre, present_udf, rmax, time, target_dir)
				pre = read_udf
				template = present_udf

		##################
		# 最終の平衡化計算
		repeat = 3
		for i in range(repeat):
			# 平衡化
			batch = self.make_title(batch, "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "_uin.udf"]
			time = [0.01, 200000, 1000]
			batch = self.make_title(batch, "Calculating-Equiv")
			present_udf, read_udf, batch = self.make_step(time, fn_ext, batch)
			self.eq_udf(template, pre, present_udf, time, target_dir)
			pre = read_udf
			template = present_udf
		
		################################
		# 評価用のパイソンスクリプトを作成
		es.evaluate_setup(self.py_mod, target_dir, self.f_eval_py)
		
		#####################
		# バッチファイルを作成
		f_batch = os.path.join(target_dir, '_Calc_all.bat')
		with open(f_batch,'w') as f:
			# f.write(batch_all)
			f.write(batch)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)

		return

	#####################################################################################
	# ターミナルのタイトルを設定
	def make_title(self, batch, title):
		if platform.system() == "Windows":
			batch += "title " + title + "\n"
		elif platform.system() == "Linux":
			batch += r'echo -ne "\033]0; ' + title + ' \007"' + '\n'
		return batch

	# ファイル名の処理
	def make_step(self, time, fn_ext, batch):
		present_udf = fn_ext[0] + self.target_name + fn_ext[1]
		out_udf = present_udf.replace("uin", "out")
		batch += self.ver_cognac + ' -I ' + present_udf + ' -O ' + out_udf + self.core + ' \n'
		batch += 'python ' + self.f_eval_py + ' ' + out_udf + '\n'
		read_udf = out_udf
		return present_udf, read_udf, batch


	################################################################################
	def make_base(self, target_dir):
		#--- create an empty UDF file ---
		target_udf = os.path.join(target_dir, self.base_udf)
		with open(target_udf, 'w') as f:
			f.write(r'\include{"%s"}' % self.blank_udf)

		u = UDFManager(target_udf)
		# goto global data
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Solver
		p = 'Simulation_Conditions.Solver.'
		u.put('Dynamics', p + 'Solver_Type')
		u.put('NVT_Kremer_Grest', p + 'Dynamics.Dynamics_Algorithm')
		u.put(0.5, p + 'Dynamics.NVT_Kremer_Grest.Friction')
		# Boundary_Conditions
		p = 'Simulation_Conditions.Boundary_Conditions'
		u.put(['PERIODIC', 'PERIODIC', 'PERIODIC', 1], p)
		#
		p = "Simulation_Conditions.Dynamics_Conditions.Moment."
		u.put(10000, p + "Interval_of_Calc_Moment")
		u.put(1, p + "Calc_Moment")
		u.put(1, p + "Stop_Translation")
		u.put(1, p + "Stop_Rotation")

		#--- Initial_Structure ---
		# Initial_Unit_Cell
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0, p + 'Density')
		u.put([self.system_size, self.system_size, self.system_size, 90.0, 90.0, 90.0], p + 'Cell_Size')

		#--- Molecular_Attributes ---
		# Atomes
		for i, atomname in enumerate(self.atom_name):
			p = 'Molecular_Attributes.Atom_Type[].'
			u.put(atomname, 	p + 'Name', [i])
			u.put(1.0, 			p + 'Mass', [i])
		# Bond
		for i, bondname in enumerate(self.bond_name):
			p = 'Molecular_Attributes.Bond_Potential[].'
			u.put(bondname, 				p + 'Name', [i])
			u.put('Harmonic', 				p + 'Potential_Type', [i])
			u.put(self.harmonic_low[0], 	p + 'R0', [i])
			u.put(self.harmonic_low[1], 	p + 'Harmonic.K', [i])
		# Angle
		for i, anglename in enumerate(self.angle_name):
			p = 'Molecular_Attributes.Angle_Potential[].'
			u.put(anglename, 			p + 'Name', [i])
			u.put(self.angle_low[0], 	p + 'Potential_Type', [i])
			u.put(self.angle_low[1], 	p + 'theta0', [i])
			u.put(self.angle_low[2], 	p + 'Theta2.K', [i])

		# Site
		for i, sitename in enumerate(self.site_name):
			p = 'Molecular_Attributes.Interaction_Site_Type[].'
			u.put(sitename, 		p + 'Name', [i])
			u.put(1, 				p + 'Num_of_Atoms', [i])
			u.put(self.lj_cond[4], 	p + 'Range', [i])

		#--- Pair_Interaction[] ---
		for i, pairname in enumerate(self.pair_name):
			p = 'Interactions.Pair_Interaction[].'
			u.put(pairname,   					p + 'Name', [i])
			u.put('Lennard_Jones', 				p + 'Potential_Type', [i])
			u.put(self.site_pair_name[i][0],	p + 'Site1_Name', [i])
			u.put(self.site_pair_name[i][1],	p + 'Site2_Name', [i])
			u.put(self.lj_cond[0],				p + 'Cutoff', [i])
			u.put(self.lj_cond[1],				p + 'Scale_1_4_Pair', [i])
			u.put(self.lj_cond[2],				p + 'Lennard_Jones.sigma', [i])
			u.put(self.lj_cond[3],				p + 'Lennard_Jones.epsilon', [i])

		#--- Write UDF ---
		u.write(target_udf)
		return

	################################################################################
	# すべてのアトムの位置座標及びボンド情報を設定
	def setup_atoms(self, target_dir):
		# 多重配置の場合の位置シフト量
		shift = 2
		#
		target_udf = os.path.join(target_dir, self.base_udf)
		u = UDFManager(target_udf)
		u.jump(-1)

		#--- Set_of_Molecules の入力
		p = 'Set_of_Molecules.molecule[].'
		pa = p + 'atom[].'
		pi = p + 'interaction_Site[].'
		pb = p + 'bond[].'
		pang = p +  'angle[].'
		#
		count = 0
		for mul in range(len(self.calcd_data_dic_list)):
			atom_all = self.calcd_data_dic_list[mul]["atom_all"]
			bond_all = self.calcd_data_dic_list[mul]["bond_all"]
			pos_all = self.calcd_data_dic_list[mul]["pos_all"]
			angle_all = self.calcd_data_dic_list[mul]["angle_all"]
			#
			u.put(self.nw_name + '_' + str(mul), p + 'Mol_Name', [count])
			# beads
			n_atom = 0
			for atom in atom_all:
				# atom
				id_shift = len(atom_all)
				u.put(atom[0] + count*id_shift, 	pa + 'Atom_ID', [count, n_atom])
				u.put(self.atom_name[atom[1]], 		pa + 'Atom_Name', [count, n_atom])
				u.put(self.atom_name[atom[1]], 		pa + 'Atom_Type_Name', [count, n_atom])
				u.put(0, 							pa + 'Chirality', [count, n_atom])
				u.put(1, 							pa + 'Main_Chain', [count, n_atom])
				# interaction site
				u.put(self.site_name[atom[2]], 		pi + 'Type_Name', [count, n_atom])
				u.put(n_atom, 						pi + 'atom[]', [count, n_atom, 0])
				n_atom += 1
			# bonds
			n_bond = 0
			ang_list=[]
			for bond in range(len(bond_all)):
				u.put(self.bond_name[bond_all[bond][0]], 	pb + 'Potential_Name', [count, n_bond])
				u.put(bond_all[bond][1][0], 				pb + 'atom1', [count, n_bond])
				ang_list.append(bond_all[bond][1][0])
				u.put(bond_all[bond][1][1], 				pb + 'atom2', [count, n_bond])
				n_bond += 1

			# angles
			n_ang = 0
			for a_list in angle_all:
				for j in range(len(a_list)-2):
					u.put(self.angle_name[0], 	pang + 'Potential_Name', [count, n_ang])
					u.put(a_list[j], 			pang + 'atom1', [count, n_ang])
					u.put(a_list[j + 1], 		pang + 'atom2', [count, n_ang])
					u.put(a_list[j + 2], 		pang + 'atom3', [count, n_ang])
					n_ang += 1

			# Draw_Attributes
			color = ["Red", "Green", "Blue", "Magenta", "Cyan", "Yellow", "White", "Black", "Gray"]
			mm = mul % 9
			u.put([self.nw_name + '_' + str(mul), color[mm], 1.0, 1.0], 'Draw_Attributes.Molecule[]', [count])
			count += 1

		# アトムの座標位置を、シフトしながら、設定
		sp = 'Structure.Position.mol[].atom[]'
		count = 0
		self.totalatom = 0
		self.expand = 2.
		#
		for mul in range(len(self.calcd_data_dic_list)):
			shift_vec = self.expand*count*shift*np.array(np.random.rand(3))
			pos_all = self.calcd_data_dic_list[mul]["pos_all"]
			for i in range(len(pos_all)):
				if self.nw_type == "KG_simple":
					mod_pos = self.expand*self.a_cell*np.array(list(pos_all[i])) + shift_vec
				else:
					mod_pos = self.a_cell*np.array(list(pos_all[i])) + shift_vec
				u.put(list(mod_pos), sp, [count, i])
				self.totalatom +=1
			count+=1

		#--- Write UDF ---
		u.write(target_udf)
		return

################################################################################
	# 初期状態を設定
	def setup_nolj_udf(self, template, present_udf, time, harmonic, target_dir):
		u = UDFManager(os.path.join(target_dir, template))
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
		u.put(0, p + 'Angle')
		u.put(0, p + 'Non_Bonding_Interchain')
		u.put(0, p + 'Non_Bonding_1_3')
		u.put(0, p + 'Non_Bonding_1_4')
		u.put(0, p + 'Non_Bonding_Intrachain')

		# Output_Flags.Statistics
		p = 'Simulation_Conditions.Output_Flags.Statistics.'
		u.put(1, p + 'Energy')
		u.put(1, p + 'Temperature')
		u.put(1, p + 'Pressure')
		u.put(0, p + 'Stress')
		u.put(0, p + 'Volume')
		u.put(1, p + 'Density')
		u.put(1, p + 'Cell')
		u.put(0, p + 'Wall_Pressure')
		u.put(0, p + 'Energy_Flow')

		# Output_Flags.Structure
		p = 'Simulation_Conditions.Output_Flags.Structure.'
		u.put(1, p + 'Position')
		u.put(0, p + 'Velocity')
		u.put(0, p + 'Force')

		# #--- Initial_Structure ---
		# # Initial_Unit_Cell
		# p = 'Initial_Structure.Initial_Unit_Cell.'
		# u.put(0.85, p + 'Density')
		# u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')

		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 		p + 'Method')
		u.put(['', -1, 0, 0], 	p + 'Restart')

		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		u.put(1, p + 'Relaxation')
		u.put('DYNAMICS', p + 'Method')
		u.put(300, p + 'Max_Relax_Force')
		u.put(10000, p + 'Max_Relax_Steps')

		#--- Simulation_Conditions ---
		# Bond
		for i, bondname in enumerate(self.bond_name):
			p = 'Molecular_Attributes.Bond_Potential[].'
			u.put(bondname, 		p + 'Name', [i])
			u.put('Harmonic', 		p + 'Potential_Type', [i])
			u.put(harmonic[0], p + 'R0', [i])
			u.put(harmonic[1], p + 'Harmonic.K', [i])

		#--- Write UDF ---
		u.write(os.path.join(target_dir, present_udf))
		return

	################################################################################
	def lj_harm(self, template, read_udf, present_udf, time, harmonic, target_dir):
		u = UDFManager(os.path.join(target_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000000.,	p + 'Max_Force')
		u.put(time[0],		p + 'Time.delta_T')
		u.put(time[1], 		p + 'Time.Total_Steps')
		u.put(time[2], 		p + 'Time.Output_Interval_Steps')
		u.put(1.0, 			p + 'Temperature.Temperature')
		u.put(0., 			p + 'Pressure_Stress.Pressure')

		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(1, p + 'Bond')
		u.put(0, p + 'Angle')
		u.put(1, p + 'Non_Bonding_Interchain')
		u.put(1, p + 'Non_Bonding_1_3')
		u.put(1, p + 'Non_Bonding_1_4')
		u.put(1, p + 'Non_Bonding_Intrachain')

		# Output_Flags.Statistics
		p = 'Simulation_Conditions.Output_Flags.Statistics.'
		u.put(1, p + 'Energy')
		u.put(1, p + 'Temperature')
		u.put(1, p + 'Pressure')
		u.put(1, p + 'Stress')
		u.put(0, p + 'Volume')
		u.put(1, p + 'Density')
		u.put(1, p + 'Cell')
		u.put(0, p + 'Wall_Pressure')
		u.put(0, p + 'Energy_Flow')

		# Output_Flags.Structure
		p = 'Simulation_Conditions.Output_Flags.Structure.'
		u.put(1, p + 'Position')
		u.put(0, p + 'Velocity')
		u.put(0, p + 'Force')

		#--- Initial_Structure ---
		# Initial_Unit_Cell
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0.85, p + 'Density')
		u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')

		if read_udf == '':
			# Generate_Method
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', 		p + 'Method')
			u.put(['', -1, 0, 0], 	p + 'Restart')

			# Relaxation
			p = 'Initial_Structure.Relaxation.'
			u.put(1, p + 'Relaxation')
			u.put('DYNAMICS', p + 'Method')
			u.put(300, p + 'Max_Relax_Force')
			u.put(10000, p + 'Max_Relax_Steps')
		else:
			# Generate_Method
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', 			p + 'Method')
			u.put([read_udf, -1, 1, 0], p + 'Restart')

		#--- Simulation_Conditions ---
		# Bond
		for i, bondname in enumerate(self.bond_name):
			p = 'Molecular_Attributes.Bond_Potential[].'
			u.put(bondname, 		p + 'Name', [i])
			u.put('Harmonic', 		p + 'Potential_Type', [i])
			u.put(harmonic[0], p + 'R0', [i])
			u.put(harmonic[1], p + 'Harmonic.K', [i])

		#--- Write UDF ---
		u.write(os.path.join(target_dir, present_udf))
		return

	################################################################################
	def ang_no_lj(self, template, read_udf, present_udf, time, angle, target_dir):
		u = UDFManager(os.path.join(target_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time[0], p + 'Time.delta_T')
		u.put(time[1], p + 'Time.Total_Steps')
		u.put(time[2], p + 'Time.Output_Interval_Steps')

		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(1, p + 'Bond')
		u.put(1, p + 'Angle')
		u.put(0, p + 'Non_Bonding_Interchain')
		u.put(0, p + 'Non_Bonding_1_3')
		u.put(0, p + 'Non_Bonding_1_4')
		u.put(0, p + 'Non_Bonding_Intrachain')

		#--- Initial_Structure ---
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 			p + 'Method')
		u.put([read_udf, -1, 1, 0], p + 'Restart')

		# Angle
		for i, anglename in enumerate(self.angle_name):
			p = 'Molecular_Attributes.Angle_Potential[].'
			u.put(anglename, 			p + 'Name', [i])
			u.put(angle[0], 	p + 'Potential_Type', [i])
			u.put(angle[1], 	p + 'theta0', [i])
			u.put(angle[2], 	p + 'Theta2.K', [i])

		#--- Write UDF ---
		u.write(os.path.join(target_dir, present_udf))
		return

	###########################
	# 斥力的な非結合相互作用を設定
	def ang_w_lj(self, template, read_udf, present_udf, time, angle, target_dir):
		u = UDFManager(os.path.join(target_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time[0], p + 'Time.delta_T')
		u.put(time[1], p + 'Time.Total_Steps')
		u.put(time[2], p + 'Time.Output_Interval_Steps')

		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(1, p + 'Bond')
		u.put(1, p + 'Angle')
		u.put(1, p + 'Non_Bonding_Interchain')
		u.put(1, p + 'Non_Bonding_1_3')
		u.put(1, p + 'Non_Bonding_1_4')
		u.put(1, p + 'Non_Bonding_Intrachain')

		#--- Initial_Structure ---
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 		p + 'Method')
		u.put([read_udf, -1, 1, 0], 	p + 'Restart')

		# Angle
		for i, anglename in enumerate(self.angle_name):
			p = 'Molecular_Attributes.Angle_Potential[].'
			u.put(anglename, 			p + 'Name', [i])
			u.put(angle[0], 	p + 'Potential_Type', [i])
			u.put(angle[1], 	p + 'theta0', [i])
			u.put(angle[2], 	p + 'Theta2.K', [i])

		#--- Write UDF ---
		u.write(os.path.join(target_dir, present_udf))

		return

	################################################################################
	def bond_fene(self, template, read_udf, present_udf, rmax, time, target_dir):
		u = UDFManager(os.path.join(target_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time[0], p + 'Time.delta_T')
		u.put(time[1], p + 'Time.Total_Steps')
		u.put(time[2], p + 'Time.Output_Interval_Steps')
		# Solver
		p = 'Simulation_Conditions.Solver.'
		u.put('Dynamics', p + 'Solver_Type')
		u.put('NVT_Kremer_Grest', p + 'Dynamics.Dynamics_Algorithm')
		u.put(0.5, p + 'Dynamics.NVT_Kremer_Grest.Friction')
		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(0, p + 'Angle')

		#--- Initial_Structure ---
		# Initial_Unit_Cell
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0.85, p + 'Density')
		u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', p + 'Method')
		u.put([read_udf, -1, 0, 0], p + 'Restart')
		#--- Simulation_Conditions ---
		# Bond
		# bond_name = u.get('Molecular_Attributes.Bond_Potential[].Name')
		for i, b_name in enumerate(self.bond_name):
			p = 'Molecular_Attributes.Bond_Potential[].'
			u.put(b_name, 		p + 'Name', [i])
			u.put('FENE_LJ', 	p + 'Potential_Type', [i])
			u.put(1.0,			p + 'R0', [i])
			u.put(rmax,			p + 'FENE_LJ.R_max', [i])
			u.put(30,			p + 'FENE_LJ.K', [i])
			u.put(1.0,			p + 'FENE_LJ.sigma', [i])
			u.put(1.0,			p + 'FENE_LJ.epsilon', [i])

		#--- Write UDF ---
		u.write(os.path.join(target_dir, present_udf))
		return
	
	################################################################################
	def kg_npt(self, template, read_udf, present_udf, time, pressure, target_dir):
		u = UDFManager(os.path.join(target_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000000.,	p + 'Max_Force')
		u.put(time[0],		p + 'Time.delta_T')
		u.put(time[1], 		p + 'Time.Total_Steps')
		u.put(time[2], 		p + 'Time.Output_Interval_Steps')
		u.put(1.0, 			p + 'Temperature.Temperature')
		u.put(0., 			p + 'Pressure_Stress.Pressure')

		# Solver
		p = 'Simulation_Conditions.Solver.'
		u.put('Dynamics', p + 'Solver_Type')
		u.put('NPT_Andersen_Kremer_Grest', p + 'Dynamics.Dynamics_Algorithm')
		u.put(self.totalatom, p + 'Dynamics.NPT_Andersen_Kremer_Grest.Cell_Mass')
		u.put(0.5, p + 'Dynamics.NPT_Andersen_Kremer_Grest.Friction')

		# Pressure
		u.put(pressure, 'Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure')

		# Calc_Potential_Flags
		p = 'Simulation_Conditions.Calc_Potential_Flags.'
		u.put(1, p + 'Bond')
		u.put(0, p + 'Angle')
		u.put(1, p + 'Non_Bonding_Interchain')
		u.put(1, p + 'Non_Bonding_1_3')
		u.put(1, p + 'Non_Bonding_1_4')
		u.put(1, p + 'Non_Bonding_Intrachain')

		# Output_Flags.Statistics
		p = 'Simulation_Conditions.Output_Flags.Statistics.'
		u.put(1, p + 'Energy')
		u.put(1, p + 'Temperature')
		u.put(1, p + 'Pressure')
		u.put(1, p + 'Stress')
		u.put(0, p + 'Volume')
		u.put(1, p + 'Density')
		u.put(1, p + 'Cell')
		u.put(0, p + 'Wall_Pressure')
		u.put(0, p + 'Energy_Flow')

		# Output_Flags.Structure
		p = 'Simulation_Conditions.Output_Flags.Structure.'
		u.put(1, p + 'Position')
		u.put(0, p + 'Velocity')
		u.put(0, p + 'Force')

		#--- Initial_Structure ---
		# Initial_Unit_Cell
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0.85/self.expand**3, p + 'Density')
		u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')

		if read_udf == '':
			# Generate_Method
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', 		p + 'Method')
			u.put(['', -1, 0, 0], 	p + 'Restart')

			# Relaxation
			p = 'Initial_Structure.Relaxation.'
			u.put(1, p + 'Relaxation')
			u.put('DYNAMICS', p + 'Method')
			u.put(300, p + 'Max_Relax_Force')
			u.put(10000, p + 'Max_Relax_Steps')
		else:
			# Generate_Method
			p = 'Initial_Structure.Generate_Method.'
			u.put('Restart', 			p + 'Method')
			u.put([read_udf, -1, 1, 0], p + 'Restart')

		#--- Simulation_Conditions ---
		# Bond
		for i, bondname in enumerate(self.bond_name):
			p = 'Molecular_Attributes.Bond_Potential[].'
			u.put(bondname, 		p + 'Name', [i])
			u.put('Harmonic', 		p + 'Potential_Type', [i])
			u.put(self.harmonic_high[0], p + 'R0', [i])
			u.put(self.harmonic_high[1], p + 'Harmonic.K', [i])
		# Angle
		for i, anglename in enumerate(self.angle_name):
			p = 'Molecular_Attributes.Angle_Potential[].'
			u.put(anglename, 			p + 'Name', [i])
			u.put(self.angle_high[0], 	p + 'Potential_Type', [i])
			u.put(self.angle_high[1], 	p + 'theta0', [i])
			u.put(self.angle_high[2], 	p + 'Theta2.K', [i])

		#--- Write UDF ---
		u.write(os.path.join(target_dir, present_udf))
		return

	################################################################################
	def eq_udf(self, template, read_udf, present_udf, time, target_dir):
		u = UDFManager(os.path.join(target_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time[0],  p+'Time.delta_T')
		u.put(time[1],  p+'Time.Total_Steps')
		u.put(time[2],  p+'Time.Output_Interval_Steps')

		# Moment
		p = "Simulation_Conditions.Dynamics_Conditions.Moment."
		u.put(0, p + "Interval_of_Calc_Moment")
		u.put(0, p + "Calc_Moment")
		u.put(0, p + "Stop_Translation")
		u.put(0, p + "Stop_Rotation")

		#--- Initial_Structure ---
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', p+'Method')
		u.put([read_udf, -1, 1, 0], p+'Restart')
		p = 'Initial_Structure.Relaxation.'
		u.put(0, p + 'Relaxation')

		#--- Write UDF ---
		u.write(os.path.join(target_dir, present_udf))
		return


################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()
