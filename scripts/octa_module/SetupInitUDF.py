#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
import sys
import numpy as np

# import platform
from UDFManager import UDFManager
import os

##########################################
# Initial_UDF の作成
##########################################
class MakeInitUDF:
	def __init__(self, sim_cond, target_cond, files_cond, names, target_name, calcd_data_dic, expand):
		#
		self.nw_model = sim_cond[0]
		self.nw_type = sim_cond[1]
		self.n_segments = sim_cond[2]
		self.n_cell = sim_cond[3]
		self.target_density = sim_cond[5]
		self.n_strand = sim_cond[6]
		self.l_bond = sim_cond[7]
		self.c_n = sim_cond[8]
		#
		self.system_size = target_cond[0]
		self.a_cell = target_cond[1]
		self.total_atom = target_cond[2]
		self.nu = target_cond[3]
		self.structure = target_cond[4]
		self.multi = target_cond[5]
		self.n_solvent = target_cond[6]
		#
		self.ver_cognac = files_cond[0]
		self.blank_udf = files_cond[1]
		self.base_udf = files_cond[2]
		self.core = ' -n ' + str(files_cond[3])
		# Cognac用の名称設定
		self.nw_name = names[0]
		self.atom_name = names[1]
		self.bond_name = names[2]
		self.angle_name = names[3]
		self.site_name = names[4]
		self.pair_name = names[5]
		self.site_pair_name = names[6]
		#
		self.target_name = target_name
		self.calcd_data_dic = calcd_data_dic
		self.expand = expand
		# 条件設定
		# [R0, K]
		self.harmonic = [0.97, 1000]
		# [Potential_Type, theta0, K]		
		self.angle = ['Theta2', 74, 10.0]
		# [Cutoff, Scale_1_4_Pair, sigma, epsilon, range]
		self.lj_cond = [2**(1/6), 1.0, 1.0, 1.0, 1.0]	
	################################################################################
	# UDFファイルを設定し、バッチ処理を作成
	def setup_baseudf(self):
		# 計算用のディレクトリーを作成
		self.make_dir()
		# base_udfを作成
		self.make_base_udf()

		return self.target_dir
	################################################################################
	# 計算用のディレクトリーを作成
	def make_dir(self):
		self.target_dir = str(self.target_name)
		os.makedirs(self.target_dir, exist_ok = True)
		with open(os.path.join(self.target_dir, "calc.dat"), "w") as f:
			f.write("# segments\tbond_length\tCN\tfunc\tnu\tNW_type\n" + str(self.n_segments) + '\t' + str(self.l_bond) + '\t' + str(self.c_n) + '\t' + str(self.n_strand) + "\t" + str(round(self.nu, 5)) + '\t' + self.structure)
		return

	############################################
	# base_udfを作成
	def make_base_udf(self):
		# 初期udfの内容を作成する
		self.make_base()
		# すべてのアトムの位置座標及びボンド情報を設定
		self.setup_atoms()
		return

	################################################################################
	def make_base(self):
		#--- create an empty UDF file ---
		target_udf = os.path.join(self.target_dir, self.base_udf)
		with open(target_udf, 'w') as f:
			f.write(r'\include{"%s"}' % self.blank_udf)

		u = UDFManager(target_udf)
		# goto global data
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Solver
		p = 'Simulation_Conditions.Solver.'
		u.put('Dynamics', p + 'Solver_Type')
		if self.nw_type == "KG_NPT":
			u.put('NPT_Andersen_Kremer_Grest', 	p + 'Dynamics.Dynamics_Algorithm')
			u.put(self.total_atom, 				p + 'Dynamics.NPT_Andersen_Kremer_Grest.Cell_Mass')
			u.put(0.5, 							p + 'Dynamics.NPT_Andersen_Kremer_Grest.Friction')
		else:
			u.put('NVT_Kremer_Grest', 			p + 'Dynamics.Dynamics_Algorithm')
			u.put(0.5, 							p + 'Dynamics.NVT_Kremer_Grest.Friction')
		# Boundary_Conditions
		p = 'Simulation_Conditions.Boundary_Conditions'
		u.put(['PERIODIC', 'PERIODIC', 'PERIODIC', 1], p)
		#
		p = "Simulation_Conditions.Dynamics_Conditions.Moment."
		u.put(10000, p + "Interval_of_Calc_Moment")
		u.put(1, p + "Calc_Moment")
		u.put(1, p + "Stop_Translation")
		u.put(1, p + "Stop_Rotation")

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
		u.put(1, p + 'Volume')
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
		if self.nw_type == "KG_NPT":
			p = 'Initial_Structure.Initial_Unit_Cell.'
			u.put(0, p + 'Density')
			u.put([self.system_size*self.expand, self.system_size*self.expand, self.system_size*self.expand, 90.0, 90.0, 90.0], p + 'Cell_Size')
		else:
			u.put(self.target_density, p + 'Density')
			u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
		
		#--- Molecular_Attributes ---
		# Atomes
		for i, atomname in enumerate(self.atom_name):
			p = 'Molecular_Attributes.Atom_Type[].'
			u.put(atomname, 	p + 'Name', [i])
			u.put(1.0, 			p + 'Mass', [i])
		# Bond
		for i, bondname in enumerate(self.bond_name):
			p = 'Molecular_Attributes.Bond_Potential[].'
			u.put(bondname, 			p + 'Name', [i])
			u.put('Harmonic', 			p + 'Potential_Type', [i])
			u.put(self.harmonic[0], 	p + 'R0', [i])
			u.put(self.harmonic[1], 	p + 'Harmonic.K', [i])
		# Angle
		for i, anglename in enumerate(self.angle_name):
			p = 'Molecular_Attributes.Angle_Potential[].'
			u.put(anglename, 			p + 'Name', [i])
			u.put(self.angle[0], 	p + 'Potential_Type', [i])
			u.put(self.angle[1], 	p + 'theta0', [i])
			u.put(self.angle[2], 	p + 'Theta2.K', [i])

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
	def setup_atoms(self):
		# 多重配置の場合の位置シフト量
		shift = 2
		#
		target_udf = os.path.join(self.target_dir, self.base_udf)
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
		for mul in range(len(self.calcd_data_dic)):
			atom_all = self.calcd_data_dic[mul]["atom_all"]
			bond_all = self.calcd_data_dic[mul]["bond_all"]
			pos_all = self.calcd_data_dic[mul]["pos_all"]
			angle_all = self.calcd_data_dic[mul]["angle_all"]
			#
			u.put(self.nw_name + '_' + str(mul), p + 'Mol_Name', [count])
			# beads
			n_atom = 0
			for atom in atom_all:
				# atom
				id_shift = len(atom_all) + self.n_solvent
				atom_id = atom[0] + count*id_shift
				u.put(atom_id, 						pa + 'Atom_ID', [count, n_atom])
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
		# tmp_set = []
		# ネットワークのセグメントをセット
		for mul in range(len(self.calcd_data_dic)):
			shift_vec = self.expand*count*shift*np.array(np.random.rand(3))
			pos_all = self.calcd_data_dic[mul]["pos_all"]
			for i in range(len(pos_all)):
				if 'NPT' in self.target_name or 'single' in self.target_name:
					mod_pos = self.expand*self.a_cell*np.array(list(pos_all[i])) + self.expand*shift_vec
				else:
					mod_pos = self.a_cell*np.array(list(pos_all[i])) + shift_vec
				# print(mod_pos)
				u.put(list(mod_pos), sp, [count, i])
				# tmp_set.append(list(mod_pos))
				self.totalatom +=1
			count+=1
	# ソルベントのアトムをセット
		if self.n_solvent != 0:
			sol_pos = [n*self.system_size for n in np.random.random_sample((self.n_solvent,3))]
			for mol in range(self.n_solvent):
				u.put('Solvent', p + 'Mol_Name', [count])
				atom_id += 1
				# atom
				u.put(atom_id,				pa + 'Atom_ID', [count, 0])
				u.put("Solvent", 			pa + 'Atom_Name', [count, 0])
				u.put("Solvent", 			pa + 'Atom_Type_Name', [count, 0])
				u.put(0, 					pa + 'Chirality', [count, 0])
				u.put(1, 					pa + 'Main_Chain', [count, 0])
				# interaction site
				u.put("site_Solvent", 		pi + 'Type_Name', [count, 0])
				u.put(0, 					pi + 'atom[]', [count, 0, 0])
				# position
				u.put(list(sol_pos[mol]), sp, [count, 0])
				count += 1

		#--- Write UDF ---
		u.write(target_udf)

		return
