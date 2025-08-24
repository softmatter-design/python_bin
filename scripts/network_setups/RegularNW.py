#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
import sys
import numpy as np

# import platform
from UDFManager import UDFManager
import os
######################################
##### ネットワークの初期設定を行う #####
######################################
class InitialSetup:
	def __init__(self, sim_cond):
		#
		self.nw_model = sim_cond[0]
		self.nw_type = sim_cond[1]
		self.n_segments = sim_cond[2]
		self.n_cell = sim_cond[3]
		self.multi_init = sim_cond[4]
		self.target_density = sim_cond[5]
		self.n_strand = sim_cond[6]
		self.l_bond = sim_cond[7]
		self.c_n = sim_cond[8]

	##########################################
	##### ネットワークポリマーの諸量を計算 ######
	##########################################
	#-----ネットワークポリマーの諸量を計算
	def calc_conditions(self):
		## 計算システムの諸量を計算して、出力
		target_cond = self.init_calc()
		target_name = self.set_target_name()
		# ネットワークの計算
		calcd_data_dic = self.calc_all()	
		return target_name, target_cond, target_name, calcd_data_dic

	################################################################################
	def init_calc(self):
		structure = "Regular_NW"
		e2e = self.l_bond*((self.n_segments + 1)*self.c_n)**0.5	# 理想鎖状態での末端間距離
		#
		if self.nw_model == "3_Chain_S":
			n_chains = 12						        # サブチェインの本数
			n_p_unit = 8 + self.n_segments*n_chains		# ユニットセル当たりの粒子数
			a_cell = (2*2**0.5)*e2e				        # 理想鎖状態でのユニットセル長
		elif self.nw_model == "3_Chain_D":
			n_chains = 24						        # サブチェインの本数
			n_p_unit = 16 + self.n_segments*n_chains	# ユニットセル当たりの粒子数
			a_cell = (2*2**0.5)*e2e				        # 理想鎖状態でのユニットセル長
		elif self.nw_model == "4_Chain":
			n_chains = 16						        # サブチェインの本数
			n_p_unit = 8 + self.n_segments*n_chains		# ユニットセル当たりの粒子数
			a_cell = (4*3**0.5)*e2e/3			        # 理想鎖状態でのユニットセル長
		elif self.nw_model == "6_Chain":
			n_chains = 3						        # サブチェインの本数
			n_p_unit = 1 + self.n_segments*n_chains		# ユニットセル当たりの粒子数
			a_cell = e2e						            # 理想鎖状態でのユニットセル長
		elif self.nw_model == "8_Chain":
			n_chains = 8						        # サブチェインの本数
			n_p_unit = 2 + self.n_segments*n_chains     # ユニットセル当たりの粒子数
			a_cell = (2*3**0.5)*e2e/3					# 理想鎖状態でのユニットセル長
		#
		n_solvent = 0
		#
		if self.nw_type == "KG_entangled" or self.nw_type == "KG_NPT":
			multi_calcd = round(self.target_density*a_cell**3/n_p_unit)		# 密度を設定値とした場合に必要な多重度
			self.multi = multi_calcd
			fin_dens = n_p_unit*self.multi/a_cell**3						# 上記の多重度での密度
			err_dens = round((fin_dens/self.target_density - 1)*100, 2) 	# 設定密度との誤差
			single_net_atom = int(n_p_unit*self.n_cell**3.)	    			# 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    			# 全システム中のネットワーク粒子数
			system = (total_net_atom/self.target_density)**(1/3)			# その際のシステムサイズ
			nu = n_chains*self.multi/a_cell**3.								# ストランドの数密度
		elif self.nw_type == "KG_single" or self.nw_type == "Sunuke_dens":
			self.multi = multi_init
			system = a_cell*self.n_cell							# e2e から決めたシステムサイズ
			single_net_atom = int(n_p_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    # 全システム中のネットワーク粒子数
			vol = total_net_atom/self.target_density			# システム体積
			mod_system = vol**(1./3.)							# 収縮後のシステムサイズ
			shrinkage = mod_system/system						# 収縮比
			mod_e2e = shrinkage*e2e								# 収縮後の末端間距離
			nu = n_chains*self.multi*self.n_cell**3./vol		# ストランドの数密度
		elif self.nw_type == "KG_gel":
			self.multi = multi_init
			system = a_cell*self.n_cell						    # システムサイズ
			single_net_atom = int(n_p_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    # 全システム中のネットワーク粒子数
			vol = system**3									    # システム体積
			density = total_net_atom/vol                        # 密度
			#
			total_atom = int(vol*self.target_density)		    # システム中の粒子総数
			n_solvent = int(total_atom - total_net_atom)		# 溶媒の粒子数
			nv = total_net_atom/total_atom						# ネットワークの体積分率
			nu = n_chains*self.multi*self.n_cell**3./vol		# ストランドの数密度
		elif self.nw_type == "Sunuke":
			self.multi = multi_init
			system = a_cell*self.n_cell						    # システムサイズ
			single_net_atom = int(n_p_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(self.multi*single_net_atom)    # 全システム中のネットワーク粒子数
			vol = system**3									    # システム体積
			density = total_net_atom/vol                        # 密度
			nu = n_chains*self.multi*self.n_cell**3./vol		# ストランドの数密度
		else:
			sys.exit("Something Wrong!!")
		#
		text = "#########################################" + "\n"
		text += "ネットワークモデル\t\t" + str(self.nw_model) + "\n"
		text += "ネットワークタイプ\t\t" + str(self.nw_type) + "\n"
		text += "ストランド中のセグメント数\t" + str(self.n_segments) + "\n"
		text += "特性比：\t\t\t" + str(round(self.c_n,2)) + "\n"
		text += "末端間距離：\t\t\t" + str(round(e2e,5)) + "\n"
		text += "一辺当たりの単位ユニット数\t" + str(self.n_cell) + "\n"
		text += "#########################################" + "\n"
		if self.nw_type == "KG_entangled" or self.nw_type == "KG_NPT":
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "セグメント数:\t\t\t" + str(total_net_atom) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "密度:\t\t\t\t" + str(round(fin_dens, 4)) + "\n"
			text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			if abs(err_dens) > 1:
				print(u"\n##### \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？")
				self.prompt()
			else:
				self.prompt()
		elif self.nw_type == "KG_single" or self.nw_type == "Sunuke_dens":
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "セグメント数:\t\t\t" + str(total_net_atom) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "収縮比:\t\t\t\t" + str(round(shrinkage, 4)) + "\n"
			text += "収縮後の末端間距離：\t\t" + str(round(mod_e2e,5)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			self.prompt()
		elif self.nw_type == "Sunuke":
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "セグメント数:\t\t\t" + str(total_net_atom) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "密度:\t\t\t\t" + str(round(density, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			self.prompt()
		elif self.nw_type == "KG_gel":
			text += "設定密度:\t\t\t" + str(self.target_density) + "\n"
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "単一ネットワークのセグメント数:\t" + str(single_net_atom) + "\n"
			text += "溶媒のセグメント数:\t\t" + str(n_solvent) + "\n"
			text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
			text += "ネットワークの体積分率:\t\t" + str(round(nv, 4)) + "\n"
			text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
			text += "#########################################" + "\n"
			print(text)
			self.prompt()
		#
		with open("calc_conditions.txt", 'w') as f:
			f.write(text)
		#
		target_cond = [system, a_cell, total_net_atom, nu, structure, self.multi, n_solvent]
		return target_cond

	###############
	# 計算条件の確認
	def prompt(self):
		dic={'y':True,'yes':True,'n':False,'no':False}
		while True:
			inp = input(u'計算を続行 ==> [Y]es >> ').lower()
			if inp in dic:
				inp = dic[inp]
				break
			print(u'##### \nもう一度入力してください。')
		if inp :
			print(u"計算を続行します")
		else:
			sys.exit(u"##### \n計算を中止します。")
		return
	
	###################
	# target_nameを決める。
	def set_target_name(self):
		target_name = str(self.nw_type ) + "_" + str(self.nw_model) + '_N_' + str(self.n_segments) + "_Cells_" + str(self.n_cell) + "_Multi_" + str(self.multi)
		return target_name

	################################################################################
	# ネットワーク設定の計算
	################################################################################
	def calc_all(self):
		# 架橋点 JP を設定
		jp_xyz, subchain_se_xyz = self.calc_jp_subChains()
		#
		calcd_data_dic = self.set_atom(jp_xyz, subchain_se_xyz)
		return calcd_data_dic
	################################################################################
	# JPおよびサブチェインの始点と終点のXYZを設定
	def calc_jp_subChains(self):
		# jp_xyz は、JPの座標のリスト
		# subchain_se_xyz は、サブチェインの出発点と終点のリスト
		if self.nw_model == "3_Chain_S":
			# JPを設定
			jp_xyz = [
			[
			[0, 0, 0],
			[0, 0.25, 0.25],
			[0.25, 0.25, 0.5],
			[0.25, 0, 0.75],
			[0.5, 0.5, 0.5],
			[0.5, 0.75, 0.75],
			[0.75, 0.5, 0.25],
			[0.75, 0.75, 0]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0, 0, 0], [0, 0.25, 0.25]],
			[[0, 0.25, 0.25], [0.25, 0.25, 0.5]],
			[[0.25, 0.25, 0.5], [0.25, 0, 0.75]],
			[[0.25, 0.25, 0.5], [0.5, 0.5, 0.5]],
			[[0.5, 0.5, 0.5], [0.5, 0.75, 0.75]],
			[[0.5, 0.5, 0.5], [0.75, 0.5, 0.25]],
			[[0.75, 0.5, 0.25], [0.75, 0.75, 0]],
			[[0.75, 0.5, 0.25], [1, 0.25, 0.25]],
			[[0.25, 0, 0.75], [0, 0, 1]],
			[[0.5, 0.75, 0.75], [0.25, 1, 0.75]],
			[[0.5, 0.75, 0.75], [0.75, 0.75, 1]],
			[[0.75, 0.75, 0], [1, 1, 0]]
			]
			]

		elif self.nw_model == "3_Chain_D":
			# JPを設定
			jp_xyz = [
			[
			[0, 0, 0],
			[0, 0.25, 0.25],
			[0.25, 0.25, 0.5],
			[0.25, 0, 0.75],
			[0.5, 0.5, 0.5],
			[0.5, 0.75, 0.75],
			[0.75, 0.5, 0.25],
			[0.75, 0.75, 0]
			],
			[	#ここから二つ目
			[0, 0.5, 0.75],
			[0, 0.75, 0.5],
			[0.25, 0.75, 0.25],
			[0.25, 0.5, 0],
			[0.5, 0.25, 0],
			[0.5, 0, 0.25],
			[0.75, 0, 0.5],
			[0.75, 0.25, 0.75]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0, 0, 0], [0, 0.25, 0.25]],
			[[0, 0.25, 0.25], [0.25, 0.25, 0.5]],
			[[0.25, 0.25, 0.5], [0.25, 0, 0.75]],
			[[0.25, 0.25, 0.5], [0.5, 0.5, 0.5]],
			[[0.5, 0.5, 0.5], [0.5, 0.75, 0.75]],
			[[0.5, 0.5, 0.5], [0.75, 0.5, 0.25]],
			[[0.75, 0.5, 0.25], [0.75, 0.75, 0]],
			[[0.75, 0.5, 0.25], [1, 0.25, 0.25]],
			[[0.25, 0, 0.75], [0, 0, 1]],
			[[0.5, 0.75, 0.75], [0.25, 1, 0.75]],
			[[0.5, 0.75, 0.75], [0.75, 0.75, 1]],
			[[0.75, 0.75, 0], [1, 1, 0]]
			],
			[
			[[0, 0.5, 0.75], [0, 0.75, 0.5]],
			[[0, 0.75, 0.5], [0.25, 0.75, 0.25]],
			[[0.25, 0.75, 0.25], [0.25, 0.5, 0]],
			[[0.25, 0.5, 0], [0.5, 0.25, 0]],
			[[0.5, 0.25, 0], [0.5, 0, 0.25]],
			[[0.5, 0, 0.25], [0.75, 0, 0.5]],
			[[0.75, 0, 0.5], [0.75, 0.25, 0.75]],
			[[0, 0.5, 0.75], [0.25, 0.5, 1]],
			[[0.25, 0.75, 0.25], [0.5, 1, 0.25]],
			[[0.75, 0.25, 0.75], [0.5, 0.25, 1]],
			[[0.75, 0.25, 0.75], [1, 0.5, 0.75]],
			[[0.75, 1, 0.5], [1, 0.75, 0.5]]
			]
			]

		elif self.nw_model == "4_Chain":
			# JPを設定
			jp_xyz = [
			[
			[0.0, 0.0, 0.0],
			[0.0, 0.5, 0.5],
			[0.5, 0.0, 0.5],
			[0.5, 0.5, 0.0],
			[0.25, 0.25, 0.25],
			[0.25, 0.75, 0.75],
			[0.75, 0.25, 0.75],
			[0.75, 0.75, 0.25]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0.25, 0.25, 0.25], [0.0, 0.0, 0.0]],	# No.1
			[[0.25, 0.25, 0.25], [0.0, 0.5, 0.5]],
			[[0.25, 0.25, 0.25], [0.5, 0.0, 0.5]],
			[[0.25, 0.25, 0.25], [0.5, 0.5, 0.0]],
			[[0.25, 0.75, 0.75], [0.0, 0.5, 0.5]],	# No.2
			[[0.25, 0.75, 0.75], [0.0, 1.0, 1.0]],
			[[0.25, 0.75, 0.75], [0.5, 0.5, 1.0]],
			[[0.25, 0.75, 0.75], [0.5, 1.0, 0.5]],
			[[0.75, 0.25, 0.75], [0.5, 0.0, 0.5]],	# No.3
			[[0.75, 0.25, 0.75], [0.5, 0.5, 1.0]],
			[[0.75, 0.25, 0.75], [1.0, 0.0, 1.0]],
			[[0.75, 0.25, 0.75], [1.0, 0.5, 0.5]],
			[[0.75, 0.75, 0.25], [0.5, 0.5, 0.0]],	# No.4
			[[0.75, 0.75, 0.25], [0.5, 1.0, 0.5]],
			[[0.75, 0.75, 0.25], [1.0, 0.5, 0.5]],
			[[0.75, 0.75, 0.25], [1.0, 1.0, 0.0]]
			]
			]

		elif self.nw_model == "6_Chain":
			# JPを設定
			jp_xyz = [
			[
			[0.,0.,0.]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0., 0., 0.], [1, 0, 0]],
			[[0., 0., 0.], [0, 1, 0]],
			[[0., 0., 0.], [0, 0, 1]]
			]
			]

		elif self.nw_model == "8_Chain":
			# JPを設定
			jp_xyz = [
			[
			[0.,0.,0.],
			[0.5,0.5,0.5]
			]
			]
			# サブチェインの出発点と終点を設定
			subchain_se_xyz = [
			[
			[[0.5, 0.5, 0.5], [0, 0, 0]],
			[[0.5, 0.5, 0.5], [1, 0, 0]],
			[[0.5, 0.5, 0.5], [0, 1, 0]],
			[[0.5, 0.5, 0.5], [1, 1, 0]],
			[[0.5, 0.5, 0.5], [0, 0, 1]],
			[[0.5, 0.5, 0.5], [1, 0, 1]],
			[[0.5, 0.5, 0.5], [0, 1, 1]],
			[[0.5, 0.5, 0.5], [1, 1, 1]]
			]
			]

		return jp_xyz, subchain_se_xyz

	#########################################################
	def set_atom(self, jp_xyz, subchain_se_xyz):
		calcd_data_dic={}
		count = 0
		for i in (range(self.multi)):
			for mol, jp in enumerate(jp_xyz):
				atom_all = []
				pos_all = {}
				# システム全体にわたるジャンクションポイントのxyzとIDの辞書を作成
				jp_id_dic, jp_xyz_dic, atom_jp = self.set_jp_id(jp, mol)
				atom_all.extend(atom_jp)
				pos_all.update(jp_xyz_dic)
				# print(jp_xyz_dic)
				# サブチェイン中の各アトムのxyzリストとボンドリストを作成
				strand_xyz, bond_all, atom_sc, angle_all = self.set_subchains(jp_id_dic, subchain_se_xyz[mol], mol)
				#
				atom_all.extend(atom_sc)
				pos_all.update(strand_xyz)
				#
				calcd_data_dic[count] = {"atom_all":atom_all, "bond_all":bond_all, "pos_all":pos_all, "angle_all":angle_all}
				count += 1
		return calcd_data_dic

	###################################################
	# システム全体にわたるJPのxyzとIDの辞書を作成
	def set_jp_id(self, jp_xyz, mol):
		jp_id_dic = {}
		jp_xyz_dic = {}
		atom_jp = []
		jp_id = 0
		for z in range(self.n_cell):
			for y in range(self.n_cell):
				for x in range(self.n_cell):
					base_xyz = np.array([x,y,z])
					for jp in jp_xyz:
						jp_id_dic[tuple(np.array(jp) + base_xyz)] = (jp_id)
						jp_xyz_dic[(jp_id)] = tuple(np.array(jp) + base_xyz)
						atom_jp.append([jp_id, 2*mol + 0, 0])
						jp_id += 1
		return jp_id_dic, jp_xyz_dic, atom_jp

		
	#########################################################
	# サブチェイン中の各アトムのxyzリストとボンドリストを作成
	def set_subchains(self, jp_id_dic, subchain_se_xyz, mol):
		strand_xyz = {}
		bond_all = {}
		atom_sc = []
		angle_all = []
		sub_id = len(jp_id_dic)
		bond_id = 0
		for z in range(self.n_cell):
			for y in range(self.n_cell):
				for x in range(self.n_cell):
					b_xyz = (x,y,z)
					for se_xyz in subchain_se_xyz:
						tmp_xyz, tmp_bond, new_sub_id, new_bond_id, tmp_atom_sc, tmp_angle = self.calc_single_subchain(jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, mol)
						strand_xyz.update(tmp_xyz)
						bond_all.update(tmp_bond)
						atom_sc.extend(tmp_atom_sc)
						angle_all.append(tmp_angle)
						sub_id = new_sub_id
						bond_id = new_bond_id
		return strand_xyz, bond_all, atom_sc, angle_all

	###############################################################
	# 一本のサブチェイン中の各アトムのxyzリストとボンドリストを作成
	def calc_single_subchain(self, jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, mol):
		tmp_xyz = {}
		tmp_bond = {}
		tmp_angle = []
		tmp_atom_sc = []
		bas_xyz = np.array(b_xyz)
		# サブチェインの末端間のベクトルを設定
		start_xyz = np.array(se_xyz[0]) + bas_xyz
		end_xyz = np.array(se_xyz[1]) + bas_xyz
		vec = end_xyz - start_xyz
		# 始点のアトムのIDを設定
		mod_xyz = list(start_xyz)[:]
		for dim in range(3):
			if mod_xyz[dim] == self.n_cell:
				mod_xyz[dim] = 0
		s_id = jp_id_dic[tuple(mod_xyz)]
		tmp_angle.append(s_id)
		# 終点のアトムのIDを周期境界条件で変更
		mod_xyz = list(end_xyz)[:]
		for dim in range(3):
			if mod_xyz[dim] == self.n_cell:
				mod_xyz[dim] = 0
		E_id = jp_id_dic[tuple(mod_xyz)]
		# サブチェインの鎖長分のループ処理
		for seg in range(self.n_segments):
			tmp_xyz[sub_id] = tuple(start_xyz + vec*(seg+1)/(self.n_segments+1.))
			if seg == 0 or seg == self.n_segments - 1:
				tmp_atom_sc.append([sub_id, 1, 1])
			else:
				tmp_atom_sc.append([sub_id, 2, 2])
			e_id = sub_id
			#
			if seg == 0:
				bond = 0
			else:
				bond = 1
			tmp_bond[bond_id] = tuple([bond, [s_id, e_id]])
			bond_id += 1
			tmp_angle.append(e_id)
			s_id = e_id
			sub_id += 1
			
			# if self.n_sc != 0:
			# 	sc_s_id = s_id
			# 	for i in range(self.n_sc):
			# 		tmp_xyz[sub_id] = tuple(np.array(pos)  +  (i + 1)*mod_o_vec*unit_len)
			# 		tmp_atom_st.append([sub_id, 2, 1])
			# 		sc_e_id = sub_id
			# 		#
			# 		bond = 2
			# 		tmp_bond[bond_id] = tuple([bond, [sc_s_id, sc_e_id]])
			# 		sc_s_id = sc_e_id
			# 		seq_atom_id += 1
			# 		bond_id += 1
			#
		e_id = E_id
		bond = 0
		tmp_bond[bond_id] = tuple([bond, [s_id, e_id]])
		tmp_angle.append(e_id)
		bond_id += 1
		return tmp_xyz, tmp_bond, sub_id, bond_id, tmp_atom_sc, tmp_angle

	######
	def find_ortho_vec(self, list):
		vec = np.array(list).reshape(-1,1)
		# 線形独立である新たな三次元ベクトルを見つける。
		rank = 0
		while rank != 2:
			a = np.array(np.random.rand(3)).reshape(-1,1)
			target = np.hstack((vec, a))
			rank = np.linalg.matrix_rank(target)
		# QR分解により
		q, r = np.linalg.qr( target )
		# print(q[:,1])
		ortho_vec = q[:,1]
		return ortho_vec

