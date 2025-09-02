#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################
import numpy as np
######################################
##### ネットワークの初期設定を行う #####
######################################
class NWSetup:
	def __init__(self, nw_cond, target_cond):
		# self.nw_model = nw_cond[0]
		self.strand = nw_cond[1]
		# self.n_strand = nw_cond[2]
		self.n_segments = nw_cond[3]
		self.n_cell = nw_cond[4]
		self.n_sc = nw_cond[5]
		# self.l_bond = nw_cond[6]
		# self.c_n = nw_cond[7]

		self.multi = target_cond[0]

	################################################################################
	# レギュラー・ネットワーク設定
	################################################################################
	def calc_all(self):
		# 架橋点 JP を設定
		jp_xyz, strand_se_xyz = self.calc_jp_strands()
		#
		calcd_data_dic = self.set_atom(jp_xyz, strand_se_xyz)
		
		return calcd_data_dic
	################################################################################
	# JPおよびサブチェインの始点と終点のXYZを設定
	def calc_jp_strands(self):
		# jp_xyz は、JPの座標のリスト
		# strand_se_xyz は、サブチェインの出発点と終点のリスト
		if self.strand == "3_Chain_S":
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
			strand_se_xyz = [
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

		elif self.strand == "3_Chain_D":
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
			strand_se_xyz = [
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

		elif self.strand == "4_Chain":
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
			strand_se_xyz = [
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

		elif self.strand == "6_Chain":
			# JPを設定
			jp_xyz = [
			[
			[0.,0.,0.]
			]
			]
			# サブチェインの出発点と終点を設定
			strand_se_xyz = [
			[
			[[0., 0., 0.], [1, 0, 0]],
			[[0., 0., 0.], [0, 1, 0]],
			[[0., 0., 0.], [0, 0, 1]]
			]
			]

		elif self.strand == "8_Chain":
			# JPを設定
			jp_xyz = [
			[
			[0.,0.,0.],
			[0.5,0.5,0.5]
			]
			]
			# サブチェインの出発点と終点を設定
			strand_se_xyz = [
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

		return jp_xyz, strand_se_xyz

	#########################################################
	def set_atom(self, jp_xyz, strand_se_xyz):
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
				strand_xyz, bond_all, atom_sc, angle_all = self.set_strands(jp_id_dic, strand_se_xyz[mol], mol)
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
	def set_strands(self, jp_id_dic, strand_se_xyz, mol):
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
					for se_xyz in strand_se_xyz:
						tmp_xyz, tmp_bond, new_sub_id, new_bond_id, tmp_atom_sc, tmp_angle = self.calc_single_strand(jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, mol)
						strand_xyz.update(tmp_xyz)
						bond_all.update(tmp_bond)
						atom_sc.extend(tmp_atom_sc)
						angle_all.append(tmp_angle)
						sub_id = new_sub_id
						bond_id = new_bond_id
		return strand_xyz, bond_all, atom_sc, angle_all

	###############################################################
	# 一本のサブチェイン中の各アトムのxyzリストとボンドリストを作成
	def calc_single_strand(self, jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, mol):
		tmp_xyz = {}
		tmp_bond = {}
		tmp_angle = []
		tmp_atom_sc = []
		bas_xyz = np.array(b_xyz)
		# サブチェインの末端間のベクトルを設定
		start_xyz = np.array(se_xyz[0]) + bas_xyz
		end_xyz = np.array(se_xyz[1]) + bas_xyz
		vec = end_xyz - start_xyz
		# ストランドの鎖長分のループ処理
		unit_len = 1./(self.n_segments + 1)
		ortho_vec = self.find_ortho_vec(vec)
		mod_o_vec = np.linalg.norm(vec)*ortho_vec
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
			pos = tuple(start_xyz + vec*(seg + 1)/(self.n_segments + 1.))
			tmp_xyz[sub_id] = pos

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
			
			if self.n_sc != 0:
				sc_s_id = s_id
				for i in range(self.n_sc):
					tmp_xyz[sub_id] = tuple(np.array(pos)  +  (i + 1)*mod_o_vec*unit_len)
					tmp_atom_sc.append([sub_id, 2, 1])
					sc_e_id = sub_id
					#
					bond = 2
					tmp_bond[bond_id] = tuple([bond, [sc_s_id, sc_e_id]])
					sc_s_id = sc_e_id
					sub_id += 1
					bond_id += 1
			
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

