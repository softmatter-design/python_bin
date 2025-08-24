#!/usr/bin/env python
# -*- coding: utf-8 -*-
##################################################
# Import Modules
##################################################
from UDFManager import *
# ##################################################
# class Init_Chain_Select:
# 	def __init__(self, t_udf):
# 		self.uobj = UDFManager(t_udf)
# 		self.atom_list = self.uobj.get("Set_of_Molecules.molecule[].atom[]")
# 	################################################################################
# 	# 架橋点およびストランドの構成アトムのリスト
# 	def make_chain_list(self):
# 		self.uobj.jump(-1)
# 		chain_list = []
# 		end_list = []
# 		for i, target_chain in enumerate(self.atom_list):
# 			tmp = []
# 			for j, atom in enumerate(target_chain):
# 				tmp.append(j)
# 			end_list.append([i, [tmp[0], tmp[-1]]])
# 			chain_list.append([i, tmp])
# 		return end_list, chain_list

##################################################
class Init_Strand_Select:
	def __init__(self, t_udf):
		self.uobj = UDFManager(t_udf)
		self.mols = self.uobj.get("Set_of_Molecules.molecule[]")
		self.bonds = self.uobj.get("Set_of_Molecules.molecule[].bond[]")
		# self.n_atom = len(self.uobj.get("Set_of_Molecules.molecule[0].atom[]"))
	################################################################################
	# 架橋点およびストランドの構成アトムのリスト
	def make_jp_pair_list(self):
		jp_list = self.make_jp_list()
		#
		jp_pair_list = []
		strand_list = []
		for target_jp in jp_list:
			jp_pair, strand = self.make_jp_pair(target_jp, jp_list)
			for i in jp_pair:
				jp_pair_list.append(i)
			if len(strand) > 0:
				for i in strand:
					strand_list.append(i)
		return jp_list, jp_pair_list, strand_list

	# 架橋点のリストを作成
	def make_jp_list(self):
		self.uobj.jump(-1)
		jp_list = []
		#
		for i, mol in enumerate(self.mols):
			for j, atom in enumerate(mol[1]):
				tmp = []
				if atom[1] == 'JP_A' or atom[1] == 'JP_B':
					jp_list.append([i, j])
		return jp_list

	# 
	def make_jp_pair(self, target_jp, jp_list):
		molecule = target_jp[0]
		start_jp = target_jp[1]
		jp_pair = []
		strand = []
		tmp_bonds = self.bonds[molecule]
		#
		for i, bond in enumerate(tmp_bonds):
			tmp = []
			if ((bond[1] == start_jp) or (bond[2] == start_jp)) and (i < len(tmp_bonds) - 1):
				if bond[1] == start_jp:
					adj = bond[2]
				else:
					adj = bond[1]
				tmp.append(start_jp)
				tmp.append(adj)
				tmp_id = i + 1
				while tmp_bonds[tmp_id][0] == "bond_Strand":
					if tmp_bonds[tmp_id][1] == adj:
						adj = tmp_bonds[tmp_id][2]
					elif tmp_bonds[tmp_id][2] == adj:
						adj = tmp_bonds[tmp_id][1]
					tmp.append(adj)
					tmp_id += 1
				#
				if tmp_bonds[tmp_id][1] == adj:
					end_jp = tmp_bonds[tmp_id][2]
				elif tmp_bonds[tmp_id][2] == adj:
					end_jp = tmp_bonds[tmp_id][1]
				if len(tmp)>2:
					tmp.append(end_jp)
					jp_pair.append([molecule, [start_jp, end_jp]])
					strand.append([molecule, tmp])
		return jp_pair, strand