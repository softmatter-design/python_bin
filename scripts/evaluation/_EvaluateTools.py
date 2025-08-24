#!/usr/bin/env python
# -*- coding: utf-8 -*-
########################################################
# Import Modules
########################################################
from UDFManager import *
import CognacUtility as CU
from CognacBasicAnalysis import *
from CognacGeometryAnalysis import CognacGeometryAnalysis
from CognacTrajectoryAnalysis import CognacTrajectoryAnalysis
#
import numpy as np
########################################################
class Initialize:
	def __init__(self, target_udf):
		self.uobj = UDFManager(target_udf)
		self.t_udf = target_udf
		self.kg_bond = 0.97
	####################################################
	# ポリマー鎖関連の特性情報
	def calc_chain(self, chain_list, rec = 1):
		# 初期化
		self.uobj.jump(rec)
		self.bound_setup()
		CU.setCell(tuple(self.uobj.get("Structure.Unit_Cell.Cell_Size")))
		# ステップの数に対応した空リストを作成
		r2_ij = [[] for i in range(len(chain_list[0][1]))]
		#
		e2e_x = []
		e2e_y = []
		e2e_z = []
		e2e_list = []
		r2_list = []
		#
		bond_list = []
		#
		cn = []
		# 
		self.uobj.jump(rec)
		ba = CognacBasicAnalysis(self.t_udf, rec)
		for chain in chain_list:
			mol = chain[0]
			c_len = len(chain[1])
			#
			atom = self.uobj.get("Set_of_Molecules.molecule[].atom[]", [mol, chain[1][2]])[1]
			#		
			for step in range(1, c_len):
				for start in range(c_len - step):
					end1 = tuple(self.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start]]))
					end2 = tuple(self.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start + step]]))

					e2e_vec = CU.distanceWithBoundary(end1, end2)

					# e2e_vec = ba.vector([mol, chain[1][start]], [mol, chain[1][start + step]])
					e2e_dist = np.linalg.norm(np.array(e2e_vec))
					r2 = e2e_dist**2
					r2_ij[step].append(r2)
					if step == 1:
						bond_list.append(e2e_dist)
					if step == c_len -1:
						e2e_x.append(e2e_vec[0])
						e2e_y.append(e2e_vec[1])
						e2e_z.append(e2e_vec[2])
						#
						e2e_list.append(e2e_dist)
						r2_list.append(r2)
		print(np.average(r2_list))
		# gr
		cg = CognacGeometryAnalysis(self.t_udf, rec)
		gr = cg.gr([atom])
		# cn
		for i in range(1, len(r2_ij)):
			cn.append([i, np.average(np.array(r2_ij[i]))/(i*self.kg_bond**2)])
		# angle
		anglename = self.uobj.get("Molecular_Attributes.Angle_Potential[].Name")
		tmp = np.array(ba.angle(anglename[0]))
		angle_list = list(tmp[~np.isnan(tmp)])
		
		return bond_list, angle_list, e2e_x, e2e_y, e2e_z, e2e_list, r2_list, gr, cn

	def calc_nm(self):
		ct = CognacTrajectoryAnalysis(self.t_udf)
		tmp = ct.normalCoordinate("Polymer")
		res = []
		for data in tmp:
			res.append([data[0], data[1]])
		return res

	# 周期境界条件の設定
	def bound_setup(self):
		#
		axis = self.uobj.get("Simulation_Conditions.Boundary_Conditions")
		boundarylist = [0,0,0]
		#
		for i in range(0,3):
			if axis[i] == "NONE" :
				boundarylist[i] = 0
			elif axis[i] == "PERIODIC" :
				boundarylist[i] = 1
			elif axis[i] == "REFLECTIVE1" :
				boundarylist[i] = 2
			elif axis[i] == "REFLECTIVE2" :
				boundarylist[i] = 3
		CU.setBoundary(tuple(boundarylist))
		return
