#!/home/octa/OCTA85/Python310/bin/python
# -*- coding: utf-8 -*-
################################################################################
# Import Modules
################################################################################
from UDFManager import UDFManager
import os
import sys
import numpy as np
import platform
import subprocess
#
import CognacUtility as CU
from CognacBasicAnalysis import *
from CognacGeometryAnalysis import CognacGeometryAnalysis
################################################################################
################################################################################
class EvaluateAll:
	def __init__(self):
		# 対象となるUDFファイルを指定
		self.target = sys.argv[1]
		u = UDFManager(self.target)
		self.rec_size = u.totalRecord()
		self.target_name = self.target.split('.')[0]
		# 計算対象の条件を読み取る
		with open('calc.dat', 'r') as f:
			self.calc_cond = f.readlines()[-1].strip().split('\t')
		if len(self.calc_cond) == 3:
			self.eval_target = "homo"
		elif len(self.calc_cond) == 6:
			self.eval_target = "network"

	def evaluate_all(self):
		if self.eval_target == "homo":
			# ホモポリマーのリストを作成
			sc = SelectChain(self.target)
			end_list, chain_list = sc.make_chain_list()
		elif self.eval_target == "network": 
			# ネットワークストランドのリストを作成
			ss = Init_Strand_Select(self.target)
			jp_list, jp_pair_list, chain_list = ss.make_strand_list()
		# ポリマー鎖関連の特性情報を計算
		ec = EvaluateChain(chain_list, self.rec_size, self.target, self.calc_cond, self.target_name)
		ec.eval_chain()

##############################################################################
##############################################################################
class SelectChain:
	def __init__(self, target):
		self.target = target

	# ホモポリマーのリストを作成
	def make_chain_list(self):
		uobj = UDFManager(self.target)
		atom_list = uobj.get("Set_of_Molecules.molecule[].atom[]")
		uobj.jump(-1)
		chain_list = []
		end_list = []
		for i, target_chain in enumerate(atom_list):
			tmp = []
			for j, atom in enumerate(target_chain):
				tmp.append(j)
			end_list.append([i, [tmp[0], tmp[-1]]])
			chain_list.append([i, tmp])
		return end_list, chain_list

##################################################
class Init_Strand_Select:
	def __init__(self, t_udf):
		self.uobj = UDFManager(t_udf)
		self.rec_size = self.uobj.totalRecord()
		self.mols = self.uobj.get("Set_of_Molecules.molecule[]")
		self.bonds = self.uobj.get("Set_of_Molecules.molecule[].bond[]")

	# 架橋点およびストランドの構成アトムのリスト
	def make_strand_list(self):
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

	# 架橋点どうしのペアを作成
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


###############################################################################
###############################################################################
class EvaluateChain:
	def __init__(self, chain_list, rec_size, target, calc_cond, target_name):
		self.chain_list = chain_list
		self.rec_size = rec_size
		self.target = target
		self.calc_cond = calc_cond
		self.n_seg = int(calc_cond[0])
		self.l_bond = float(calc_cond[1])
		self.target_name = target_name

	def eval_chain(self):
		bond_list = []
		angle_list = []
		Rx_list = []
		Ry_list = []
		Rz_list = []
		R_list = []
		#
		gr_list = []
		cn_list = []
		#
		if self.target.split('_')[0] == 'GK':
			corr = self.calc_corr()
			mm = MakeMulti(["Corr_stress", corr, ['Time', 'sigma']], self.target_name)
			mm.make_all()
			return

		for rec in range(1, self.rec_size):
			print("Reading Rec=", rec, '/', self.rec_size)
			bond, angle, e2e_x, e2e_y, e2e_z, e2e, r2, gr, cn = self.read_chain(rec)
			bond_list.extend(bond)
			angle_list.extend(angle)
			Rx_list.extend(e2e_x)
			Ry_list.extend(e2e_y)
			Rz_list.extend(e2e_z)
			R_list.extend(e2e)
			#
			gr_list.append(gr)
			cn_list.append(cn)
		# 鎖に沿ったセグメント間距離の平均を計算
		cn_ave, cn_part = self.calc_cn(cn_list)
		#
		self.make_output(bond_list, angle_list, Rx_list, Ry_list, Rz_list, R_list, gr_list, cn_list, cn_ave, cn_part)
		return

	# 鎖に沿ったセグメント間距離の平均を計算
	def calc_cn(self, cn_list):
		cn_ave = []
		cn_part = []
		#
		l_part = len(cn_list)//10
		# データの分割
		multi = 0
		part_cn = []
		tmp = []
		for i, part in enumerate(cn_list):
			if i < l_part*(multi + 1):
				tmp.append(part)
			else:
				part_cn.append(tmp)
				tmp = []
				tmp.append(part)
				multi += 1
		# 各パートごとに平均
		for part in part_cn:
			tmp = [ [i + 1, 0] for i in range(len(cn_list[0]))]
			count = 0
			cn_part_ave = []
			for data in part:
				for i, el in enumerate(data):
					tmp[i][1] += el[1]
				count += 1
			for data in tmp:
				cn_part_ave.append([data[0], data[1]/count])
			cn_part.append(cn_part_ave)
		# パートごとの平均をさらに平均
		tmp = [ [i + 1, 0] for i in range(len(cn_list[0]))]
		count = 0
		for data in cn_part:
			for i, el in enumerate(data):
				tmp[i][1] += el[1]
			count += 1
		for data in tmp:
			cn_ave.append([data[0], data[1]/count])
		return cn_ave, cn_part

	def make_output(self, bond_list, angle_list, Rx_list, Ry_list, Rz_list, R_list, gr_list, cn_list, cn_ave, cn_part):
		# 結果をヒストグラムで出力 
		hist_list = [
				["bond", bond_list, 200, "True", ['bond length', 'Freq.'], 'box'],
				["angle", angle_list, 200, "True", ['angle [deg]', 'Freq.'], 'box'],
				["Rx", Rx_list, 200, "True", ['|Rx|', 'Freq.'], [self.n_seg - 1, self.l_bond] ],
				["Ry", Ry_list, 200, "True", ['|Ry|', 'Freq.'], [self.n_seg - 1, self.l_bond] ],
				["Rz", Rz_list, 200, "True", ['|Rz|', 'Freq.'], [self.n_seg - 1, self.l_bond] ],
				["R", R_list, 200, "True", ['|R|', 'Freq.'], [self.n_seg - 1, self.l_bond] ]
				]
		for cond in hist_list:
			mh = MakeHist(cond, self.target_name)
			mh.make_hist_all()

		# マルチ形式での出力
		multi_list = [
				["gr", gr_list, ['Distance', 'g(r)']],
				["CN", cn_list, ['|i-j|', 'C_{|i-j|}']],
				["CN_part", cn_part, ['|i-j|', 'C_{|i-j|}']],
				["CN_ave", cn_ave, ['|i-j|', 'C_{|i-j|}']]
				]
		for cond in multi_list:
			mm = MakeMulti(cond, self.target_name)
			mm.make_all()
		return

	# ポリマー鎖関連の特性情報
	def read_chain(self, rec):
		# 初期化
		uobj = UDFManager(self.target)
		uobj.jump(rec)
		self.bound_setup()
		CU.setCell(tuple(uobj.get("Structure.Unit_Cell.Cell_Size")))
		# ステップの数に対応した空リストを作成
		r2_ij = [[] for i in range(len(self.chain_list[0][1]))]
		#
		e2e_x = []
		e2e_y = []
		e2e_z = []
		e2e_list = []
		r2_list = []
		bond_list = []
		cn = []
		#
		xp = [[] for i in range(len(self.chain_list[0][1]))]
		# 
		ba = CognacBasicAnalysis(self.target, rec)
		for chain in self.chain_list:
			mol = chain[0]
			c_len = len(chain[1])
			atom = uobj.get("Set_of_Molecules.molecule[].atom[]", [mol, chain[1][2]])[1]
			#		
			for step in range(1, c_len):
				for start in range(c_len - step):
					if len(self.calc_cond) == 3: # ポリマー鎖の場合
						e2e_vec = ba.vector([mol, chain[1][start]], [mol, chain[1][start + step]])
					elif len(self.calc_cond) == 6: # ストランドの場合
						end1 = tuple(uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start]]))
						end2 = tuple(uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start + step]]))
						e2e_vec = CU.distanceWithBoundary(end1, end2)
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
			#
			# for p in range(c_len):
			# 	xp[p].append(np.linalg.norm(ba.Xp(mol, p)))
		#
		# xp_list = []
		# for i in range(c_len):
		# 	xp_list.append([i+1, np.average(np.array(xp[i]))])
		# print(xp_list)

		# gr
		cg = CognacGeometryAnalysis(self.target, rec)
		gr = cg.gr([atom])
		# cn
		for i in range(1, len(r2_ij)):
			cn.append([i, np.average(np.array(r2_ij[i]))/(i*self.l_bond**2)])
		# angle
		anglename = uobj.get("Molecular_Attributes.Angle_Potential[].Name")
		tmp = np.array(ba.angle(anglename[0]))
		angle_list = list(tmp[~np.isnan(tmp)])
		# print(cn[:3])
		# print(r2_list[:3])
		return bond_list, angle_list, e2e_x, e2e_y, e2e_z, e2e_list, r2_list, gr, cn
	
	# 周期境界条件の設定
	def bound_setup(self):
		uobj = UDFManager(self.target)
		axis = uobj.get("Simulation_Conditions.Boundary_Conditions")
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

	def calc_corr(self):
		uobj = UDFManager(self.target)
		uobj.jump(uobj.totalRecord() -1)
		tmp = uobj.get('Correlation_Functions.Stress.Correlation[]')
		corr = []
		for data in tmp:
			corr.append([data[0], data[1]])
		return corr

##############################################################################################
##############################################################################################
class MakeHist:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, bins, normalize, Legend, option]
		self.list = cond_list[1]
		self.bins = cond_list[2]

		self.dir = os.path.join(target_name, cond_list[0])

		self.base = cond_list[0]
		self.norm = cond_list[3]
		#
		self.f_dat = cond_list[0] + "_hist.dat"
		self.f_plt = cond_list[0] + ".plt"
		self.f_png = cond_list[0] + ".png"
		self.leg = cond_list[4]
		self.option = cond_list[5]

	# ヒストグラムのグラフの作成
	def make_hist_all(self):
		# ヒストグラムのデータ作成
		bin_width, hist_data = self.make_hist_data()
		# ヒストグラムのデータを書き出し 
		self.write_data(hist_data, bin_width)
		# グラフを作成
		self.make_graph(bin_width)
		return

	# ヒストグラムのデータ作成
	def make_hist_data(self):
		# ヒストグラムを作成
		weights = np.ones(len(self.list))/float(len(self.list))
		if self.norm:
			val, x = np.histogram(self.list, bins=self.bins, weights= weights)
		else:
			val, x = np.histogram(self.list, bins=self.bins)
		# グラフ用にデータを変更
		bin_width = (x[1]-x[0])
		mod_x = (x + bin_width/2)[:-1]
		hist_data = np.stack([mod_x, val], axis = 1)
		return bin_width, hist_data

	# ヒストグラムのデータを書き出し 
	def write_data(self, hist_data, bin_width):
		os.makedirs(self.dir, exist_ok=True)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# Histgram data:\n\n")
			for line in hist_data:
				f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

	# グラフを作成
	def make_graph(self, bin_width):
		self.make_script(bin_width)
		cwd = os.getcwd()
		os.chdir(self.dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return

	# 必要なスクリプトを作成
	def make_script(self, bin_width):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content(bin_width)
			f.write(script)
		return

	# スクリプトの中身
	def script_content(self, bin_width):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		#
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		#
		script += '#\nset size square\n#set xrange [0:]\n#set yrange [0:100]\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		#
		if self.base == "Rx" or self.base == "Ry" or self.base == "Rz":
			if type(self.option) == list:
				n_seg = self.option[0]
				bond = self.option[1]
			script += 'N = ' + str(n_seg) + '\n'
			script += 'bond = ' + str(bond) + '\n'
			script += 'CN=1.7\n'
			script += 'C=1\n\n'
			#
			script += 'f(x) = C*(3/(2*pi*N*CN*bond**2))**(3/2)*exp(-3*x**2/(2*N*CN*bond**2))\n\n'
			script += 'fit f(x) data via C, CN\n\n'
			script += '#\nset label 1 sprintf("C_N=%.3f", CN) at graph 0.7, 0.8\n\n'
			#
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			script += ', \\\n f(x)'
		#
		if self.base == "R":
			if (type(self.option) == list) and len(self.option) == 4:
				n_seg = self.option[0]
				bond = self.option[1]
				cn = self.option[2]
				func = self.option[3]
			elif (type(self.option) == list) and len(self.option) == 2:
				n_seg = self.option[0]
				bond = self.option[1]
				cn = 1.7
				func = 0
			else:
				n_seg = 39
				bond = 0.97
				cn = 1.7
				func = 4
			script += 'N = ' + str(n_seg) + '\n'
			script += 'bond = ' + str(bond) + '\n'
			script += 'CN = ' + str(cn) + '\n'
			script += 'f = ' + str(func) + '\n'
			script += 'C = 0.02\n\n'
			script += 'f(x, CN) = C*4.*pi*x**2.*(3./(2.*pi*N*CN*bond**2.))**(3./2.)*exp(-3.*x**2./(2.*N*CN*bond**2.))\n'	
			script += 'fit f(x, CN) data via CN, C\n\n'
			script += '#\nset label 1 sprintf("C_N=%.3f", CN) at graph 0.7, 0.8\n\n'
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			script += ', \\\n f(x, CN)'
		#
		if self.base == "angle":
			if self.option != "box":
				script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w l noti'
			else:
				script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
				script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w boxes noti'
		#
		if self.base == "bond":
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
		#
		elif self.option == "box":
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			
		return script


#############################################################################################
class MakeMulti:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, Legend]
		self.list = cond_list[1]

		self.dir = os.path.join(target_name, cond_list[0])

		self.base = cond_list[0]
		self.repeat = len(cond_list[1])
		#
		self.f_dat = cond_list[0] + ".dat"
		self.f_plt = cond_list[0] + ".plt"
		self.f_png = cond_list[0] + ".png"
		self.leg = cond_list[2]
	
	############################################################
	# マルチリストのグラフの作成
	def make_all(self):
		# データを書き出し 
		self.write_data()
		# グラフを作成
		self.make_graph()
		return

	##############################
	# データを書き出し 
	def write_data(self):
		os.makedirs(self.dir, exist_ok=True)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# data:\n")
			if self.base == 'CN_ave' or self.base == 'Corr_stress':
				for line in self.list:
					f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
			else:
				for i, data in enumerate(self.list):
					f.write("\n\n# " + str(i) +":\n\n")
					for line in data:
						f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

	###########################################################
	# グラフを作成
	def make_graph(self):
		self.make_script()
		cwd = os.getcwd()
		os.chdir(self.dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return

	#######################
	# 必要なスクリプトを作成
	def make_script(self):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content()
			f.write(script)
		return

	#################
	# スクリプトの中身
	def script_content(self):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		#
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		#
		script += '#\nset size square\n#set xrange [1:]\n#set yrange [1:]\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		#
		if self.base == "CN" or self.base == "CN_ave" or self.base == "CN_part":
			script += '#\nset xrange [1:]\nset yrange [1:]\n'
			script += 'set key bottom\n\n'
			script += 'ct = 0.274\n'
			script += "f(x) = (1+ct)/(1-ct) -(2*ct*(1-ct**x))/(1-ct)**2/x\n\n"
			script += 'plot '
			if self.base == "CN":
				for i in range(self.repeat):
					script += 'data ind ' + str(i) + ' w l lc ' + str(i) + ' noti, \\\n'
			elif self.base == "CN_part":
				for i in range(self.repeat):
					script += 'data ind ' + str(i) + ' w l lc ' + str(i) + ' ti "part:' + str(i) + '", \\\n'
			else:
				script += 'data w l ti "averaged", \\\n'
			script += 'f(x) w l lw 2 ti "FreeRotationalModel"'
		elif self.base == 'Corr_stress':
			script += 'set logscale xy \n\nset format x "10^{%L}" \nset format y "10^{%L}"\n\n'
			script += 'plot '
			script += 'data w l ti "Stress" \\\n'
		else:
			script += 'plot '
			for i in range(self.repeat):
				script += 'data ind ' + str(i) + ' w l lc ' + str(i) + 'noti, \\\n'

		return script





















# ########################################################
# class InitEval:
# 	def __init__(self, target_udf):
# 		self.uobj = UDFManager(target_udf)
# 		self.t_udf = target_udf
# 		self.kg_bond = 0.97
# 	# ポリマー鎖関連の特性情報
# 	def calc_chain(self, chain_list, rec = 1):
# 		# 初期化
# 		self.uobj.jump(rec)
# 		self.bound_setup()
# 		CU.setCell(tuple(self.uobj.get("Structure.Unit_Cell.Cell_Size")))
# 		# ステップの数に対応した空リストを作成
# 		r2_ij = [[] for i in range(len(chain_list[0][1]))]
# 		#
# 		e2e_x = []
# 		e2e_y = []
# 		e2e_z = []
# 		e2e_list = []
# 		r2_list = []
# 		#
# 		bond_list = []
# 		#
# 		cn = []
# 		# 
# 		ba = CognacBasicAnalysis(self.t_udf, rec)
# 		for chain in chain_list:
# 			mol = chain[0]
# 			c_len = len(chain[1])
# 			#
# 			atom = self.uobj.get("Set_of_Molecules.molecule[].atom[]", [mol, chain[1][2]])[1]
# 			#		
# 			for step in range(1, c_len):
# 				for start in range(c_len - step):

# 					end1 = tuple(self.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start]]))
# 					end2 = tuple(self.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start + step]]))
# 					e2e_vec = CU.distanceWithBoundary(end1, end2)

# 					# e2e_vec = ba.vector([mol, chain[1][start]], [mol, chain[1][start + step]])

# 					e2e_dist = np.linalg.norm(np.array(e2e_vec))
# 					r2 = e2e_dist**2
# 					r2_ij[step].append(r2)
# 					if step == 1:
# 						bond_list.append(e2e_dist)
# 					if step == c_len -1:
# 						e2e_x.append(e2e_vec[0])
# 						e2e_y.append(e2e_vec[1])
# 						e2e_z.append(e2e_vec[2])
# 						#
# 						e2e_list.append(e2e_dist)
# 						r2_list.append(r2)
# 		# print(np.average(r2_list))
# 		# gr
# 		cg = CognacGeometryAnalysis(self.t_udf, rec)
# 		gr = cg.gr([atom])
# 		# cn
# 		for i in range(1, len(r2_ij)):
# 			cn.append([i, np.average(np.array(r2_ij[i]))/(i*self.kg_bond**2)])
# 		# angle
# 		anglename = self.uobj.get("Molecular_Attributes.Angle_Potential[].Name")
# 		tmp = np.array(ba.angle(anglename[0]))
# 		angle_list = list(tmp[~np.isnan(tmp)])
		
# 		return bond_list, angle_list, e2e_x, e2e_y, e2e_z, e2e_list, r2_list, gr, cn


	
