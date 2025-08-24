#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
# Import Modules
################################################################################
# from IObasics import MakeGraph
# from evaluation import SelectStrand
# from evaluation import EvaluateTools
#
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
from CognacTrajectoryAnalysis import CognacTrajectoryAnalysis
################################################################################
# # 必要なスクリプトを作成
# def evaluate_setup(py_mod, target_dir, f_name):
# 	script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n'
# 	script += '################################\n'
# 	script += 'import os \nimport sys \n'
# 	script += 'py_mod = "' + py_mod + '"\n'
# 	script += 'sys.path.append(py_mod)\n'
# 	script += 'from evaluation import EvaluateStrand\n'
# 	script += '################################\n'
# 	script += 'EvaluateStrand.evaluate_all()\n\n'
# 	#
# 	f_script = os.path.join(target_dir, f_name)
# 	with open(f_script, 'w') as f:
# 		f.write(script)
# 	return

###########################################################################################
def evaluate_all():
	# 対象となるUDFファイルを指定し、架橋点、そのペア、および、ストランドの構成アトムのリストを作成
	target, target_name, jp_list, strand_list, rec_size, calc_cond = select_strand()
	# ポリマー鎖関連の特性情報を計算
	eval_strand(target, target_name, strand_list, rec_size, calc_cond)
############################################################################################
def select_strand():
	# 計算対象のhonya_out.udfを選択
	target = select_outudf()
	target_name = target.split('.')[0]
	
	# 計算対象のレコードを決める
	u = UDFManager(target)
	rec_size = u.totalRecord()

	# 計算対象の条件を読み取る
	with open('network.dat', 'r') as f:
		calc_cond = f.readlines()[-1].strip().split('\t')

	# 架橋点、そのペア、および、ストランドの構成アトムのリストを作成
	ss = Init_Strand_Select(target)
	jp_list, jp_pair_list, strand_list = ss.make_jp_pair_list()

	return target, target_name, jp_list, strand_list, rec_size, calc_cond

def select_outudf():
	if len(sys.argv) == 2:
		if sys.argv[1].split("_")[-1] == "out.udf":
			target = sys.argv[1]
		else:
			sys.exit("Select hnya_out.udf !")
	else:
		sys.exit("Input file_name !")
	return target

def eval_strand(target, target_name, strand_list, rec_size, calc_cond):
	# ポリマー鎖関連の特性情報を計算
	et = InitEval(target)
	#
	n_seg = calc_cond[0]
	l_bond = calc_cond[1]
	cn_set = calc_cond[2]
	func = calc_cond[3]
	nu = calc_cond[4]
	#
	bond_list = []
	angle_list = []
	Rx_list = []
	Ry_list = []
	Rz_list = []
	R_list = []
	gr_list = []
	cn_list = []
	cn_ave = []
	cn_part = []

	if target.split('_')[0] == 'GK':
		corr = et.calc_corr()
		mm = MakeMulti(["Corr_stress", corr, ['Time', 'sigma']], target_name)
		mm.make_all()
		return

	for rec in range(1, rec_size):
		print("Reading:", rec, '/', rec_size)
		bond, angle, e2e_x, e2e_y, e2e_z, e2e, r2, gr, cn = et.calc_chain(strand_list, rec)
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
	if len(cn_list) >= 50:
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
			multi += 1
		# パートごとの平均をさらに平均
		tmp = [ [i + 1, 0] for i in range(len(cn_list[0]))]
		count = 0
		for data in cn_part:
			for i, el in enumerate(data):
				tmp[i][1] += el[1]
			count += 1
		for data in tmp:
			cn_ave.append([data[0], data[1]/count])


	# 結果をヒストグラムで出力 
	# cond_list = [base_name, data_list, bins, normalize, Legend, option]
	cond_list = [
			["bond", bond_list, 100, "True", ['bond length', 'Freq.'], 'box'],
			["angle", angle_list, 100, "True", ['angle [deg]', 'Freq.'], 'box'],
			["Rx", Rx_list, 100, "True", ['|Rx|', 'Freq.'], [n_seg, l_bond, cn_set, func] ],
			["Ry", Ry_list, 100, "True", ['|Ry|', 'Freq.'], [n_seg, l_bond, cn_set, func] ],
			["Rz", Rz_list, 100, "True", ['|Rz|', 'Freq.'], [n_seg, l_bond, cn_set, func] ],
			["R", R_list, 100, "True", ['|R|', 'Freq.'], [n_seg, l_bond, cn_set, func] ]
			]
	for cond in cond_list:
		mg = MakeHist(cond, target_name)
		mg.make_hist_all()
	# 結果をマルチグラフで出力
	cond_list = [
				["gr", gr_list, ['Distance', 'g(r)']],
				["CN", cn_list, ['|i-j|', 'C_{|i-j|}']],
				["CN_part", cn_part, ['|i-j|', 'C_{|i-j|}']],
				["CN_ave", cn_ave, ['|i-j|', 'C_{|i-j|}']]
				]
	for cond in cond_list:
		mm = MakeMulti(cond, target_name)
		mm.make_all()
	
	return

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

########################################################
class InitEval:
	def __init__(self, target_udf):
		self.uobj = UDFManager(target_udf)
		self.t_udf = target_udf
		self.kg_bond = 0.97
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
		# print(np.average(r2_list))
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

	def calc_corr(self):
		self.uobj.jump(self.uobj.totalRecord() -1)
		tmp = self.uobj.get('Correlation_Functions.Stress.Correlation[]')
		corr = []
		for data in tmp:
			corr.append([data[0], data[1]])
		return corr

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

################################################################################
class MakeHist:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, bins, normalize, Legend, option]
		self.list = cond_list[1]
		self.bins = cond_list[2]
		if target_name != '':
			self.dir = os.path.join(target_name, cond_list[0])
		else:
			self.dir = cond_list[0]

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
			script += 'C=0.1\n\n'
			#
			script += 'f(x) = C*(3/(2*pi*(N+1)*CN*bond**2))**(3/2)*exp(-3*x**2/(2*(N+1)*CN*bond**2))\n\n'
			script += 'fit f(x) data via C, CN\n\n'
			script += '#\nset label 1 sprintf("CN=%.3f", CN) at graph 0.7, 0.8\n\n'
			#
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			script += ', \\\n f(x)'

		# if self.base == "Rx":
		# 	if (type(self.option) == list) and len(self.option) == 4:
		# 		n_seg = self.option[0]
		# 		bond = self.option[1]
		# 		cn = self.option[2]
		# 		func = self.option[3]
		# 	elif (type(self.option) == list) and len(self.option) == 2:
		# 		n_seg = self.option[0]
		# 		bond = self.option[1]
		# 		cn = 1.7
		# 		func = 0
		# 	else:
		# 		n_seg = 39
		# 		bond = 0.97
		# 		cn = 1.7
		# 		func = 4
		# 	script += 'N = ' + str(n_seg) + '\n'
		# 	script += 'bond = ' + str(bond) + '\n'
		# 	script += 'CN = ' + str(cn) + '\n'
		# 	script += 'f = ' + str(func) + '\n'
		# 	#
		# 	script += 'R1 = CN*(N**0.5)*bond\n'
		# 	script += 'C=0.25\n\n'
		# 	#
		# 	if func == 3:
		# 		script += 'Pos = R1/2**0.5\ndelta = Pos*(1. - 2./f)**0.5\n\n'
		# 		script += 'f(x) = C*(1./2.)*(1./(delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*delta**2)) + exp(-1.*((x+Pos)**2)/(2.*delta**2)))\n\n'
		# 		script += 'fit f(x) data via C\n\n'
		# 	elif func == 4:
		# 		script += 'Pos = R1/3**0.5\ndelta = Pos*(1. - 2./f)**0.5\n\n'
		# 		script += 'f(x) = C*(1./2.)*(1./(delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*delta**2)) + exp(-1.*((x+Pos)**2)/(2.*delta**2)))\n\n'
		# 		script += 'fit f(x) data via C\n\n'
		# 	#
		# 	script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
		# 	script += '#\nplot data w boxes noti'
		# 	script += ', \\\n f(x)'
		#
		if self.base == "R":
			if type(self.option) == list:
				n_seg = self.option[0]
				bond = self.option[1]
				cn = 1.7
				func = 0
			else:
				n_seg = 100
				bond = 0.97
				cn = 1.7
				func = 4
			script += 'N = ' + str(n_seg) + '\n'
			script += 'bond = ' + str(bond) + '\n'
			script += 'CN = ' + str(cn) + '\n'
			script += 'f = ' + str(func) + '\n'
			script += 'C = 0.02\n\n'
			script += 'f(x, CN) = C*4.*pi*x**2.*(3./(2.*pi*(N+1)*CN*bond**2.))**(3./2.)*exp(-3.*x**2./(2.*(N+1)*CN*bond**2.))\n'	
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

		elif self.option == "box":
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			
		return script


#############################################################################################
class MakeMulti:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, Legend]
		self.list = cond_list[1]
		if target_name != '':
			self.dir = os.path.join(target_name, cond_list[0])
		else:
			self.dir = cond_list[0]

		self.base = cond_list[0]
		self.repeat = len(cond_list[1])
		#
		self.f_dat = cond_list[0] + ".dat"
		self.f_plt = cond_list[0] + ".plt"
		self.f_png = cond_list[0] + ".png"
		self.leg = cond_list[2]
	
	# マルチリストのグラフの作成
	def make_all(self):
		# データを書き出し 
		self.write_data()
		# グラフを作成
		self.make_graph()
		return

	# データを書き出し 
	def write_data(self):
		os.makedirs(self.dir, exist_ok=True)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# data:\n")
			if self.base == "CN" or self.base == "CN_part" or self.base == 'gr':
				for i, data in enumerate(self.list):
					f.write("\n\n# " + str(i) +":\n\n")
					for line in data:
						f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
			else:
				for line in self.list:
					f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

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

	# 必要なスクリプトを作成
	def make_script(self):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content()
			f.write(script)
		return

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
			if self.base == "CN" or self.base == "CN_part":
				for i in range(self.repeat):
					script += 'data ind ' + str(i) + ' w l lc ' + str(i) + ' ti "ave:' + str(i) + '", \\\n'
			else:
				script += 'data w l ti "averaged", \\\n'
			script += 'f(x) w l lw 2 ti "theory"'
		elif self.base == 'Corr_stress':
			script += 'set logscale xy \n\nset format x "10^{%L}" \nset format y "10^{%L}"\n\n'
			script += 'plot '
			script += 'data w l ti "Stress" \\\n'
		else:
			script += 'plot '
			for i in range(self.repeat):
				script += 'data ind ' + str(i) + ' w l lc ' + str(i) + 'noti, \\\n'

		return script