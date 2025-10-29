#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
# Import Modules
################################################################################
import os
import sys
import math
import cmath
import numpy as np
import platform
import subprocess
import scipy.signal as signal
#
from UDFManager import UDFManager
import CognacUtility as CU
from CognacBasicAnalysis import *
from CognacGeometryAnalysis import CognacGeometryAnalysis
#
from mod_chain_evaluation import variables as var
################################################################################
# MAIN
################################################################################
def evaluate_nw():
	# ネットワークからそれぞれのストランドに対応するポリマー鎖を抽出
	chain_select()
	# ポリマー鎖関連の特性情報を計算
	eval_chain()
	# 計算結果を出力
	make_output()
	return
#
def evaluate_exc():
	file_select()
	var.exchange_list = eval_exc_strand()
	make_exc_output()
	return
################################################################################
# ネットワークからそれぞれのストランドに対応するポリマー鎖を抽出
################################################################################
def chain_select():
	# Select target UDF
	file_select()
	# Read conditions from 'target_condition.udf' and strand data from target UDF
	read_all()
	return

# 対象となる udf ファイルを選択
def file_select():
	param = sys.argv
	if len(param) == 1:
		print("usage: python", param[0], "Honya_out.udf")
		exit(1)
	elif not os.access(param[1],os.R_OK):
		print('Target file of ', param[1], "is not exist here.")
		exit(1)
	else:
		var.target = param[1]
		var.target_name = var.target.split('.')[0]
		var.uobj = UDFManager(var.target)
	return

# 計算条件から、ホモポリマーとネットワークを判断し、chain_list を読み出す。
def read_all():
	# 計算対象の条件を読み取る
	if not os.access('target_condition.udf', os.R_OK):
		print("'target_condition.udf' is not exists.")
		exit(1)
	else:
		cond_u = UDFManager('target_condition.udf')
		var.nw_type = cond_u.get('TargetCond.Model.TargetModel')
		var.func = cond_u.get('TargetCond.NetWork.N_Strands')
		var.n_seg = cond_u.get('TargetCond.NetWork.N_Segments')
		var.l_bond= cond_u.get('SimulationCond.l_bond')
		var.cn = cond_u.get('TargetCond.Strand.Characteristic_Ratio')
		var.nu = cond_u.get('TargetCond.System.Nu')
	# ネットワークストランドのリストを作成
	strand, dangling, roop, jp_pair = make_chain_list(0)
	var.chain_list = strand
	var.chain_len = len(strand[0][1][0])
	return

# 架橋点およびストランドの構成アトムのリスト
def make_chain_list(record):
	if record != 0:
		var.uobj.jump(record)
	else:
		var.uobj.jump(-1)
	strand = []
	dangling = []
	roop = []
	jp_pair =[]
	N_strand = N_dangling = N_roop = N_jp_pair = 0
	#
	mols = var.uobj.get("Set_of_Molecules.molecule[]")
	atoms = var.uobj.get("Set_of_Molecules.molecule[].atom[]")
	bonds = var.uobj.get("Set_of_Molecules.molecule[].bond[]")
	#
	prev_mol_atoms = 0
	for mol, mol_list in enumerate(mols):
		tmp_strand = []
		tmp_dangling = []
		tmp_roop = []
		tmp_jp_pair =[]
		tmp_atoms = atoms[mol]
		tmp_jp_list = [i for i in tmp_atoms if i[1] == 'JP_A']
		jp_list = []
		for line in tmp_jp_list:
			jp_list.append(line[0])
		mod_jp_list = [x-prev_mol_atoms for x in jp_list]
		# JP毎にストランドを数え上げる
		for target_jp in mod_jp_list:
			print(f'\rmol=', mol, "/", len(mols), ', jp=', mod_jp_list.index(target_jp), "/", len(mod_jp_list), "   ", end="")
			tmp_bonds = bonds[mol]
			for select in [[1,2], [2,1]]:
				targetlist = [i for i in tmp_bonds if i[select[0]] == target_jp]
				if targetlist != []:
					for each in targetlist:
						next = each[select[1]]
						reduced_bonds = [i for i in tmp_bonds if i != each]
						chain = [target_jp, next]
						while next>=0:
							next, reduced_bonds = find_next(next, reduced_bonds)
							if next != -1:
								chain.append(next)
								if next in mod_jp_list:
									break
						# chainの分類
						if chain[0] in mod_jp_list and chain[-1] in mod_jp_list:
							if chain[0] != chain[-1]:
								tmp_strand.append(chain)
								tmp_jp_pair.append([chain[0], chain[-1]])
							else:
								tmp_roop.append(chain)
						else:
							tmp_dangling.append(chain)
		# 重複を除く
		# print(len(tmp_strand), len(tmp_jp_pair))
		strand.append([mol, remove_duplicate(tmp_strand)])
		N_strand += len(strand[-1][1])
		dangling.append([mol, remove_duplicate(tmp_dangling)])
		N_dangling += len(dangling[-1][1])
		roop.append([mol, remove_duplicate(tmp_roop)])
		N_roop += len(roop[-1][1])
		jp_pair.append([mol, remove_duplicate(tmp_jp_pair)])
		N_jp_pair += len(jp_pair[-1][1])
		# 2つ目以降のmoleculeのatom数を調整
		prev_mol_atoms = len(tmp_atoms)
	print()
	print("N_strand=", N_strand, ", N_Dangling=", N_dangling, ", N_Roop=", N_roop, ", N_JP_Pair=", N_jp_pair)
	return strand, dangling, roop, jp_pair

# 隣接するアトムを探す
def find_next(target, tmp_bond_list):
	reduced_bond = tmp_bond_list
	next = target
	for select in [[1,2], [2,1]]:
		adj_list = [i for i in tmp_bond_list if i[select[0]] == target]
		if adj_list != []:
			next = adj_list[0][select[1]]
			reduced_bond = [i for i in tmp_bond_list if i != adj_list[0]]
			break
	if next == target:
		next = -1
	return next, reduced_bond

# 重複を除く
def remove_duplicate(list):
	seen = []
	tmp = [x for x in list if sorted(x) not in seen and not seen.append(sorted(x))]

	# count=0
	# tmp=[]
	# for x in list:
	# 	print(x)
	# 	if sorted(x) not in seen:
	# 		seen.append(sorted(x))
	# 		tmp.append(x)
	# 		count+=1
	# 		print(count, x)
	return tmp

###############################################################################
# ポリマー鎖関連の特性情報を計算
###############################################################################
def eval_chain():
	rec_size = var.uobj.totalRecord()
	for record in range(1, rec_size):
		print("Reading Rec=", record, '/', rec_size - 1)
		var.uobj.jump(record)
		chains = make_chains()
		make_r2_ij(chains, record)
		# read_chain()

	# 鎖に沿ったセグメント間距離の平均を計算
	calc_cn()
	#
	if var.target.split('_')[0] == 'GK':
		calc_gk()
	return

def eval_exc_strand():
	result=[]
	rec_size = var.uobj.totalRecord()
	for record in range(1, rec_size):
		print("Reading Rec=", record, '/', rec_size - 1)
		strand, dangling, roop, jp_pair = make_chain_list(record)
		n_strand=0
		for each in strand:
			n_strand+=len(each[1])
		n_dangling=0
		for each in dangling:
			n_dangling+=len(each[1])
		n_roop=0
		for each in roop:
			n_roop+=len(each[1])
		result.append([n_strand, n_dangling, n_roop])
	return result

#
def make_chains():
	chains = []
	for chain in var.chain_list:
		mol = chain[0]
		for each in chain[1]:
			tmp = []
			for atom in each:
				# print(mol, each, atom)
				segment = var.uobj.get("Structure.Position.mol[].atom[]", [mol, atom])
				tmp.append(segment)
			chains.append(tmp)
	return chains

#
def make_r2_ij(chains, record):
	bound_setup()
	CU.setCell(tuple(var.uobj.get("Structure.Unit_Cell.Cell_Size")))
	# ステップの数に対応した空リストを作成
	r2_ij = [[] for i in range(var.chain_len)]

	for chain in chains:
		for step in range(1, var.chain_len):
			for start in range(var.chain_len - step):
				end1 = tuple(chain[start])
				end2 = tuple(chain[start + step])
				e2e_vec =  CU.distanceWithBoundary(end1, end2)
				e2e_dist = np.linalg.norm(np.array(e2e_vec))
				r2 = e2e_dist**2
				r2_ij[step].append(r2)
				if step == 1:
					if e2e_dist > 2:
						print(e2e_dist)
					var.bond_list.append(e2e_dist)
				if step == var.chain_len -1:
					var.Rx_list.append(e2e_vec[0])
					var.Ry_list.append(e2e_vec[1])
					var.Rz_list.append(e2e_vec[2])
					#
					var.R_list.append(e2e_dist)
	
	# cn
	# print(var.chain_len)
	# print(len(r2_ij))
	# for i, line in enumerate(r2_ij):
	# 	if i>10:
	# 		print(i)
	# 		print(line[0:5])
	# 		print(len(line))

	cn = []
	# このシミュレーションでの平均ボンド長
	var.l_bond = np.average(np.array(var.bond_list))
	for i in range(1, len(r2_ij)):
		# 伸び切り長さで割り込んで、特性比にする。
		cn.append([i, np.average(np.array(r2_ij[i]))/(i*var.l_bond**2)])
	# print(cn)
	var.cn_list.append(cn)

	# angle
	ba = CognacBasicAnalysis(var.target, record)	
	anglename = var.uobj.get("Molecular_Attributes.Angle_Potential[].Name")
	tmp = np.array(ba.angle(anglename[0]))
	var.angle_list.extend(list(tmp[~np.isnan(tmp)]))
	return

# ポリマー鎖関連の特性情報
# def read_chain():
# 	# 初期化
# 	# ステップの数に対応した空リストを作成
# 	r2_ij = [[] for i in range(var.chain_len)]
# 	xp = [[] for i in range(var.chain_len)]
# 	# 
# 	ba = CognacBasicAnalysis(var.target, var.record)
# 	for chain in var.chain_list:
# 		mol = chain[0]
		
# 		atom = var.uobj.get("Set_of_Molecules.molecule[].atom[]", [mol, chain[1][2]])[1]

# 		#		
# 		for step in range(1, var.chain_len):
# 			for start in range(var.chain_len - step):
# 				end1 = tuple(var.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start]]))
# 				end2 = tuple(var.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start + step]]))
# 				e2e_vec = CU.distanceWithBoundary(end1, end2)
# 				e2e_dist = np.linalg.norm(np.array(e2e_vec))
# 				r2 = e2e_dist**2
# 				r2_ij[step].append(r2)
# 				if step == 1:
# 					var.bond_list.append(e2e_dist)
# 				if step == var.chain_len -1:
# 					var.Rx_list.append(e2e_vec[0])
# 					var.Ry_list.append(e2e_vec[1])
# 					var.Rz_list.append(e2e_vec[2])
# 					#
# 					var.R_list.append(e2e_dist)
# 					# r2_list.append(r2)
# 		# gr
# 		# cg = CognacGeometryAnalysis(var.target, var.record)
# 		# var.gr_list.append(cg.gr([atom]))
		
# 		# xp
# 		pos = []
# 		for i in range(var.chain_len):
# 			segment = np.array(var.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][i]]))
# 			pos.append(segment)
# 		for p in range(var.chain_len):
# 			tmp = np.zeros(3)
# 			end0 = np.array(var.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][0]]))
# 			end1 = np.array(var.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][var.chain_len - 1]]))
# 			k = np.pi*p/(var.chain_len-1)
# 			for i in range(var.chain_len):
# 				segment = np.array(var.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][i]]))
# 				tmp += segment*np.cos(k*i)
# 			tmp2 = (tmp - (end0 + end1)/2.)/var.chain_len
# 			xp[p].append(tmp2)
			

# 	xp_ave = []
# 	for p in range(var.chain_len):
# 		xp_ave.append([p, np.average(np.array(xp[p]), axis = 0)])
# 		print('cog', ba.Xp(pos, p))
		
# 	var.xp_list.append(xp_ave)
# 	print(var.xp_list)

		
# 	# cn
# 	cn = []
# 	for i in range(1, len(r2_ij)):
# 		cn.append([i, np.average(np.array(r2_ij[i]))/(i*var.l_bond**2)])
# 	var.cn_list.append(cn)
# 	# angle
# 	anglename = var.uobj.get("Molecular_Attributes.Angle_Potential[].Name")
# 	tmp = np.array(ba.angle(anglename[0]))
# 	var.angle_list.extend(list(tmp[~np.isnan(tmp)]))
# 	return











# 周期境界条件の設定
def bound_setup():
	axis = var.uobj.get("Simulation_Conditions.Boundary_Conditions")
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

















###################################################################
# 鎖に沿ったセグメント間距離の平均を計算
def calc_cn():
	l_part = len(var.cn_list)//10
	# データの分割
	multi = 0
	part_cn = []
	tmp = []
	for i, part in enumerate(var.cn_list):
		if i < l_part*(multi + 1):
			tmp.append(part)
		else:
			part_cn.append(tmp)
			tmp = []
			tmp.append(part)
			multi += 1
	# 各パートごとに平均
	for part in part_cn:
		tmp = [ [i + 1, 0] for i in range(len(var.cn_list[0]))]
		count = 0
		cn_part_ave = []
		for data in part:
			for i, el in enumerate(data):
				tmp[i][1] += el[1]
			count += 1
		for data in tmp:
			cn_part_ave.append([data[0], data[1]/count])
		var.cn_part.append(cn_part_ave)
	# パートごとの平均をさらに平均
	tmp = [ [i + 1, 0] for i in range(len(var.cn_list[0]))]
	count = 0
	for data in var.cn_part:
		for i, el in enumerate(data):
			tmp[i][1] += el[1]
		count += 1
	for data in tmp:
		var.cn_ave.append([data[0], data[1]/count])
	return



















###############################################################################
# 計算結果を出力
###############################################################################
def make_output():
	# 結果をヒストグラムで出力 
	hist_list = [
			["bond", var.bond_list, 200, "True", ['bond length', 'Freq.'], 'box'],
			["angle", var.angle_list, 200, "True", ['angle [deg]', 'Freq.'], 'box'],
			["Rx", var.Rx_list, 200, "True", ['Rx', 'Freq.'], '' ],
			["Ry", var.Ry_list, 200, "True", ['Ry', 'Freq.'], '' ],
			["Rz", var.Rz_list, 200, "True", ['Rz', 'Freq.'], '' ],
			["R", var.R_list, 200, "True", ['|R|', 'Freq.'], '' ]
			]
	for cond in hist_list:
		make_hist_all(cond)
	# マルチ形式での出力
	multi_list = [
			["gr", var.gr_list, ['Distance', 'g(r)']],
			["CN", var.cn_list, ['|i-j|', 'C_{|i-j|}']],
			["CN_part", var.cn_part, ['|i-j|', 'C_{|i-j|}']],
			["CN_ave", var.cn_ave, ['|i-j|', 'C_{|i-j|}']]
			]
	for cond in multi_list:
		make_multi(cond)
	return

# 結合交換の結果を出力
def make_exc_output():
	cond = ["Exchange", var.exchange_list, ['record', 'results']]
	make_multi(cond)
	return

##########################
# ヒストグラムのグラフの作成
def make_hist_all(cond_list):
	var.base_name = cond_list[0]
	var.data_list = cond_list[1]
	var.n_bins = cond_list[2]
	var.norm = cond_list[3]
	var.leg = cond_list[4]
	var.option = cond_list[5]
	var.target_dir = os.path.join(var.target_name, var.base_name)
	var.f_dat = var.base_name + "_hist.dat"
	var.f_plt = var.base_name + ".plt"
	var.f_png = var.base_name + ".png"

	# ヒストグラムのデータ作成
	var.bin_width, hist_data = make_hist_data()
	# ヒストグラムのデータを書き出し 
	write_data(hist_data)
	# グラフを作成
	make_graph()
	return

# ヒストグラムのデータ作成
def make_hist_data():
	# ヒストグラムを作成
	data_weights = np.ones(len(var.data_list))/float(len(var.data_list))
	if var.norm:
		value, x = np.histogram(var.data_list, bins = var.n_bins, weights = data_weights)
	else:
		value, x = np.histogram(var.data_list, bins = var.n_bins)
	# グラフ用にデータを変更
	width = (x[1] - x[0])
	mod_x = (x + width/2)[:-1]
	hist_data = np.stack([mod_x, value], axis = 1)
	return width, hist_data

# ヒストグラムのデータを書き出し 
def write_data(hist_data):
	os.makedirs(var.target_dir, exist_ok = True)
	with open(os.path.join(var.target_dir, var.f_dat), 'w') as f:
		f.write("# Histgram data:\n\n")
		for line in hist_data:
			f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
	return

# グラフを作成
def make_graph():
	make_script()
	cwd = os.getcwd()
	os.chdir(var.target_dir)
	if platform.system() == "Windows":
		subprocess.call(var.f_plt, shell = True)
	elif platform.system() == "Linux":
		subprocess.call('gnuplot ' + var.f_plt, shell = True)
	os.chdir(cwd)
	return

# 必要なスクリプトを作成
def make_script():
	with open(os.path.join(var.target_dir , var.f_plt), 'w') as f:
		script = script_content()
		f.write(script)
	return

# スクリプトの中身
def script_content():
	script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
	script += '# \ndata = "' + var.f_dat + '" \nset output "' + var.f_png + ' "\n'
	script += '#\nset size square\n'
	script += '#\nset xlabel "' + var.leg[0] + '"\nset ylabel "' + var.leg[1] + '"\n\n'
	if var.base_name == "Rx" or var.base_name == "Ry" or var.base_name == "Rz":
		script += 'N = ' + str(var.n_seg) + '\n'
		script += 'bond = ' + str(var.l_bond) + '\n'
		script += 'CN = ' + str(var.cn) + '\n'
		script += 'func = ' + str(var.func) + '\n\n'
		script += 'R1 = bond*(CN*(N+1))**0.5\n'
		script += 'C=0.1\n'
		script += 'delta=R1\n'
		script += 'frc = 1.0\n\n'
		#
		if var.nw_type == 'Regular' and var.func == 3:
			script += 'set xrange [0:]\n#set yrange [0:100]\n'
			script += 'Pos = R1/2**0.5\n\n'
			script += 'f(x) = C*(1./2.)*(1./(frc*delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*(frc*delta)**2)) + exp(-1.*((x+Pos)**2)/(2.*(frc*delta)**2)))\n\n'
			script += 'fit f(x) data via C, frc\n\n'
			script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
		elif var.nw_type == 'Regular' and var.func == 4:
			script += 'set xrange [0:]\n#set yrange [0:100]\n'
			script += 'Pos = R1/3**0.5\ndelta = Pos*(1. - 2./func)**0.5\n\n'
			script += 'f(x) = C*(1./2.)*(1./(frc*delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*(frc*delta)**2)) + exp(-1.*((x+Pos)**2)/(2.*(frc*delta)**2)))\n\n'
			script += 'fit f(x) data via C, frc\n\n'
			script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
		elif var.nw_type == 'Regular' and var.func == 6:
			script += 'set xrange [0:]\n#set yrange [0:100]\n'
			script += 'Pos = R1\ndelta = Pos*(1. - 2./func)**0.5\n\n'
			script += 'f(x) = C*(1./2.)*(1./(frc*delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*(frc*delta)**2)) + exp(-1.*((x+Pos)**2)/(2.*(frc*delta)**2)))\n\n'
			script += 'fit f(x) data via C, frc\n\n'
			script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
		elif var.nw_type == 'Regular' and var.func == 8:
			script += 'set xrange [0:]\n#set yrange [0:100]\n'
			script += 'Pos = R1/3**0.5\ndelta = Pos*(1. - 2./func)**0.5\n\n'
			script += 'f(x) = C*(1./2.)*(1./(frc*delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*(frc*delta)**2)) + exp(-1.*((x+Pos)**2)/(2.*(frc*delta)**2)))\n\n'
			script += 'fit f(x) data via C, frc\n\n'
			script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
		else:
			script += '#set xrange [0:]\n#set yrange [0:100]\nfrc=1.0\n\n'
			script += 'f(x) = C*(3/(2*pi*N*CN*bond**2))**(3/2)*exp(-3*x**2/(2*N*CN*bond**2))\n\n'
			script += 'fit f(x) data via C, CN\n\n'
			script += '#\nset label 1 sprintf("CN =%.3f", CN) at graph 0.7, 0.8\n\n'
		#
		script += 'set style fill solid 0.5\nset boxwidth ' + str(var.bin_width) + '\n'
		script += '#\nplot data w boxes noti'
		script += ', \\\n f(x)'

	if var.base_name == "R":
		script += 'N = ' + str(var.n_seg) + '\n'
		script += 'bond = ' + str(var.l_bond) + '\n'
		script += 'CN = ' + str(var.cn) + '\n'
		script += 'f = ' + str(var.func) + '\n'
		script += 'Rsq = (N+1)*bond**2\n'
		script += 'C=0.1\nfrc=1.0\n\n'
		script += '#f(x) = C*exp(-1.*(x-frc*R1)**2./(2.*sigma**2.))/(2.*pi*sigma**2.)**(1/2)\n\n'
		script += 'f(x) = C*4.*pi*x**2.*(3./(2.*pi*CN*Rsq))**(3./2.)*exp(-3.*x**2./(2.*CN*Rsq))\n'	
		script += 'fit f(x) data via C, CN\n\n'
		script += '#\nset label 1 sprintf("CN=%.3f", CN) at graph 0.7, 0.8\n'
		script += 'set style fill solid 0.5\nset boxwidth ' + str(var.bin_width) + '\n'
		script += '#\nplot data w boxes noti'
		script += ', \\\n f(x)'
	#
	if var.base_name == "angle":
		if var.option != "box":
			script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w l noti'
		else:
			script += 'set style fill solid 0.5\nset boxwidth ' + str(var.bin_width) + '\n'
			script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w boxes noti'
	#
	if var.base_name== "bond":
		script += 'set style fill solid 0.5\nset boxwidth ' + str(var.bin_width) + '\n'
		script += '#\nplot data w boxes noti'
	#
	elif var.option == "box":
		script += 'set style fill solid 0.5\nset boxwidth ' + str(var.bin_width) + '\n'
		script += '#\nplot data w boxes noti'
		
	return script


##########################
# マルチリストのグラフの作成
def make_multi(cond_list):
	var.base_name = cond_list[0]
	var.data_list = cond_list[1]
	var.leg = cond_list[2]
	var.target_dir = os.path.join(var.target_name, var.base_name)
	var.f_dat = var.base_name + "_hist.dat"
	var.f_plt = var.base_name + ".plt"
	var.f_png = var.base_name + ".png"

	# データを書き出し 
	write_multi_data()
	# グラフを作成
	make_multi_graph()
	return

# データを書き出し 
def write_multi_data():
	os.makedirs(var.target_dir, exist_ok=True)
	with open(os.path.join(var.target_dir, var.f_dat), 'w') as f:
		f.write("# data:\n")
		if var.base_name == 'CN_ave' or var.base_name == 'Corr_stress' or var.base_name == 'Corr_stress_semi' or var.base_name == 'Corr_stress_mod' or var.base_name == 'Corr_stress_all':
			for line in var.data_list:
				for data in line:
					f.write(str(data) + '\t')
				f.write('\n')
		elif var.base_name == 'Exchange':
			for i, line in enumerate(var.data_list):
				f.write(str(i) + '\t')
				for data in line:
					f.write(str(data) + '\t')
				f.write('\n')
		else:
			for i, data in enumerate(var.data_list):
				f.write("\n\n# " + str(i) +":\n\n")
				for line in data:
					f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
	return

# グラフを作成
def make_multi_graph():
	make_multi_script()
	cwd = os.getcwd()
	os.chdir(var.target_dir)
	if platform.system() == "Windows":
		subprocess.call(var.f_plt, shell=True)
	elif platform.system() == "Linux":
		subprocess.call('gnuplot ' + var.f_plt, shell=True)
	os.chdir(cwd)
	return

# 必要なスクリプトを作成
def make_multi_script():
	with open(os.path.join(var.target_dir, var.f_plt), 'w') as f:
		script = multi_script_content()
		f.write(script)
	return

# スクリプトの中身
def multi_script_content():
	repeat = len(var.data_list )
	#
	script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
	script += '# \ndata = "' + var.f_dat + '" \nset output "' + var.f_png + ' "\n'
	script += '#\nset size square\n#set xrange [1:]\n#set yrange [1:]\n'
	script += '#\nset xlabel "' + var.leg[0] + '"\nset ylabel "' + var.leg[1] + '"\n\n'
	#
	if var.base_name == "CN" or var.base_name == "CN_ave" or var.base_name == "CN_part":
		script += '#\nset xrange [1:]\nset yrange [1:]\n'
		script += 'set key bottom\n\n'
		script += 'ct = 0.274\n'
		script += "f(x) = (1+ct)/(1-ct) -(2*ct*(1-ct**x))/(1-ct)**2/x\n\n"
		script += 'plot '
		if var.base_name == "CN":
			for i in range(repeat):
				script += 'data ind ' + str(i) + ' w l lc ' + str(i) + ' noti, \\\n'
		elif var.base_name == "CN_part":
			for i in range(repeat):
				script += 'data ind ' + str(i) + ' w l lc ' + str(i) + ' ti "part:' + str(i) + '", \\\n'
		else:
			script += 'data w l ti "averaged", \\\n'
		script += 'f(x) w l lw 2 ti "FreeRotationalModel"'
	elif var.base_name == 'Corr_stress' or var.base_name == 'Corr_stress_mod':
		script += 'set logscale xy \n\nset format x "10^{%L}" \nset format y "10^{%L}"\n\n'
		script += f'G={var.nu:}\nfunc={var.func:}\nf1 = (func - 1.)/(func + 1.)\nf2 = 1. - 2./func\n\n'
		script += 'tau = 10000\nst = 0.1\neq = 0.1\ns = 1000\ne = 100000\n'
		script += 'g(x) = eq +(st-eq)*exp(-x/tau)\nfit [s:e] g(x) data via st, eq, tau\n\n'
		script += 'set label 1 sprintf("{/Symbol t} = %.2e", tau) at graph 0.15, 0.3\n'
		script += 'set label 2 sprintf("{/Symbol s}_{pre} = %.2e", st) at graph 0.15, 0.2\n'
		script += 'set label 3 sprintf("{/Symbol s}_{eq} = %.2e", eq) at graph 0.15, 0.1\n'
		script += 'set label 4 "{/Symbol s}_{nom}(t) = {/Symbol s}_{eq}+({/Symbol s}_{pre}-{/Symbol s}_{eq})*exp(-t/{/Symbol t})" at graph 0.25, 0.6\n'
		script += 'set label 5 sprintf("fitted: %.d to %.d", s, e) at graph 0.5, 0.5\n'
		script += 'set label 6 sprintf("{/Symbol n}k_BT = %.2e", G) at graph 0.5, 0.4\n\n'
		script += 'plot '
		script += 'data w l ti "Stress", \\\n'
		script += '[s:e] g(x) w l lw 2 lt 2 ti "fit", \\\n'
		script += '[1000:] G w l lw 3 dt (10, 5) lt 8 ti "Affin", \\\n'
		script += '[1000:] G*f1 w l lw 3 dt (10, 5) lt 9 ti "Q. Pht.", \\\n'
		script += '[1000:] G*f2 w l lw 3 dt (10, 5) lt 7 ti "Phantom"\n\n'
	elif var.base_name == 'Corr_stress_semi':
		script += 'set logscale y \n\n#set format x "10^{%L}" \nset format y "10^{%L}"\n\n'
		script += 'a = 1\ntau =1000\n\ns = 100\ne = 1000\n\n'
		script += 'f(x) = a*exp(-1*x/tau) \n'
		script += 'fit [s:e] f(x) data usi 1:2 via a,tau\n\n'
		script += 'set label 1 sprintf("Fitted \\nA = %.1e \\n{/Symbol t} = %.1e \\nFitting Region: %d to %d", a, tau, s, e) at graph 0.35, 0.75\n\n'
		script += 'plot '
		script += 'data w l ti "Stress", \\\n'
		script += '[s:e] f(x) noti'
	elif var.base_name == 'Corr_stress_all':
		script += 'set logscale xy \n\nset format x "10^{%L}" \nset format y "10^{%L}"\n\n'
		script += 'plot data u 1:2 w l ti "G_t", \\\n'
		script += 'data u 1:3 w l ti "xy", \\\n'
		script += 'data u 1:4 w l ti "yz", \\\n'
		script += 'data u 1:5 w l ti "zx", \\\n'
		script += 'data u 1:6 w l ti "xx-yy", \\\n'
		script += 'data u 1:7 w l ti "yy-zz"'
	elif var.base_name == 'Exchange':
		script += 'plot data u 1:2 w l ti "Strands", \\\n'
		script += 'data u 1:3 w l ti "Danglings", \\\n'
		script += 'data u 1:4 w l ti "Roops"'
	else:
		script += 'plot '
		for i in range(repeat):
			script += 'data ind ' + str(i) + ' w l lc ' + str(i) + 'noti, \\\n'

	return script


















def msd():
	samples = 10
	records = 7

	target = np.array([range(i, i+records) for i in np.arange(samples)])

	base = np.zeros(records)

	tlist = [[[-1. if i == number else 1.0 if i == number + step else x for i, x in enumerate(base)] if number < records - step else base for number in range(records - 1)] for step in range(1, records)]
	# tlist = []
	# for step in range(1, records):
	# 	tmp2 = []
	# 	for number in range(records-1):
	# 		tmp = []
	# 		for i, elm in enumerate(base):
	# 			if i == number:
	# 				tmp.append(-1.)
	# 			elif i == number+step:
	# 				tmp.append(1.)
	# 			else:
	# 				tmp.append(elm)
	# 		if number < records-step:
	# 			tmp2.append(tmp)
	# 		else:
	# 			tmp2.append(base)
	# 	tlist.append(tmp2)


	modar = np.array(tlist)


	norm_ar = np.array([1./x for x in reversed(range(1, records))])
	print(norm_ar)

	abs_d = np.abs(np.matmul(modar, np.transpose(target)))
	sum_data = np.sum(np.average(abs_d, axis = 2), axis = 1)
	ave_abs = np.multiply(sum_data, norm_ar)
	print(ave_abs)

	sqred_d = np.square(np.matmul(modar, np.transpose(target)))
	sum_sq_data = np.sum(np.average(sqred_d, axis = 2), axis = 1)
	ave_sq = np.multiply(sum_sq_data, norm_ar)
	print(ave_sq)












####################################################################################
# Green Kubo での計算を処理
def calc_gk():
	corr, corr_all = calc_corr()
	cond_list = ["Corr_stress", corr, ['Time', 'sigma']]
	make_multi(cond_list)
	cond_list = ["Corr_stress_all", corr_all, ['Time', 'sigma', 'ave']]
	make_multi(cond_list)

	irheo(corr)
	return

def calc_corr():
	var.uobj.jump(var.uobj.totalRecord() - 1)
	#
	vol = var.uobj.get('Statistics_Data.Volume.Total_Average')
	corr_all = var.uobj.get('Correlation_Functions.Stress.Correlation[]')
	corr = []
	prev = 0.
	for data in corr_all:
		time = data[0]
		ave = vol*np.average(np.array(data[2:]))
		if data[1] > 0:
			g = data[1]
			prev = data[1]
		else:
			g = prev
		corr.append([time, g, ave])
	# g_mod = signal.savgol_filter(g, 3, 2)
	# corr_mod = np.stack([time, g_mod], 1)
	return corr, corr_all


##################################
# 
def irheo(data_list):
	minmax = [1e-5, 1e2]
	div = 10
	#
	# mod = self.modify(data_list)
	# gt, mod_gt = self.modify_data(mod)
	#
	gw = calcgw(data_list, minmax, div)
	save_gw_data(gw, 'gw.dat')
	#
	# self.save_data(data_list, 'modified.dat')
	# self.plotgtgw('modified.dat')
	# cmd = "corr2gw < modified.dat > gw.dat"
	# subprocess.call(cmd, shell=True)

	plotgtgw('gw.dat')
	#
	return

def modify_data(data_list):
	fine_div = 100
	#
	glist = []
	timelist = []
	for data in data_list:
		time = data[0]
		g = data[1]
		if time == 0.0:
			timelist.append(time)
			glist.append(g)
		else:
			for i in range(1, fine_div + 1):
				timelist.append(pre_time + i*(time-pre_time)/fine_div)
				glist.append(pre_g + i*(g-pre_g)/fine_div)
		pre_time = time
		pre_g = g
		#
	mod_g = signal.savgol_filter(glist, 5, 3)
	#
	gt = np.stack([timelist, glist], 1)
	mod_gt = np.stack([timelist, mod_g], 1)
	#
	return gt, mod_gt

def calcgw(gt, minmax, div):
	gw = []
	mag = math.log10(minmax[0])
	while mag < math.log10(minmax[1]):
		for i in range(div):
			omega = 10**(mag+i/div)
			gstar = gs(gt, omega)
			gw.append([omega, gstar.real, abs(gstar.imag)])
		mag += 1
	#
	return gw

def gs(gt, omega):
	gstar = gt[0][1] + (1 - cmath.exp(-1j*omega*gt[1][0]))*(gt[1][1] - gt[0][1])/gt[1][0]/(1j*omega)
	for k in range(len(gt) - 2):
		gstar += (gt[k+2][1] - gt[k+1][1])*(cmath.exp(-1j*omega*gt[k+1][0]) - cmath.exp(-1j*omega*gt[k+2][0]))/(gt[k+2][0] - gt[k+1][0])/(1j*omega)
	#
	return gstar 

#----- 計算結果をターゲットファイル名で保存
def save_gw_data(target, f_data):
	with open(f_data,'w') as f:
		for line in target:
			for data in line:
				f.write(str(data) + '\t')
			f.write('\n')
	return

#----- 結果をプロット
def plotgtgw(f_data):
	plt = make_gtgw(f_data)
	#
	if platform.system() == "Windows":
		subprocess.call([plt], shell=True)
	elif platform.system() == "Linux":
		subprocess.call(['gnuplot ' + plt], shell=True)
	return

# 必要なスクリプトを作成
def make_gtgw(f_data):
	script = gtgw_content(f_data)
	plt = f_data.replace('dat', 'plt')
	with open(plt, 'w') as f:
		f.write(script)
	return plt

# スクリプトの中身
def gtgw_content(f_data):
	out_png = f_data.replace('dat', 'png')
	script = 'set term pngcairo font "Arial,14"\n\n'
	script += 'set colorsequence classic\n\n'
	script += 'data = "' + f_data + '"\n\n'
	script += 'set output "' + out_png + '"\n\n'
	script += 'set key left\nset size square\n'
	script += '#set xrange [1:4]\n#set yrange [0:0.2]\n#set xtics 1\n#set ytics 0.1\n'

	if f_data == 'modified.dat' or f_data == 'ave_all_stress.dat':
		script += 'set logscale xy\n'
		script += 'set format x "10^{%L}" \nset format y "10^{%L}"\n'
		script += 'set xlabel "Time"\nset ylabel "Stress"\n'	
		script += 'plot	data u 1:2 axis x1y1 w l lw 2 lt 1 ti "Stress"'
	elif f_data == 'gw.dat' or f_data == 'freq_mod.dat':
		script += 'set xrange [:1e2]\nset yrange [1e-4:]\nset y2range [1e-1:1e1]\nset y2tics\n'
		script += 'set logscale xyy2\n'
		script += '# 斜辺の傾きが -2 の三角形の準備\n'
		script += 'a = 30; # グラフの中に入るように三角形の高さを調整\n'
		script += 'x1=5e-4; x2=1e-3;\n'
		script += 'y1=a*x1**(1);y2=a*x2**(1);\n'
		script += 'set object 1 polygon from x1,y1 to x2,y1 to x2,y2 to x1,y1 fs empty border\n\n'
		script += 'set format x "10^{%L}" \nset format y "10^{%L}"\nset format y2 "10^{%L}"\n'
		# script += 'set label 1 sprintf("{/Symbol l} = %.1f", deform) at graph 0.6, 0.9\n\n'
		script += 'set xlabel "Frequency"\nset ylabel "G' + "', G''" + '"\nset y2label "tan{/Symbol d}"\n\n'
		script += 'plot	data u 1:2 w lp lt 1 ti "G' + "'" + '", \\\n'
		script += 'data u 1:3 w lp lt 2 ti "G' + "''" + '", \\\n'
		script += 'data u 1:($3/$2) axis x1y2 w lp lt 3 ti "tan{/Symbol d}"'
	script += '\n\nreset'

	return script

def modify(data_list):
	a = 0.057
	tau = 190
	fitstart = 500
	mod_gt = []
	for data in data_list:
		time = float(data[0])
		g = float(data[1])
		if time < fitstart:
			# if g > 0:
			mod_gt.append([time, g])
		else:
			break
	time = fitstart
	while time < 1e5:
		tmp = a*np.exp(-time/tau)
		if tmp > 1e-10:
			mod_gt.append([time, tmp])
			time += 10**int(np.log10(time))/100
		else:
			break
	# save_data(mod_gt, 'mod_gt.dat')
	return mod_gt









def tttt():
	
	samples = 10
	records = 7

	target = np.array([range(i, i+records) for i in np.arange(samples)])

	base = np.zeros(records)

	tlist = [[[-1. if i == number else 1.0 if i == number + step else x for i, x in enumerate(base)] if number < records - step else base for number in range(records - 1)] for step in range(1, records)]
	# tlist = []
	# for step in range(1, records):
	# 	tmp2 = []
	# 	for number in range(records-1):
	# 		tmp = []
	# 		for i, elm in enumerate(base):
	# 			if i == number:
	# 				tmp.append(-1.)
	# 			elif i == number+step:
	# 				tmp.append(1.)
	# 			else:
	# 				tmp.append(elm)
	# 		if number < records-step:
	# 			tmp2.append(tmp)
	# 		else:
	# 			tmp2.append(base)
	# 	tlist.append(tmp2)


	modar = np.array(tlist)


	norm_ar = np.array([1./x for x in reversed(range(1, records))])
	print(norm_ar)

	abs_d = np.abs(np.matmul(modar, np.transpose(target)))
	sum_data = np.sum(np.average(abs_d, axis = 2), axis = 1)
	ave_abs = np.multiply(sum_data, norm_ar)
	print(ave_abs)

	sqred_d = np.square(np.matmul(modar, np.transpose(target)))
	sum_sq_data = np.sum(np.average(sqred_d, axis = 2), axis = 1)
	ave_sq = np.multiply(sum_sq_data, norm_ar)
	print(ave_sq)


def array_test():
	t=5
	tlist = []
	base = np.zeros(t)
	for i, elm in enumerate(base):
		if i == 0:
			tlist.append(-1.)
		elif i == 1:
			tlist.append(1.)
		else:
			tlist.append(elm)

	tlist2 = [-1. if i ==0 else x for i, x in enumerate(base)]

	print(tlist)
	print(tlist2)
	# test = np.array()

	return

