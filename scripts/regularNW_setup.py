#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################
from octa_module import RegularNW
from octa_module import SetupInitUDF
from octa_module import EquivCalcSetup
#######################################################
## 計算条件
##################
# 使用するCognacのバージョンを入れてください。
ver_cognac = "cognac101"
blank_udf = 'cognac101.udf'
# ベースとするUDFの名前
base_udf = "base_uin.udf"
# 計算に使用するコア数
core = 2
###########
files_cond = [ver_cognac, blank_udf, base_udf, core]
#######################################################
## 計算ターゲット
###################
## Networkモデルの設定
# nw_model = "3_Chain_S"	# 三分岐　シングル
# nw_model = "3_Chain_D"	# 三分岐　ダブル
nw_model = "4_Chain"		# 四分岐
# nw_model = "6_Chain"		# 六分岐
# nw_model = "8_Chain"		# 八分岐
###################
## ポリマー鎖の設定
nw_type = "KG_entangled"		# 密度、末端間距離を設定値に合わせるように多重度を変化。 絡み合いが入るように初期化
# nw_type = "KG_simple_NPT"		# 密度、末端間距離を設定値に合わせるように多重度を変化。 絡み合いが入らないようにNPTで縮める
#
# nw_type = "KG_single"			# 設定した多重度で、密度を設定値になるように、システムサイズを縮める
# nw_type = "KG_gel"			# 設定した多重度で、密度、末端間距離を設定値に合わせるように、溶剤を添加
#
# nw_type = "Sunuke"			# 設定した多重度で末端間距離を設定値に合わせるように、密度を無視してスヌケ鎖を設定
# nw_type = "Sunuke_dens"		# 設定した多重度、設定した密度となるようにシステムサイズを縮めてスヌケ鎖を設定

#################################
# ストランドの長さ
n_segments = 20
# 一辺当たりの単位ユニット数
n_cell = 3
# ネットワークの多重度
multi_init = 1
# 設定密度
target_density = 0.85

########################################################
if nw_model == "3_Chain_S" or nw_model == "3_Chain_D":
	n_strand = 3
elif nw_model == "4_Chain":
	n_strand = 4
elif nw_model == "6_Chain":
	n_strand = 6
elif nw_model == "8_Chain":
	n_strand = 8
#
if nw_type == "KG_entangled" or nw_type == "KG_simple_NPT" or nw_type == "KG_single" or nw_type == "KG_gel":
	l_bond = 0.97
	if n_segments <= 10:
		c_n = 1.5
	elif n_segments <= 20:
		c_n = 1.65
	elif n_segments <= 40:
		c_n = 1.7
	else:
		c_n = 1.75
elif nw_type == "Sunuke" or nw_type == "Sunuke_dens":
	l_bond = 0.97
	c_n = 1.0
#
sim_cond = [nw_model, nw_type, n_segments, n_cell, multi_init, target_density, n_strand, l_bond, c_n]

####################################################
# Cognac用の名称設定
nw_name = "Network"
atom_name = ["JP_A", "End_A", "Strand_A", "Side_A", "Solvent"]
bond_name = ["bond_JP-Chn", "bond_Strand", "bond_Side"]
angle_name = ["angle_AAA"]
site_name = ["site_JP", "site_End", "site_Strand", "site_Solvent"]
pair_name = ["site_JP-site_JP", "site_Strand-site_JP", "site_Strand-site_Strand", 
				"site_JP-site_End", "site_Strand-site_End", "site_End-site_End",
				"site_Solvent-site_Solvent", "site_Solvent-site_JP", "site_Solvent-site_End",
				"site_Solvent-site_Strand"]
site_pair_name = [ 
				["site_JP", "site_JP"], 
				["site_Strand", "site_JP"], 
				["site_Strand", "site_Strand"],
				["site_JP", "site_End"], 
				["site_Strand", "site_End"], 
				["site_End", "site_End"],
				["site_Solvent", "site_Solvent"],
				["site_Solvent", "site_JP"],
				["site_Solvent", "site_End"],
				["site_Solvent", "site_Strand"],
				]

names = [nw_name, atom_name, bond_name, angle_name, site_name, pair_name, site_pair_name]

######################
##### Mainの計算 #####
######################
# import SR_Setup_UDF as Setup
import numpy as np
import sys
from UDFManager import UDFManager
import os
############
def main():
############
	###############################
	# ネットワークポリマーの諸量を計算
	init = RegularNW.InitialSetup(sim_cond)
	target_name, target_cond, target_name, calcd_data_dic = init.calc_conditions()
	##################
	# baseUDF の作成
	baseudf = SetupInitUDF.MakeInitUDF(sim_cond, target_cond, files_cond, names, target_name, calcd_data_dic)
	target_dir = baseudf.setup_baseudf()
	###############
	# 
	setup = EquivCalcSetup.SetUpUDF(nw_type, files_cond, target_name, py_mod, target_dir, names)
	setup.setup_udf()

################################################################################
if __name__=='__main__':
	main()
