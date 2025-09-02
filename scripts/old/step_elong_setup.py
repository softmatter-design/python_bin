#!/home/octa/OCTA85/Python310/bin/python
# -*- coding: utf-8 -*-
################################################################
import sys
import os
import platform
from network_setups import Step_Elong_setup

################################################################

##### Set Target #####
# 使用するCognacのバージョンを入れてください。
Ver_Cognac = "cognac112"
# ネットワーク条件のファイル名
f_data = "calc.dat"

# シミュレーションに使用するコア数
core = 1
# 応力評価スクリプト名
stress_eval = "calc_ss.py"
# 応力緩和スクリプト名
gt_eval = "calc_gt.py"

##### Conditions #####
# ステップ変形レート
step_rate = 1e-1 
# シミュレーションの時間分割
step_time_div = 0.01
# 伸長伸度
lamd_max = 2.0
# これは１ステップ計算での伸長度　Res = lambda/1_step
res = 0.1
# 
quench_list = [
			[1e-2, 10, 1],
			[1e-2, 100, 10],
			[1e-2, 1e3, 1e2],
			[1e-2, 1e4, 1e3],
			[1e-2, 1e5, 1e4],
			[1e-2, 1e6, 1e5],
			[1e-2, 1e7, 1e6],
			[1e-2, 1e8, 1e7],
		]
#
cond_list = [step_rate, step_time_div, lamd_max, res, quench_list]
#
# 計算で使用するディレクトリ
calc_dir = "Step_Elong_lambda_" + str(lamd_max ).replace('.', '_')
##### Main #####
def main():
	setup = Step_Elong_setup.Setup(py_mod, Ver_Cognac, calc_dir, f_data, core, stress_eval, gt_eval, cond_list)
	setup.calc_step()

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()
