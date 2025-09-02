#!/home/octa/OCTA85/Python310/bin/python
# -*- coding: utf-8 -*-
################################################################
from network_setups import Elong_setup
################################################################

##### Set Target #####
# 使用するCognacのバージョンを入れてください。
Ver_Cognac = "cognac112"
# ネットワーク条件のファイル名
f_data = "calc.dat"
# 計算で使用するディレクトリ
calc_dir = "Elong_calc"
# シミュレーションに使用するコア数
core = 6
# 応力評価スクリプト名
stress_eval = "read_stress.py"

##### Conditions #####
# これらは変形レートのリストであり、rate=lambda/tau
rate_list = [5e-4, 2e-4] # 
# シミュレーションの時間分割
time_div = 0.01
# 伸長伸度
Lamd_max = 7
# これは１ステップ計算での伸長度　Res = lambda/1_step
res = 0.02

##### Main #####
def main():
	setup = Elong_setup.Setup(Ver_Cognac, calc_dir, f_data, core, stress_eval, rate_list, time_div, Lamd_max, res)
	setup.make_all()

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()
