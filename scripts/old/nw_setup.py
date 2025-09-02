#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################
import modules_1
################################################################
################################################################
def main():
################################################################
# 設定条件を読み込み、ネットワークポリマーの諸量を計算
    basic_cond, nw_cond, sim_cond, rnd_cond, target_cond, target_dir = modules_1.ReadNWConditions.setupcondition()

    ###################
    # ネットワークを設定
    nwsetup = modules_1.NWSetup.SelectSet(nw_cond, sim_cond, rnd_cond, target_dir)
    calcd_data_dic = nwsetup.select_set()

    ##################
    # baseUDF の作成
    baseudf = modules_1.SetupInitUDF.MakeInitUDF(basic_cond, sim_cond, target_cond, calcd_data_dic)
    baseudf.setup_baseudf(target_dir)
    
    ###############
    # シミュレーションを設定
    setup = modules_1.EquivCalcSetup.SetUpUDF(basic_cond, sim_cond, target_dir)
    setup.setup_udf()

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()