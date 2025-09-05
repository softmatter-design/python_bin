#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################
from mod_nw_setup import ReadNWConditions
from mod_nw_setup import NWSetup
from mod_nw_setup import SetupInitUDF
from mod_nw_setup import EquivCalcSetup
################################################################
################################################################
def main():
    ################################################################
    # 設定条件を読み込み、ネットワークポリマーの諸量を計算
    ReadNWConditions.setupcondition()

    ###################
    # ネットワークを設定
    calcd_data_dic = NWSetup.select_set()

    ##################
    # baseUDF の作成
    SetupInitUDF.setup_baseudf(calcd_data_dic)

    ###############
    # シミュレーションを設定
    EquivCalcSetup.setup_udf()

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()