#!/usr/bin/env python
# -*- coding: utf-8 -*-
##########################
from mod_deform_setup import DeformRead
from mod_deform_setup import DeformSetup
##### Main #####
def main():
	# 変形に関する各種条件を読み取り
	DeformRead.read_all()
	# 上記条件に基づきそれぞれの設定を実行
	DeformSetup.setup()
	return
##########################
# Main
if __name__=='__main__':
	main()

