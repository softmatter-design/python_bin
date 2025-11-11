#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
import os

import mod_deform_setup.DeformSetupSimple as simple
import mod_deform_setup.DeformSetupCyclic as cyclic
import mod_deform_setup.DeformSetupStep as step
import mod_deform_setup.variables as var

#######################################
# 各種UDFファイルを設定
def setup():
	print("\n\nSetting UP progress !!\n")
	if var.simple_def_mode != 'none':
		simple.setup_simple_deform()
	if var.cyclic_deform != 'none':
		cyclic.setup_cyclic_deform()
	if var.step_deform != 'none':
		step.setup_step_deform()
	return
