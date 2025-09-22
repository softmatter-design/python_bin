#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *

import numpy as np
import os
import platform

import mod_deform_setup.variables as var

##################
# アトムのポジションを回転
def rotate_position(u, axis):
	R = rotate(axis, np.pi/2.)
	u.jump(u.totalRecord() - 1)
	pos = u.get('Structure.Position.mol[].atom[]')
	for i, mol in enumerate(pos):
		for j, atom in enumerate(mol):
			tmp = list(np.dot(np.array(R), np.array(atom)))
			u.put(tmp, 'Structure.Position.mol[].atom[]', [i, j])
	return

def rotate(axis, deg):
	if axis == 'x':
		R = [
			[1., 0., 0.],
			[0., np.cos(deg), -1*np.sin(deg)],
			[0., np.sin(deg), np.cos(deg)]
		]
	elif axis == 'y':
		R = [
			[np.cos(deg), 0., np.sin(deg)],
			[0., 1., 0.],
			[-1*np.sin(deg), 0., np.cos(deg)]
		]
	elif axis == 'z':
		R = [
			[np.cos(deg), -1*np.sin(deg), 0.],
			[np.sin(deg), np.cos(deg), 0.],
			[0., 0., 1.]
		]
	elif axis == 'yx':
		Ry = [
			[np.cos(deg), 0., np.sin(deg)],
			[0., 1., 0.],
			[-1*np.sin(deg), 0., np.cos(deg)]
		]
		Rx = [
			[1., 0., 0.],
			[0., np.cos(deg), -1*np.sin(deg)],
			[0., np.sin(deg), np.cos(deg)]
		]
		R = list(np.dot(np.array(Rx), np.array(Ry)))
	elif axis == 'zx':
		Rz = [
			[np.cos(deg), -1*np.sin(deg), 0.],
			[np.sin(deg), np.cos(deg), 0.],
			[0., 0., 1.]
		]
		Rx = [
			[1., 0., 0.],
			[0., np.cos(deg), -1*np.sin(deg)],
			[0., np.sin(deg), np.cos(deg)]
		]
		R = list(np.dot(np.array(Rx), np.array(Rz)))
	return R

###########################################
# ターミナルのタイトルを設定
def make_title(title):
	if platform.system() == "Windows":
		var.batch += "title " + title + "\n"
	elif platform.system() == "Linux":
		var.batch += r'echo -ne "\033]0; ' + title + ' \007"' + '\n'
	return

#####################
# バッチファイルを作成
def write_batchfile(dir, filename, batch_file):
	# バッチファイルを作成
	f_batch = os.path.join(dir, filename)
	with open(f_batch, 'w') as f:
		f.write(batch_file)
	if platform.system() == "Linux":
		os.chmod(f_batch, 0o744)
	return

#######################################
# サブディレクトリを使うバッチファイルを作成
def make_batch_series(subdir_list, dir, task, filename, option):
	batch_series = '#!/bin/bash\n'
	for subdir in subdir_list:
		if platform.system() == "Windows":
			batch_series += 'cd /d %~dp0\\' + subdir +'\n'
			batch_series += task
			batch_series += 'cd /d %~dp0\n'
		elif platform.system() == "Linux":
			batch_series += 'cd ./' + subdir +'\n'
			batch_series += task
			batch_series += 'cd ../\n'
	if option != '':
		batch_series += option
	write_batchfile(dir, filename, batch_series)
	return