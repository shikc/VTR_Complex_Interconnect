# -*- coding: utf-8 -*-

import os
import shutil
import re
import sys


class ArgvError(Exception):
    pass

work_dir = os.getcwd()
if len(sys.argv) != 2:
    raise ArgvError("Please give the src_dir runxxx")

src_dir = sys.argv[1]
source_path = work_dir + "/vtr_flow/tasks/before_pack/" + src_dir
target_path = work_dir + "/vtr_flow/benchmarks/VTR_8_blif/"

if not os.path.exists(target_path):
    os.makedirs(target_path)

if os.path.exists(source_path):
    # root 所指的是当前正在遍历的这个文件夹的本身的地址
    # dirs 是一个 list，内容是该文件夹中所有的目录的名字(不包括子目录)
    # files 同样是 list, 内容是该文件夹中所有的文件(不包括子目录)
    re_blif = re.compile("(.*).pre-vpr.blif")
    for root, dirs, files in os.walk(source_path):
        for file in files:
            if file.find('.pre-vpr.blif') > -1:
                src_file = os.path.join(root, file)
                target_file = os.path.join(target_path, re.match(re_blif, file).group(1)+".blif")
                shutil.copy(src_file, target_file)
                print(src_file)
                continue
            
            if file.find('.net') > -1:
                src_file = os.path.join(root, file)
                shutil.copy(src_file, target_path)
                print(src_file)

print('copy files finished!')
