import os

file = "./test.py"
dirname=current_dir = os.path.basename(os.path.dirname(os.path.abspath(file)))
print(dirname)