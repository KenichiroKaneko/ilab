import shutil
import os
import glob


dirnames = ["z2033_r602", "z1000_r602", "z0960_r602", "z0920_r602" ,
             "z0880_r602", "z0840_r602", "z0800_r602", "z0760_r602",
             "z0720_r602", "z0680_r602", "z0640_r602", "z0600_r602",
             "z0560_r602", "z0520_r602", "z0480_r602", "z0440_r602",
             "z0400_r602", "z0360_r602", "z0320_r602", "z0280_r602",
             "z0240_r602", "z0200_r602", "z0160_r602", "z0120_r602",
             "z0080_r602"]

# GSPOLTからCCSへのコピー
# filenames = ["CCS.txt", "FluxProfile_2D.txt", "jeddy.txt", "Sensor_B.txt", "Sensor_Flux.txt"]
# for i in dirnames:
#     for file in filenames:
#         copyfrom = "i/GSPLOT_OUTPUT/" + i + "/" + file
#         copyto = "CCS_Set_NEWNEW/UTST_numel_" + i[1:5] + "/" + file

#         shutil.copy(copyfrom, copyto)

# inputファイルの作成
for i in dirnames:
    path1 = "CCS_Set_NEWNEW/CCS_input/input_0000.txt"
    path2 = "CCS_Set_NEWNEW/CCS_input/input_{0}.txt".format(i[1:5])

    with open(path1, mode="r", encoding="utf-8")as f:
        data_lines = f.read()
        data_lines = data_lines.replace("UTST_numel_0000", "UTST_numel_" + i[1:5])

    with open(path2, mode="w", encoding="utf-8")as f:
        f.write(data_lines)