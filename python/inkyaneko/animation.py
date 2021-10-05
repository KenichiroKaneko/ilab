# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 21:38:33 2021

@author: maeda youhei
"""

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


plt.rcParams["font.size"] = 12

date = ["210720"]
figsize = (38.4, 21.6)  # 4Kサイズ
figsize2k = (19.2, 10.8)  # 2Kサイズ
pnum = 10  # 移動平均点

path = r"python/inkyaneko/data/210720/"
shots = ["023"]  # 特定のショットだけ見るとき用、複数ショット指定可能、コメントアウトで全ショット

# センサー位置の読み込み
pos = r"python/inkyaneko/CCS_FLXLP_sensor_position.txt"

with open(pos, mode='r', encoding='utf-8') as f:
    f.readline()
    l_strip = [s.strip() for s in f.readlines()]

r = []
z = []

for row in l_strip:
    list = row.split('\t')
    r.append(float(list[0]))
    z.append(float(list[1]))


def readNf(Shot, Fn, format_file):  # Fn: NI, NJ
    if format_file == ".hdr":
        f = open(path + Fn + Shot + format_file,
                 'r', encoding="utf-8")  # ファイルを開く
        f = f.readlines()
        return np.loadtxt(f, dtype='str', skiprows=1)
    elif format_file == ".dat":
        if Fn == "NI":
            return np.fromfile(path + Fn + Shot + format_file, dtype=">h")[2:].reshape(128, 30000)
        elif Fn == "NJ":
            return np.fromfile(path + Fn + Shot + format_file, dtype=">h")[2:].reshape(96, 30000)


def offsetcal(Data, Range1, Range2):  # Data:行列形式、[t]x[ch]
    offsetdata = np.zeros((1, len(Data.T)))
    for i in np.arange(0, len(Data.T)):
        offsetdata[0, i] = np.mean(Data[Range1:Range2, i])
    for i in np.arange(0, len(Data.T)):
        Data[:, i] = Data[:, i] - offsetdata[0, i]
    return Data


def conv(Data, Point):  # Data:行列形式、[t]x[ch]、Point:移動平均点
    p = np.ones(Point)/Point
    for i in np.arange(0, len(Data.T)):
        Data[:, i] = np.convolve(Data[:, i], p, mode='same')
    return Data


# =============================================================================
"""
ファイルから情報読み込み, 必要データの抽出
"""
fdataNI = []
fdataNJ = []
for shot in shots:
    fdataNI.append(readNf(shot, "NI", ".hdr"))
    fdataNJ.append(readNf(shot, "NJ", ".hdr"))
filedata_PMT = fdataNI[0][121, :]
filedata_FL = fdataNJ[0][17, :]
filedata_MS = fdataNJ[0][53, :]
filedatas_FL = np.zeros((36, len(filedata_FL)))
filedatas_FL = fdataNJ[0][17:53, :]

chname = np.copy(filedatas_FL[:, 1])
chname[25] = chname[0]
chname = chname[1:]
# para : [point][freq [1/us]][offset][fullrange][factor]
para_PMT = np.array(filedata_PMT[2:7], dtype="float")
para_FL = np.array(filedata_FL[2:7], dtype="float")
para_MS = np.array(filedata_MS[2:7], dtype="float")
time = np.arange(0, para_PMT[0] * para_PMT[1], para_PMT[1])  # [us]
# =============================================================================


data_PMT = readNf(shot, "NI", ".dat")[121:].T
PMT_sig = data_PMT/2**15 * para_PMT[3]*para_PMT[4]
PMT_sig = offsetcal(PMT_sig, 0, 2000)

data_FL = readNf(shot, "NJ", ".dat")[17:53].T
data_FL[:, 25] = data_FL[:, 0]  # FL25が先頭に来てるので本来あるべき場所(25)に代入、
data_FL = data_FL[:, 1:]  # 先頭は代入したので削除、配列サイズ36->35
FL_sig = data_FL/2**15 * para_FL[3]*para_FL[4]
FL_sig[:, 18:] *= 10
FL_sig[:, 24] /= 10
FL_sig = offsetcal(FL_sig, 0, 200)
FL_sig = conv(FL_sig, pnum)

data_MS = readNf(shot, "NJ", ".dat")[53:93].T
MS_sig = data_MS/2**15 * para_MS[3]*para_MS[4]
MS_sig = offsetcal(MS_sig, 0, 200)
MS_sig = conv(MS_sig, pnum)

"""
FLまとめてプロット
"""
# fig, ax = plt.subplots(1, 1, figsize=figsize2k)
# ax.plot(time, FL_sig[:, 0])
# ax.set_ylim(-5, 4)  # FL26~は死んでる?から範囲はこの程度にしようかな
# ax.set_ylabel("Voltage [V]")
# ax.set_xlabel("time [$\mu$s]")
# # ax.set_xlim(8000, 10000)
# ax.set_xlim(3800, 12000)
# plt.tight_layout()
# plt.show()
# plt.close(fig)

# 玉の正規化
FL_max = max(FL_sig[:, 0])
FL_min = min(FL_sig[:, 0])
dot_max = 1000
# x = dot_max ^ (1/FL_max)
print(FL_max, FL_min)
fig, ax = plt.subplots(1, 1, figsize=figsize2k)


def update(frame):

    ax.cla()  # ax をクリア

    # ax.set_xlim(time[7600], time[24000])
    # ax.set_ylim(4, -4)
    # ax.plot(time[frame], FL_sig[frame, 0], "o")
    # ax.plot(time[7600:24000], FL_sig[7600:24000, 0])

    for ch in range(35):

        size = (FL_sig[frame, ch]+5) ** 5
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, 1)
        ax.scatter(r[ch], z[ch],  s=size, color='r')

    ax.set_xlabel('time={}[µs] data={}, size={}'.format(
        time[frame], FL_sig[frame, 0], size))
    # ax.legend('data={}, size={}'.format())


anim = FuncAnimation(fig, update, frames=range(
    18000, 20000, 100), interval=100)

plt.show()
plt.close(fig)

"""
MSまとめてプロット
"""
# fig, ax = plt.subplots(1, 1, figsize=figsize2k)
# ax.plot(time, MS_sig)
# ax.set_ylim(-5, 5)
# # ax.set_xlim(8000, 10000)
# ax.set_ylabel("Voltage [V]")
# ax.set_xlabel("time [$\mu$s]")
# plt.tight_layout()
# plt.show()
# plt.close(fig)
