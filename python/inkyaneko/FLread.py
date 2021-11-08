# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 21:38:33 2021

@author: maeda youhei
"""


import os
import glob
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 20

utstpath = r"\\TEETH\utst"
# dates = os.listdir("data/")
dates = ["211018"]  # 実験日程指定、とめても可能
# dataposi = "Server"
dataposi = "Local"
# dataposi = "Storage"    #dataの場所、Localの場合、フォルダ内のdataを参照
# figsave = "ys"
figsave = "no"  # figをsaveするときは"y"、任意場所の場合は"ys"
# figsave = "n"
figsize = (38.4, 21.6)  # 4Kサイズ
figsize2k = (19.2, 10.8)  # 2Kサイズ
pnum = 10  # 移動平均点
for date in dates:
    if dataposi == "Local":
        path = "python/inkyaneko/data/" + date + "/"
    elif dataposi == "Storage":
        path = r"D:\東大修士研究関係\実験データ過去\data_storage/utst/" + date + '/'
    elif dataposi == "Server":
        path = utstpath + "/{}/".format(date)
    files = sorted(glob.glob(path + '/*.log'))
    shotnum = int(len(files))
    # shotnum = 13
    # shots = ["007", "008", "010", "011", "012"]
    shots = []
    for i in np.arange(1, shotnum+1):
        if i < 10:
            shots.append("00{}".format(i))
        if i >= 10:
            shots.append("0{}".format(i))

    shots = ["004", "005", "006"]  # 特定のショットだけ見るとき用、複数ショット指定可能、コメントアウトで全ショット

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

    # =============================================================================

    if figsave == "y":
        fig_FLMS = "python/inkyaneko/fig/FLMSfig"
        if not os.path.exists(fig_FLMS):
            os.mkdir(fig_FLMS)
        fig_FL = fig_FLMS + "/FL"
        if not os.path.exists(fig_FL):
            os.mkdir(fig_FL)
        fig_MS = fig_FLMS + "/MS"
        if not os.path.exists(fig_MS):
            os.mkdir(fig_MS)
    if figsave == "ys":
        fig_FLMS = r"D:\東大修士研究関係\実験データ過去\fig_storage" + "/FLMSfig"  # 任意保存場所
        if not os.path.exists(fig_FLMS):
            os.mkdir(fig_FLMS)
        fig_FL = fig_FLMS + "/FL"
        if not os.path.exists(fig_FL):
            os.mkdir(fig_FL)
        fig_MS = fig_FLMS + "/MS"
        if not os.path.exists(fig_MS):
            os.mkdir(fig_MS)
    # =============================================================================
    for shot in shots:

        data_PMT = readNf(shot, "NI", ".dat")[121:].T
        PMT_sig = data_PMT/2**15 * para_PMT[3]*para_PMT[4]
        PMT_sig = offsetcal(PMT_sig, 0, 200)

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
        FLチャンネル分複数プロット
        """
        # one = ax.ravel()  # ベクトル一次元配列
        # fig, ax = plt.subplots(5, 7, figsize=figsize, sharex=True, sharey=True)
        # fig.suptitle('fuckin flux loop {}#{} '.format(date, shot))
        # for ax, i in zip(one, np.arange(0, len(FL_sig.T))):
        #     ax.plot(time, FL_sig[:, i],
        #             color=cm.viridis(i/len(FL_sig.T)), label="{}".format(chname[i]))
        #     # ax.set_xlim(8800, 10000)
        #     ax.set_ylim(-5, 4)  # FL26~は死んでる?から範囲はこの程度にしようかな
        #     ax.legend(loc="upper left")
        #     if i == 14:
        #         ax.set_ylabel("Voltage [V]")
        #     if i == 31:
        #         ax.set_xlabel("time [$\mu$s]")
        # plt.tight_layout()
        # plt.show()
        # if figsave == "y":
        #     plt.savefig(fig_FL + "/{}#{}_sig35".format(date, shot))
        # elif figsave == "ys":
        #     plt.savefig(fig_FL + "/{}#{}_sig35".format(date, shot))
        # plt.close(fig)

        # """
        # FLまとめてプロット
        # """
        fig, ax = plt.subplots(1, 1, figsize=figsize2k)
        ax.plot(time, FL_sig)
        ax.set_ylim(-5, 4)  # FL26~は死んでる?から範囲はこの程度にしようかな
        ax.set_ylabel("Voltage [V]")
        ax.set_xlabel("time [$\mu$s]")
        ax.set_xlim(8000, 10000)
        plt.tight_layout()
        plt.title("shot#{}".format(shot))
        plt.show()
        if figsave == "y":
            plt.savefig(fig_FL + "/{}#{}_all".format(date, shot))
        elif figsave == "ys":
            plt.savefig(fig_FL + "/{}#{}_all".format(date, shot))
        plt.close(fig)

        # continue

        """
        FL５個ずつプロット
        """
        fig, ax = plt.subplots(4, 2, figsize=figsize, sharex=True, sharey=True)
        one = ax.ravel()  # ベクトル一次元配列

        for ax, i in zip(one, np.arange(0, len(FL_sig.T)//5)):
            j = i*5
            ax.plot(time, FL_sig[:, j], time, FL_sig[:, j+1],
                    time, FL_sig[:, j+2], time, FL_sig[:, j+3],
                    time, FL_sig[:, j+4])
            ax.set_xlim(8000, 12000)
            ax.set_ylim(-4, 4)  # FL26~は死んでる?から範囲はこの程度にしようかな
            ax.legend(["{0}".format(chname[j]), "{0}".format(chname[j+1]),
                      "{0}".format(chname[j+2]), "{0}".format(chname[j+3]),
                       "{0}".format(chname[j+4])], fontsize=12, loc="upper left")
            # ax.legend()
        plt.show()

        # for ax, i in zip(one, np.arange(0, len(FL_sig.T))):
        #     ax.plot(time, FL_sig[:, i], label="{}".format(chname[i]))
        #     # ax.set_xlim(8800, 10000)
        #     ax.set_ylim(-5, 4)  # FL26~は死んでる?から範囲はこの程度にしようかな
        #     ax.legend(loc="upper left")

        #     # if i == 14:
        #     #     ax.set_ylabel("Voltage [V]")
        #     # if i == 31:
        #     #     ax.set_xlabel("time [$\mu$s]")

        # """
        # MSチャンネル分複数プロット
        # """
        # fig, ax = plt.subplots(5, 8, figsize=figsize, sharex=True, sharey=True)
        # one = ax.ravel()  # ベクトル一次元配列
        # fig.suptitle('fuckin ms loop {}#{} '.format(date, shot))
        # for ax, i in zip(one, np.arange(0, len(MS_sig.T))):
        #     ax.plot(time, MS_sig[:, i],
        #             color=cm.viridis(i/len(MS_sig.T)), label="MS{}".format(i+1))
        #     # ax.set_xlim(8800, 10000)
        #     ax.set_ylim(-5, 5)
        #     ax.legend(loc="upper left")
        #     if i == 16:
        #         ax.set_ylabel("Voltage [V]")
        #     if i == 35:
        #         ax.set_xlabel("time [$\mu$s]")
        # plt.tight_layout()
        # plt.show()
        # if figsave == "y":
        #     plt.savefig(fig_MS + "/{}#{}_sig40".format(date, shot))
        # elif figsave == "ys":
        #     plt.savefig(fig_MS + "/{}#{}_sig40".format(date, shot))
        # plt.close(fig)

        # """
        # MSまとめてプロット
        # """
        # fig, ax = plt.subplots(1, 1, figsize=figsize2k)
        # ax.plot(time, MS_sig)
        # ax.set_ylim(-5, 5)
        # ax.set_xlim(8000, 10000)
        # ax.set_ylabel("Voltage [V]")
        # ax.set_xlabel("time [$\mu$s]")
        # plt.tight_layout()
        # plt.show()
        # if figsave == "y":
        #     plt.savefig(fig_MS + "/{}#{}_all".format(date, shot))
        # elif figsave == "ys":
        #     plt.savefig(fig_MS + "/{}#{}_all".format(date, shot))
        # plt.close(fig)
        """
        MS５個ずつプロット
        """
        fig, ax = plt.subplots(4, 2, figsize=figsize, sharex=True, sharey=True)
        one = ax.ravel()  # ベクトル一次元配列

        for ax, i in zip(one, np.arange(0, len(MS_sig.T)//5)):
            if i == len(MS_sig.T)//5-1:
                n = len(MS_sig) % 5
            else:
                n = 5
            j = i*5

            for k in range(n):
                ax.plot(time, MS_sig[:, j+k], label="MS{:0=2}".format(j+k+1))
            # ax.plot(time, MS_sig[:, j], time, MS_sig[:, j+1],
            #         time, MS_sig[:, j+2], time, MS_sig[:, j+3],
            #         time, MS_sig[:, j+4])
            ax.set_xlim(8000, 12000)
            ax.set_ylim(-4, 4)  # FL26~は死んでる?から範囲はこの程度にしようかな
            # ax.legend(["{0}".format(chname[j]), "{0}".format(chname[j+1]),
            #           "{0}".format(chname[j+2]), "{0}".format(chname[j+3]),
            #            "{0}".format(chname[j+4])], fontsize=12, loc="upper left")
            ax.legend(fontsize=12, loc="upper left")
            # ax.legend()
        plt.show()
