import cv2
import glob
import os
import shutil
import time

"""
parameta
"""
frame_rate = 10.0  # FPS

"""
end
"""


def timelaps(images):

    frame = cv2.imread(images[0])

    width = int(frame.shape[1])
    height = int(frame.shape[0])

    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    video = cv2.VideoWriter('timelaps.mp4', fourcc,
                            frame_rate, (width, height))

    print("動画変換中...")

    for i in range(len(images)):
        img = cv2.imread(images[i])
        img = cv2.resize(img, (width, height))
        video.write(img)

    video.release()
    print("動画変換完了")


if __name__ == '__main__':

    start = time.time()

    images = sorted(
        glob.glob('/Users/kenichiroukaneko/ilab/python/inkyaneko/img/*.jpg'))
    print("画像の総枚数{0}".format(len(images)))
    timelaps(images)

    elapsed_time = time.time() - start
    print("処理にかかった時間は:{0}".format(elapsed_time) + "[sec]")
