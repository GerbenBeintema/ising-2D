import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import cm
import cv2
from tqdm.auto import tqdm
from PIL import Image


def get_array(name, base_shape=None, resize=None):
    a = (np.fromfile(name, dtype=np.uint8)==255).astype(np.uint8)
    if base_shape is not None:
        a = Image.fromarray(a.reshape(*base_shape).T)
    else: #assume square grid
        a = Image.fromarray(a.reshape(int(len(a)**0.5),-1))
    if resize:
        a = a.resize(resize,resample=Image.NEAREST)
    return np.array(a)

def cal_means(folder_name, lowT=2.0, highT=2.4, up=True):
    means = []
    Tlist = []
    Mrefs = []
    names = sorted(os.listdir(folder_name))
    for i,name in enumerate(names):
        a = get_array(f'./{folder_name}/' + name)
        means.append(-((np.mean(a)-0.5)*2))
        if up:
            Tlist.append(lowT + i/(len(names)-1)*(highT-lowT))
        else:
            Tlist.append(highT - i/(len(names)-1)*(highT-lowT)) 
        Inner = (1-np.sinh(2/Tlist[-1])**(-4))
        if Inner < 0:
            Mrefs.append(0)
        else:
            Mrefs.append(Inner**(1/8))
    if not up:
        means = np.array(means)*np.sign(means[-1])
        
    return np.array(means), np.array(Tlist), np.array(Mrefs)

def make_video(folder_name, base_shape=None, resize=None, video_name='this_video2.mp4', subs=1):
    L = sorted(os.listdir(folder_name))
    assert len(L)>0
    
    NX, NY = get_array(folder_name+'/'+L[0], base_shape=base_shape, resize=resize).shape
    cmap = cm.get_cmap('hot')
    fps = 30
    video = cv2.VideoWriter(video_name, 0, fps, (NY, NX))
    for namei in tqdm(L[::subs]):
        ai = get_array(folder_name+'/'+namei, base_shape=base_shape, resize=resize)
        frame = (cmap(ai/1)*255).astype(np.uint8)
        video.write(frame)
    cv2.destroyAllWindows()
    video.release()

def make_plots(folder_name,base_shape=None, resize=None):

    L = sorted(os.listdir(folder_name))
    assert len(L)>0
    cmap = cm.get_cmap('hot')

    for namei in tqdm(L):
        ai = get_array(folder_name+'/'+namei, base_shape=base_shape, resize=resize)
        frame = (cmap(ai/1)*255).astype(np.uint8)
        plt.imshow(frame)
        plt.show()


def make_frames(source_video, target_folder):
    # source_video = './badapple.mp4'
    # target_folder = './sourceframes'
    import os
    try:
        os.mkdir(target_folder)
    except FileExistsError:
        pass

    vidcap = cv2.VideoCapture(source_video) 
    nframes = vidcap.get(cv2.CAP_PROP_FRAME_COUNT)
    fps = vidcap.get(cv2.CAP_PROP_FPS)
    print(f'nframe={nframes}')
    print(f'fps={fps}')

    success, image = vidcap.read()
    print(image.shape)
    print(f'shape={image.shape}')
    inc = 0

    while success:
        image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        a = (((image/255)*2-1)*127).astype(np.int8)
        a.tofile(target_folder + f'/frame{str(inc).zfill(4)}.data')
        inc += 1
        success, image = vidcap.read()