from mypythontools import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from tqdm import tqdm
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
fps = 30
subs = 5
writer = FFMpegWriter(fps=fps, metadata=metadata)

matplotlib.use("Agg")

fig = plt.figure()
means, Tlist, refs = cal_means('imagesup',lowT=2.0, highT=2.4, up=True)
# means, Tlist, refs = cal_means('imagesdown',up=False)
# plt.plot(Tlist, means,label='down')
plt.plot(Tlist, refs,label='Analytical Solution')
plt.plot(Tlist, means,label='Computation')
l, = plt.plot([Tlist[0]], [means[0]],'ro',label='Current T')
xlim = plt.xlim()
ylim = plt.ylim()
plt.legend()
plt.xlabel('T')
plt.ylabel(r'$\langle M \rangle$')
# plt.show()


with writer.saving(fig, "figureM.mp4", 300):
    
    for i,(T,M) in tqdm(enumerate(zip(Tlist[::subs],means[::subs]))):
        l.set_data(T,M)
        writer.grab_frame()

make_video('./imagesup/',video_name='videouplong.mp4',subs=subs)
