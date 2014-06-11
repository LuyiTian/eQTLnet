import matplotlib.pyplot as plt
from matplotlib import cm as CM
from matplotlib import axes
from matplotlib import colors as co
def plot_heatmap(datas,gene_class):
    min_v=min(min(j) for j in datas)
    max_v=max(max(j) for j in datas)
    print "max_min",max_v,min_v
    new_datas=[]
    y_segment=[0]
    for i in range(max(gene_class)+1):
        new_datas+=[da for ind,da in enumerate(datas) if gene_class[ind]==i]
        y_segment.append(len(new_datas))
    fig = plt.figure(facecolor='w')
    ax1 = fig.add_subplot(111,position=[0.1,0.15,0.9,0.8])
    cmap = CM.get_cmap('RdYlBu_r', 100)
    map1 = ax1.imshow(new_datas, interpolation="nearest", aspect='auto',cmap=cmap, vmin=min_v,vmax=max_v)
    cb = plt.colorbar(mappable=map1, cax=None, ax=None,shrink=0.3)
    cb.set_label('(Enrichment level)')
    #ax1.set_xticks([0,10,15,20,25,30,40])
    #ax1.set_xticklabels(['-20%','TSS','25%','50%','75%','TTS','120%'])
    ax1.set_xlabel(" No. Sample \n")
    ax1.set_ylabel("No. genes")
    set_y=[]
    set_y_text=[]
    for i in range(1,len(y_segment)-1):
        ax1.plot(range(len(datas[0])),[y_segment[i]]*len(datas[0]),color='black',linewidth=1)
    for i in range(len(y_segment)-1):
        set_y.append(y_segment[i]+(y_segment[i+1]-y_segment[i])/2)
        set_y_text.append('class:'+str(i+1))
    ax1.set_yticks(set_y+[y_segment[k] for k in range(1,len(y_segment))])
    ax1.set_yticklabels(set_y_text+[str(y_segment[k]) for k in range(1,len(y_segment))])
    y_ax=plt.gca().yaxis
    for line in y_ax.get_ticklines():
        line.set_markersize(1)
        line.set_markeredgewidth(1)
    plt.ylim(0,y_segment[-1])
    plt.show()
    del new_datas

if __name__ == "__main__":
    from Pycluster import kcluster
    from parameters import processed_data_dir
    data = []
    k = 100
    repeat_times = 5
    for line in open(processed_data_dir+'Norm_raw_exp_value.txt'):
        vals = line.strip().split(' ')
        data.append([float(i) for i in vals])
    belong,sum_val,nfound=kcluster(data,k,npass=repeat_times)
    print sum_val,nfound
    plot_heatmap(data,belong)