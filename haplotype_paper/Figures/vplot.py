import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def Vasarely(obs, exp):
    # First make sure normalized
    obs_total = sum(obs.values())*1.0
    exp_total = sum(exp.values())*1.0
    for gt in obs.keys():
        obs[gt] = obs[gt]/obs_total
    for gt in exp.keys():
        exp[gt] = exp[gt]/exp_total
    
    # Get alleles
    box_w = 1; box_h = box_w
    alleles = set()
    for gt in list(obs.keys()) + list(exp.keys()):
        alleles.add(gt[0])
        alleles.add(gt[1])
    alleles = sorted(list(alleles))
    cm = plt.cm.Greys.from_list("freq",["white","black"])
    n_alleles = len(alleles)
    fig, ax = plt.subplots()

    patches = []
    colors = []
    
    # Get expected (square)
    for i in range(len(alleles)):
        for j in range(len(alleles)):
            x = i*box_w
            y = j*box_h
            rect = mpatches.Rectangle([x,y], box_w, box_h)
            patches.append(rect)
            colors.append(cm(exp.get((alleles[i], alleles[j]), 0)))
    # Get observed (circle)
    for i in range(len(alleles)):
        for j in range(len(alleles)):
            x = i*box_w
            y = j*box_h
            circ = mpatches.Circle([x+box_w/2.0,y+box_h/2.0], box_w*0.4)
            patches.append(circ)
            colors.append(cm(obs.get((alleles[i], alleles[j]), 0)))

    # Plot shapes
    collection = PatchCollection(patches, color=colors)
    ax.add_collection(collection)
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
    # Plot labels
    for i in range(len(alleles)):
        xaxis_coord = [i*box_w+box_w/2.0, -1*box_h/4.0]
        yaxis_coord = [-1*box_w/5.0, i*box_h+box_h/2.0]
        for coord in [xaxis_coord, yaxis_coord]:
            plt.text(coord[0], coord[1], alleles[i], ha="center", family='sans-serif', size=14)
    
    plt.axis('equal')
    plt.axis('off')
