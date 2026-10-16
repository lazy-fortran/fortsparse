"""Optional benchmark visualization; numeric data only, no source matrices needed.

Usage: python plot.py trials.csv OUTPUT_DIRECTORY
CPU is the C observer's completed wait4 child user+system time. Peak RSS is
that same child's ru_maxrss; no separate peaks are summed. All requested
trials are plotted. No external-library/Pareto frontier claim is made.
"""
import argparse,csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('trials',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    with args.trials.open() as stream:rows=list(csv.DictReader(stream))
    args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':180,'svg.fonttype':'none'})
    experiments=['arbitrary-mask','uniform-full-mask','contiguous-active-run']
    cases=['kin6d-5525','drift-real-32']
    labels={'kin6d-5525':'Captured matrix, n=5525 / frozen ND',
            'drift-real-32':'Public drift fixture, n=1024 / natural'}
    # Okabe–Ito blue plus black, with redundant triangle/circle distinctions.
    styles={1:dict(color='#222222',marker='o',label='Scalar width 1'),
            8:dict(color='#0072B2',marker='^',label='Experimental width 8')}
    fig,axes=plt.subplots(3,2,figsize=(11.5,11.0),constrained_layout=True)
    for j,case in enumerate(cases):
        subset=[r for r in rows if r['case']==case]
        xs=[float(r['cpu_s'])*1000 for r in subset]
        ys=[float(r['maxrss_kib'])/1024 for r in subset]
        padx=(max(xs)-min(xs))*.12 or 1
        pady=(max(ys)-min(ys))*.18 or .1
        for i,exp in enumerate(experiments):
            ax=axes[i,j]
            group=[r for r in subset if r['experiment']==exp]
            for pair in range(1,6):
                pairrows=sorted([r for r in group if int(r['pair'])==pair],key=lambda r:int(r['width']))
                if len(pairrows)!=2:raise ValueError('missing requested pair')
                ax.plot([float(r['cpu_s'])*1000 for r in pairrows],
                        [float(r['maxrss_kib'])/1024 for r in pairrows],color='.72',linewidth=.8,zorder=1)
            for width in [1,8]:
                points=[r for r in group if int(r['width'])==width]
                if any(r['check_ok']!='True' for r in points):raise ValueError('failed trial requires explicit plot mark')
                ax.scatter([float(r['cpu_s'])*1000 for r in points],
                           [float(r['maxrss_kib'])/1024 for r in points],s=46,zorder=2,**styles[width])
            ax.set(xlabel='Complete child CPU time (ms)',ylabel='Same-child peak RSS (MiB)',
                   title=labels[case]+'\n'+['First: arbitrary-mask schedule','Second: uniform full-active path','Third: contiguous active-run path'][i],
                   xlim=(max(0,min(xs)-padx),max(xs)+padx),ylim=(min(ys)-pady,max(ys)+pady))
            ax.grid(alpha=.18)
    handles,legend=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,legend,loc='outside lower center',ncol=2,frameon=False)
    fig.suptitle('Experimental real LU panels: all 60 ordinary trials\nPaired internal controls; no external-frontier claim',fontsize=14)
    fig.savefig(args.output/'panel-cpu-rss.png')
    fig.savefig(args.output/'panel-cpu-rss.svg')
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True,sharey=True)
    for ax,case in zip(axes,cases):
        for exp,marker,color,label in [('arbitrary-mask','s','#D55E00','First: arbitrary mask'),
                                      ('uniform-full-mask','o','#0072B2','Second: uniform full-active path'),
                                      ('contiguous-active-run','^','#009E73','Third: contiguous active run')]:
            values=[]
            for pair in range(1,6):
                group=[r for r in rows if r['case']==case and r['experiment']==exp and int(r['pair'])==pair]
                bywidth={int(r['width']):r for r in group}
                values.append(float(bywidth[8]['cpu_s'])/float(bywidth[1]['cpu_s']))
            ax.plot(range(1,6),values,marker=marker,color=color,label=label,linewidth=1.1)
        ax.axhline(1,color='.3',linestyle='--',linewidth=1)
        ax.set(title=labels[case],xlabel='Alternating fresh pair',ylabel='Child CPU ratio: width 8 / scalar',xticks=range(1,6))
        ax.grid(alpha=.18)
    axes[0].legend(frameon=False,fontsize=9)
    fig.suptitle('Paired CPU observations — below 1 favors width 8',fontsize=13)
    fig.savefig(args.output/'panel-paired-cpu.png')
    fig.savefig(args.output/'panel-paired-cpu.svg')
    plt.close(fig)
if __name__=='__main__':main()
