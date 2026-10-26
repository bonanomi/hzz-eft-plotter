#!/usr/bin/env python3
"""Combine spectra from plot_pth_eft.py yields.csv with two Combine scan contours.

Dependencies: numpy matplotlib mplhep uproot. No PyROOT required.
Use CSV produced by --quantity counts/density for events; for xsec quantities,
use CSV produced by --quantity xsec. CSV always stores integrated bin contents.
All requested benchmarks must have been evaluated by the workspace script.

Use --coefficients chwb chwbtil (default: chg chgtil).
Scan assumptions: rectangular --algo grid scans of the selected coefficients. Grid rows
have quantileExpected >= 0; special best-fit/generation rows are excluded.
Contours use 2*deltaNLL as stored, relative to each scan's fitted reference;
no renormalization to the sampled grid minimum. Duplicate coordinates use the
lowest valid NLL. Nonfinite/failed points are masked, not extrapolated. Missing
whole grid rows/columns cannot be reconstructed; check the original scan logs.
The script cannot diagnose convergence from these four branches alone.

Default: 95% joint two-parameter asymptotic contours (q=5.99146).
Rate/shape totals include all CSV bins, even overflow omitted from density plots.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import warnings
import numpy as np


def read_scan(path, coefficients=('chg','chgtil')):
    import uproot
    with uproot.open(path) as f:
        tree=f['limit']
        a={n:tree[n].array(library='np') for n in (*coefficients,'deltaNLL','quantileExpected')}
    x,y,d,q=(a[n] for n in (*coefficients,'deltaNLL','quantileExpected'))
    rows=(q>=0)&np.isfinite(x)&np.isfinite(y)
    x,y,d=x[rows],y[rows],d[rows]
    if not len(x): raise ValueError(f'{path}: no grid rows')
    # Round float branch noise without changing practical scan resolution.
    x=np.round(x,9);y=np.round(y,9)
    xs=np.unique(x);ys=np.unique(y)
    if len(xs)<2 or len(ys)<2: raise ValueError(f'{path}: not a 2D grid')
    if len(xs)*len(ys)>5_000_000:
        raise ValueError('Scan is not a recognizable rectangular grid or is too large')
    z=np.full((len(ys),len(xs)),np.nan)
    for xx,yy,dd in zip(x,y,d):
        if not np.isfinite(dd) or dd>=1e6: continue
        if dd < -1e-3:
            raise ValueError(f'{path}: negative deltaNLL={dd}; check reference fit before plotting')
        i=np.searchsorted(ys,yy);j=np.searchsorted(xs,xx)
        zz=max(0.,2*dd)
        if not np.isfinite(z[i,j]) or zz<z[i,j]: z[i,j]=zz
    if not np.isfinite(z).any(): raise ValueError(f'{path}: no finite scan values')
    missing=int(np.isnan(z).sum())
    if missing: warnings.warn(f'{path}: masking {missing} missing/invalid grid points')
    return xs,ys,np.ma.masked_invalid(z)


def read_spectra(path,pairs,quantity,coefficients=('chg','chgtil')):
    with open(path) as f: rows=list(csv.DictReader(f))
    column='cross_section' if quantity in ('xsec','integrated-xsec') else 'events'
    if not rows or column not in rows[0]:
        raise ValueError(f'CSV needs {column!r}; generate it with the corresponding workspace quantity')
    if any(c not in rows[0] for c in coefficients):
        raise ValueError('CSV coefficient columns do not match --coefficients: '+', '.join(coefficients))
    smrows=[r for r in rows if r['group']=='SM']
    if not smrows: raise ValueError('CSV has no SM rows')
    obs={r['observable'] for r in smrows}
    if len(obs)!=1: raise ValueError('CSV must contain exactly one observable')
    def key(r): return (float(r['bin_low']),float(r['bin_high']),r['bin_label'])
    smrows.sort(key=key);bins=[key(r) for r in smrows]
    if len(set(bins))!=len(bins): raise ValueError('Duplicate SM bins')
    sm=np.array([float(r[column]) for r in smrows])
    curves=[]
    for g,t in pairs:
        matching=[r for r in rows if r['group']=='joint' and
                  np.isclose(float(r[coefficients[0]]),g,rtol=1e-8,atol=1e-10) and
                  np.isclose(float(r[coefficients[1]]),t,rtol=1e-8,atol=1e-10)]
        matching.sort(key=key)
        if [key(r) for r in matching]!=bins:
            raise ValueError(f'Missing/duplicate joint benchmark ({g},{t}); rerun workspace script with --pairs')
        curves.append(np.array([float(r[column]) for r in matching]))
    if any(not np.isfinite(v).all() or np.any(v<0) or v.sum()<=0 for v in [sm]+curves):
        raise ValueError('Spectra must have finite nonnegative bins and positive totals')
    return next(iter(obs)),bins,sm,curves


def render(a,scans,obs,bins,sm,curves,effect):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import mplhep as hep
    from matplotlib.lines import Line2D
    hep.style.use('CMS')
    base=sm.copy();values=[v.copy() for v in curves]
    if effect=='rate_only': values=[sm*(v.sum()/sm.sum()) for v in curves]
    elif effect=='shape_only': base=sm/sm.sum();values=[v/v.sum() for v in curves]
    density=a.quantity in ('density','xsec')
    widths=np.array([hi-lo if np.isfinite(hi) else (a.overflow_width or np.nan) for lo,hi,_ in bins])
    keep=np.isfinite(widths)&(widths>0) if density else np.ones(len(bins),dtype=bool)
    if not keep.any(): raise ValueError('No finite-width bins; supply --overflow-width or use integrated quantity')
    if density:
        base=base[keep]/widths[keep];values=[v[keep]/widths[keep] for v in values]
    labels=[b[2] for b,k in zip(bins,keep) if k]
    if len(labels)==1 and not density: labels=['Inclusive']
    edges=np.arange(len(labels)+1)-.5;x=np.arange(len(labels))
    fig,(ax,ratio)=plt.subplots(2,1,figsize=(12,10),sharex=True,
        gridspec_kw={'height_ratios':[3.5,1],'hspace':.05})
    colors=['#E69F00','#009E73','#CC79A7','#D55E00','#56B4E9','#9467bd']
    hep.histplot(base,bins=edges,ax=ax,histtype='step',yerr=False,color='black',lw=2,label='SM')
    ratio.axhline(1,color='black',lw=1)
    for i,((g,t),v) in enumerate(zip(a.pairs,values)):
        c=colors[i%len(colors)]
        hep.histplot(v,bins=edges,ax=ax,histtype='step',yerr=False,color=c,lw=2,
                     label=f'P{i+1}: {a.coefficients[0]}={g:g}, {a.coefficients[1]}={t:g}')
        ratio.plot(x,np.divide(v,base,out=np.full_like(v,np.nan),where=base>0),'-o',color=c,ms=4)
    if a.ratio_ylim: ratio.set_ylim(*a.ratio_ylim)
    unit=a.unit if a.unit is not None else ('GeV' if obs in ('PTH','massZ2','massZ1','mass4l','mZ2','mZ1') else '')
    if effect=='shape_only': ylabel='Normalized density'+(f' [1/{unit}]' if density and unit else '') if density else 'Bin probability'
    elif a.quantity=='integrated-xsec': ylabel=f'Fiducial cross section [{a.xsec_unit}]'
    elif a.quantity=='xsec': ylabel=r'$\Delta\sigma / \Delta x$ ['+a.xsec_unit+('/'+unit if unit else '')+']'
    else: ylabel='Signal events'+(' / '+(unit or 'unit '+obs) if density else ' / bin')
    ax.set_ylabel(ylabel,fontsize=22);ratio.set_ylabel('EFT / SM',fontsize=18)
    known={'PTH':r'$p_T^H$ bin [GeV]','massZ2':r'$m_{Z_2}$ bin [GeV]','mass4l':r'$m_{4\ell}$ bin [GeV]'}
    ratio.set_xlabel(a.xlabel or ('Inclusive fiducial region' if len(labels)==1 and not density else known.get(obs,obs+' bin')),fontsize=21)
    ratio.set_xticks(x);ratio.set_xticklabels(labels,rotation=40 if len(labels)>4 else 0,ha='right' if len(labels)>4 else 'center',fontsize=13)
    ratio.set_xlim(-.5,len(labels)-.5);ratio.grid(axis='y',alpha=.2)
    ratio.tick_params(axis='x',which='minor',bottom=False,top=False)
    hep.cms.label(llabel=a.cms_label,rlabel=a.rlabel,data=True,ax=ax,fontsize=18)
    # Reserve the top half for the legend and inset, so neither hides spectra.
    highest=max(v.max() for v in [base]+values)
    if a.log_y:
        lowest=min(v.min() for v in [base]+values)
        if lowest<=0: raise ValueError('--log-y requires strictly positive plotted values')
        bottom=lowest*.5
        ax.set_yscale('log');ax.set_ylim(bottom,bottom*(highest/bottom)**2.15)
    else: ax.set_ylim(0,highest*2.35)
    ax.text(.025,.97,effect.replace('_',' + ' if effect=='rate_shape' else ' ').capitalize(),transform=ax.transAxes,va='top',fontsize=17)
    ax.legend(loc='upper left',bbox_to_anchor=(.01,.91),fontsize=12,frameon=False)
    if scans is not None:
        inset=ax.inset_axes(a.inset_bounds)
        level=-2*math.log(1-a.cl)  # chi-square with 2 degrees of freedom
        handles=[]
        for (xx,yy,z),name,col in zip(scans,a.scan_names,['blue','black']):
            if not z.min()<=level<=z.max(): warnings.warn(f'{name}: requested contour is outside sampled NLL range')
            else: inset.contour(xx,yy,z,levels=[level],colors=[col],linewidths=1.5,corner_mask=False)
            handles.append(Line2D([],[],color=col,label=name,lw=1.5))
        inset.plot(0,0,'*',color='gray',ms=8,zorder=5)
        for i,(g,t) in enumerate(a.pairs):
            inset.plot(g,t,'D',color=colors[i%len(colors)],mec='black',mew=.5,ms=6,zorder=6)
            inset.annotate(str(i+1),(g,t),xytext=(4,4),textcoords='offset points',fontsize=8)
        inset.set_xlabel(a.coefficient_labels[0],fontsize=13,labelpad=1)
        inset.set_ylabel(a.coefficient_labels[1],fontsize=13,labelpad=1)
        inset.minorticks_off()
        inset.tick_params(which='both',top=False,right=False,labelsize=9,length=3,width=.8)
        inset.legend(handles=handles,fontsize=8,title=f'{a.scan_label} {100*a.cl:g}% CL (2D)',title_fontsize=8,loc='upper right',ncol=2,frameon=False,columnspacing=.8,handlelength=1.8)
        if a.scan_xlim: inset.set_xlim(a.scan_xlim)
        if a.scan_ylim: inset.set_ylim(a.scan_ylim)
        else:
            ymin,ymax=inset.get_ylim();inset.set_ylim(ymin,ymax+.25*(ymax-ymin))
    fig.subplots_adjust(left=.14,right=.97,top=.90,bottom=.17)
    stem=Path(a.outdir)/(getattr(a,'plot_name','joint')+'_'+effect)
    for ext in ('png','pdf'): fig.savefig(str(stem)+'.'+ext,dpi=180)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--spectra',required=True,help='yields.csv from plot_pth_eft.py')
    p.add_argument('--scan-inclusive','--scan-1',dest='scan_inclusive',required=True)
    p.add_argument('--scan-pth','--scan-2',dest='scan_pth',required=True)
    p.add_argument('--scan-names',nargs=2,default=['Inclusive',r'$p_T^{4\ell}$'])
    p.add_argument('--coefficients',nargs=2,default=['chg','chgtil'])
    p.add_argument('--coefficient-labels',nargs=2,help='Optional display labels (mathtext allowed)')
    p.add_argument('--ratio-ylim',nargs=2,type=float)
    p.add_argument('--inset-bounds',nargs=4,type=float,default=[.57,.52,.40,.43])
    p.add_argument('--pairs',help='Values in --coefficients order; required for nondefault coefficients')
    p.add_argument('--quantity',choices=['counts','density','xsec','integrated-xsec'],default='density')
    p.add_argument('--effects',nargs='+',choices=['rate_shape','rate_only','shape_only'],default=['rate_shape','rate_only','shape_only'])
    p.add_argument('--cl',type=float,default=.95);p.add_argument('--scan-label',default='Expected')
    p.add_argument('--scan-xlim',nargs=2,type=float);p.add_argument('--scan-ylim',nargs=2,type=float)
    p.add_argument('--overflow-width',type=float);p.add_argument('--unit');p.add_argument('--xsec-unit',default='fb')
    p.add_argument('--xlabel');p.add_argument('--log-y',action='store_true')
    p.add_argument('--cms-label',default='Simulation Work in progress');p.add_argument('--rlabel',default='13.6 TeV')
    p.add_argument('--outdir',default='eft_scan_inset')
    a=p.parse_args()
    if a.coefficients[0]==a.coefficients[1]: p.error('Choose two distinct coefficients')
    a.coefficient_labels=a.coefficient_labels or a.coefficients
    if a.pairs is None:
        if a.coefficients!=['chg','chgtil']: p.error('Supply --pairs for the selected coefficients')
        a.pairs='-0.138,0;-0.07,0.058;-0.07,-0.058'
    if a.ratio_ylim and a.ratio_ylim[0]>=a.ratio_ylim[1]: p.error('Ratio limits must increase')
    if not 0<a.cl<1:p.error('--cl must be between 0 and 1')
    if a.overflow_width is not None and (not math.isfinite(a.overflow_width) or a.overflow_width<=0):p.error('Invalid overflow width')
    try:
        a.pairs=[tuple(float(v) for v in pair.split(',')) for pair in a.pairs.split(';')]
        if not a.pairs or any(len(v)!=2 or not all(math.isfinite(x) for x in v) for v in a.pairs):raise ValueError()
    except ValueError:p.error('--pairs must be finite coefficient1,coefficient2 pairs separated by semicolons')
    scans=[read_scan(a.scan_inclusive,a.coefficients),read_scan(a.scan_pth,a.coefficients)]
    obs,bins,sm,curves=read_spectra(a.spectra,a.pairs,a.quantity,a.coefficients)
    Path(a.outdir).mkdir(parents=True,exist_ok=True)
    for effect in a.effects:render(a,scans,obs,bins,sm,curves,effect)
    settings=vars(a).copy();settings['contour_2deltaNLL']=-2*math.log(1-a.cl)
    settings['total_scope']='All CSV bins including overflow'
    (Path(a.outdir)/'inset_settings.json').write_text(json.dumps(settings,indent=2))
    print('Saved plots to',a.outdir)

if __name__=='__main__':main()

