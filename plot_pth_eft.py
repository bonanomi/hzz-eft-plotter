#!/usr/bin/env python3
"""Evaluate any two Wilson coefficients and plot rate/shape effects.

Keep plot_eft_scan_inset.py beside this file (shared plotting code).
Requires PyROOT/Combine, numpy, matplotlib, mplhep. No uproot needed here.
--coefficients chwb chwbtil chooses variables, CSV columns, and plot labels.
Other ModelConfig POIs are set to zero at each benchmark; nuisance parameters
remain at the loaded workspace/snapshot values. No fit is performed.

--quantity: counts, density (events/width, default), xsec (cross section/width),
or integrated-xsec. Xsec modes evaluate xs_hzz_smH_<observable>_<bin> times
full_scaling_hzz_<bin>, once per generator bin. Units default to fb; no unit
conversion. --group-by truth groups SELECTED events, not unfolded cross sections.

Arbitrary names: --bin-map is a JSON object mapping yield function names to
{observable,bin,low,high}. --xsec-map is a JSON list of the same bin metadata
plus factors: [workspace variable/function names to multiply]. Entries are summed
within each bin. Cross-section mappings must exclude luminosity/efficiency.
Open-ended bins (GT... or upper edge 10000) are included in totals but omitted
from density plots unless --overflow-width supplies a display width.

Outputs yields.csv with INTEGRATED contents, settings.json, and three effect
subdirectories with SM, single-coefficient and joint PNG/PDF plots.
For default chg/chgtil the old --chg-values / --chgtil-values options still work.
"""
import argparse,csv,json,math,re
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from plot_eft_scan_inset import render

TOKEN=r'(?:m?\d+(?:p\d+)?_m?\d+(?:p\d+)?|GT\d+(?:p\d+)?)'
BIN_END=re.compile(r'^(?P<observable>.+)_(?P<bin>'+TOKEN+r')$')
PATTERN=re.compile(r'^n_exp_final_bin(?P<channel>.+)_proc_(?P<process>.+)$')
def number(s):return float(s.replace('m','-').replace('p','.'))
def edges(token):
    if token.startswith('GT'):return number(token[2:]),math.inf
    lo,hi=map(number,token.split('_'));return lo,math.inf if hi==10000 else hi
def infer_bin(stem):
    m=BIN_END.fullmatch(stem)
    if not m:return stem,0.,1.,'Inclusive'
    lo,hi=edges(m['bin']);return m['observable'],lo,hi,(f'>{lo:g}' if math.isinf(hi) else f'{lo:g}–{hi:g}')
def parse_yield_name(name,group_by='reco',include_extra=False,signal_pattern=r'^smH_'):
    m=PATTERN.fullmatch(name)
    if not m:return None
    proc=m['process']
    if not re.search(signal_pattern,proc) and not (include_extra and proc in ('OutsideAcceptance','nonResH')):return None
    stem=m['channel'].split('_cat',1)[0] if group_by=='reco' else proc
    if stem.startswith('hzz_'):stem=stem[4:]
    if group_by=='truth' and stem.startswith('smH_'):stem=stem[4:]
    return infer_bin(stem)
def mapped_bin(e):
    lo,hi=float(e.get('low',0)),float(e.get('high',1))
    if not math.isfinite(lo) or math.isnan(hi) or hi<=lo:raise ValueError('Invalid mapped bin edges')
    return str(e['observable']),lo,hi,str(e['bin'])
def collection_items(coll):
    it=coll.createIterator()
    while True:
        obj=it.Next()
        if not obj:break
        yield obj

def get_records(w,a):
    rec=[]
    if a.quantity in ('xsec','integrated-xsec'):
        if a.xsec_map:entries=json.loads(Path(a.xsec_map).read_text())
        else:
            entries=[]
            for obj in list(collection_items(w.allVars()))+list(collection_items(w.allFunctions())):
                n=obj.GetName()
                if not n.startswith('xs_hzz_smH_'):continue
                stem=n[len('xs_hzz_smH_'):];key=infer_bin(stem)
                if a.observable and key[0]!=a.observable:continue
                m=BIN_END.fullmatch(stem)
                if not m:raise ValueError('Use --xsec-map for '+n)
                entries.append(dict(observable=key[0],low=key[1],high=key[2],bin=key[3],factors=[n,'full_scaling_hzz_'+m['bin']]))
        for e in entries:
            key=mapped_bin(e)
            if a.observable and key[0]!=a.observable:continue
            factors=[]
            for name in e['factors']:
                obj=w.function(name) or w.var(name)
                if not obj:raise ValueError('Missing cross-section factor '+name)
                factors.append(obj)
            if not factors:raise ValueError('Empty cross-section factors')
            rec.append((key,tuple(factors)))
    else:
        mapping=json.loads(Path(a.bin_map).read_text()) if a.bin_map else None
        found=set()
        for obj in collection_items(w.allFunctions()):
            n=obj.GetName()
            if mapping is not None:
                key=mapped_bin(mapping[n]) if n in mapping else None
                if n in mapping:found.add(n)
            else:key=parse_yield_name(n,a.group_by,a.include_extra,a.signal_pattern)
            if key and (not a.observable or key[0]==a.observable):rec.append((key,(obj,)))
        if mapping is not None and set(mapping)-found:raise ValueError('Missing mapped yields: '+str(set(mapping)-found))
    return rec

def floats(s):
    v=[float(x) for x in s.split(',')]
    if not all(math.isfinite(x) for x in v):raise ValueError('Nonfinite benchmark')
    return v

def evaluate(w,pois,records,bins,coefficients,point):
    for var in pois:
        if not var.getMin()<=0<=var.getMax():raise ValueError(var.GetName()+' does not allow zero')
        var.setVal(0)
    for name,value in zip(coefficients,point):
        var=w.var(name)
        if not var.getMin()<=value<=var.getMax():raise ValueError(f'{name}={value} outside workspace range')
        var.setVal(value)
    result=np.zeros(len(bins));index={b:i for i,b in enumerate(bins)}
    for key,factors in records:
        value=math.prod(float(o.getVal()) for o in factors)
        if not math.isfinite(value) or value<0:raise ValueError('Invalid prediction: '+' * '.join(o.GetName() for o in factors))
        result[index[key]]+=value
    if result.sum()<=0:raise ValueError('Nonpositive total prediction')
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('root_file');p.add_argument('--workspace',default='w');p.add_argument('--snapshot')
    p.add_argument('--library',action='append',default=[])
    p.add_argument('--coefficients',nargs=2,default=['chg','chgtil'])
    p.add_argument('--values-1','--chg-values',dest='values_1',default='-0.01,0.01')
    p.add_argument('--values-2','--chgtil-values',dest='values_2',default='-0.01,0.01')
    p.add_argument('--pairs',default='0.01,0.01;-0.01,0.01',help='Illustrative benchmarks in --coefficients order')
    p.add_argument('--observable');p.add_argument('--group-by',choices=['reco','truth'],default='reco')
    p.add_argument('--signal-pattern',default=r'^smH_');p.add_argument('--include-extra',action='store_true')
    p.add_argument('--bin-map');p.add_argument('--xsec-map')
    p.add_argument('--quantity',choices=['counts','density','xsec','integrated-xsec'],default='density')
    p.add_argument('--counts',action='store_true');p.add_argument('--overflow-width',type=float)
    p.add_argument('--xsec-unit',default='fb');p.add_argument('--unit');p.add_argument('--xlabel')
    p.add_argument('--ratio-ylim',nargs=2,type=float);p.add_argument('--log-y',action='store_true')
    p.add_argument('--cms-label',default='Simulation Work in progress');p.add_argument('--rlabel',default='')
    p.add_argument('--outdir',default='pth_eft_plots')
    a=p.parse_args()
    if len(set(a.coefficients))!=2:p.error('Choose two distinct coefficients')
    if a.counts:
        if a.quantity in ('xsec','integrated-xsec'):p.error('--counts conflicts with cross-section mode')
        a.quantity='counts'
    if a.include_extra and a.group_by=='truth':p.error('--include-extra requires reco grouping')
    if a.quantity in ('xsec','integrated-xsec') and (a.bin_map or a.include_extra):p.error('Use --xsec-map for cross sections')
    if a.overflow_width is not None and (not math.isfinite(a.overflow_width) or a.overflow_width<=0):p.error('Invalid overflow width')
    if a.ratio_ylim and a.ratio_ylim[0]>=a.ratio_ylim[1]:p.error('Ratio limits must increase')
    pairs=[tuple(floats(s)) for s in a.pairs.split(';')]
    if any(len(v)!=2 for v in pairs):p.error('--pairs needs two values per benchmark')
    if len(set(pairs))!=len(pairs):p.error('Duplicate --pairs benchmarks')
    groups={a.coefficients[0]:[(v,0.) for v in floats(a.values_1)],a.coefficients[1]:[(0.,v) for v in floats(a.values_2)],'joint':pairs}
    import ROOT
    ROOT.gROOT.SetBatch(True)
    for lib in ['libHiggsAnalysisCombinedLimit']+a.library:
        if ROOT.gSystem.Load(lib)<0:raise RuntimeError('Cannot load '+lib)
    f=ROOT.TFile.Open(a.root_file,'READ')
    if not f or f.IsZombie():raise RuntimeError('Cannot open ROOT file')
    w=f.Get(a.workspace)
    if not w:raise RuntimeError('Missing workspace')
    if a.snapshot:
        if not w.getSnapshot(a.snapshot):raise RuntimeError('Missing snapshot')
        w.loadSnapshot(a.snapshot)
    mc=w.obj('ModelConfig')
    if not mc:raise RuntimeError('Missing ModelConfig')
    pois=list(collection_items(mc.GetParametersOfInterest()))
    if any(not v.InheritsFrom('RooRealVar') for v in pois):raise ValueError('POIs must be RooRealVars')
    if not set(a.coefficients)<={v.GetName() for v in pois}:raise ValueError('Selected coefficients must be ModelConfig POIs')
    records=get_records(w,a);bins=sorted({b for b,_ in records})
    if not bins:raise ValueError('No matching signal objects')
    if len({b[0] for b in bins})!=1:raise ValueError('Select one --observable')
    saved=[v.getVal() for v in pois]
    try:
        sm=evaluate(w,pois,records,bins,a.coefficients,(0.,0.))
        results={k:[evaluate(w,pois,records,bins,a.coefficients,pt) for pt in pts] for k,pts in groups.items()}
        if not np.allclose(sm,evaluate(w,pois,records,bins,a.coefficients,(0.,0.)),rtol=1e-10,atol=1e-12):raise RuntimeError('SM reset closure failed')
    finally:
        for v,val in zip(pois,saved):v.setVal(val)
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    xs=a.quantity in ('xsec','integrated-xsec');column='cross_section' if xs else 'events'
    with (out/'yields.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['group','benchmark',*a.coefficients,'observable','bin_label','bin_low','bin_high',column,'EFT_over_SM','bin_width','differential_cross_section' if xs else 'event_density'])
        rows=[('SM',(0.,0.),sm)]+[(k,pt,v) for k,pts in groups.items() for pt,v in zip(pts,results[k])]
        for k,pt,values in rows:
            for (obs,lo,hi,label),v,ref in zip(bins,values,sm):
                width=hi-lo if math.isfinite(hi) else a.overflow_width
                writer.writerow([k,str(pt),*pt,obs,label,lo,hi,v,v/ref if ref>0 else '',width or '',v/width if width else ''])
    settings=vars(a).copy();settings.update(zeroed_POIs=[v.GetName() for v in pois],source_factors=[[v.GetName() for v in factors] for _,factors in records],total_scope='All bins including overflow',nuisance_policy='Loaded workspace/snapshot values; no fit')
    (out/'settings.json').write_text(json.dumps(settings,indent=2))
    obs=bins[0][0];plotbins=[b[1:] for b in bins]
    for effect in ('rate_shape','rate_only','shape_only'):
        dest=out/effect;dest.mkdir(exist_ok=True)
        for k,pts,values in [('SM',[],[])]+[(k,groups[k],results[k]) for k in groups]:
            opts=SimpleNamespace(**{**vars(a),'outdir':str(dest),'pairs':pts,'plot_name':k})
            render(opts,None,obs,plotbins,sm,values,effect)
    print('Saved figures, yields.csv and settings.json to',out)
if __name__=='__main__':main()

