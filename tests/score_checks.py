#!/usr/bin/env python3
"""Audit unified Rust output against count inputs and independent legacy Perl output."""
import csv
import math
from pathlib import Path


def rows(path):
    with Path(path).open() as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames
        if not fields or len(fields)!=len(set(fields)):raise ValueError('Missing/duplicate score columns')
        data=list(reader)
        if any(None in r or None in r.values() for r in data):raise ValueError('Malformed scored row')
        return fields,data


def require(ok,message):
    if not ok:raise ValueError(message)


def near(a,b,tolerance=1e-12):
    require(math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=tolerance),f'Score value differs: {a} != {b}')


def check_mrna(unified,oracle):
    fields,data=rows(unified);legacy,other=rows(oracle)
    require('ScoreMode' in fields and fields[:fields.index('ScoreMode')]==legacy,'mRNA canonical prefix differs')
    key=lambda r:tuple(r[k] for k in legacy[:3]);a={key(r):r for r in data};b={key(r):r for r in other}
    require(len(a)==len(data) and len(b)==len(other) and a.keys()==b.keys(),'mRNA candidate keys differ')
    for k,r in a.items():
        require(r['ScoreMode']=='mrna' and r['CoordinateSystem']=='genome_1based' and
                r['EvaluationStatus']=='evaluated' and r['Candidate']=='true' and
                r['FC_criteria']==r['Score_criteria']=='meet','Invalid mRNA candidate status')
        for i,field in enumerate(legacy):
            if i<7:require(r[field]==b[k][field],'mRNA annotation differs')
            else:near(r[field],b[k][field])
    return dict(status='passed',candidate_sites=len(data),scope='Independent Perl candidate keys, metadata and numeric score prefix')


def check_rrna(unified,oracle,inputs,reference):
    _,data=rows(unified);_,other=rows(oracle)
    require(len(data)==len(other),'rRNA row counts differ')
    names=[Path(p).parent.name for p in inputs];counts=[];bases=[]
    for path in inputs:
        with Path(path).open() as f:values=[r.rstrip().split('\t') for r in f]
        counts.append([int(r[3]) for r in values]);bases.append([r[1] for r in values])
    require(all(b==bases[0] for b in bases),'rRNA reference bases differ')
    totals=[sum(c) for c in counts];valid=all(t>0 for t in totals);selected=[]
    require(len(data)==len(bases[0])-1,'rRNA lacks complete candidate positions')
    for p,(r,old) in enumerate(zip(data,other),1):
        require(int(r['Position'])==int(old['#Position'])==p and r['Base']==old['Base']==bases[0][p-1] and
                r['Chromosome']==reference and r['Direction']=='+' and r['ScoreMode']=='rrna' and
                r['CoordinateSystem']=='rrna_reference_1based','rRNA candidate coordinate differs')
        boundary=p<=20 or p>len(counts[0])-20
        status='normalization_unavailable' if not valid else 'boundary_excluded' if boundary else 'evaluated'
        require(r['EvaluationStatus']==status,'rRNA evaluation status differs')
        require(int(r[names[0]+'_score_numerator'])==counts[0][p] and
                int(r[names[0]+'_score_denominator'])==max(counts[0][p-1],1),'NC ratio counts differ')
        for i,name in enumerate(names):
            signal=name+('_count(NC)' if i==0 else '_count')
            require(int(r[signal])==counts[i][p] and int(r[name+'_adjacent_count'])==counts[i][p-1] and
                    int(old[name+'_count'])==counts[i][p-1],'Signal/adjacent count semantics differ')
            if i:
                require(int(r[name+'_NC_TOTAL'])==totals[0] and int(r[name+'_TR_TOTAL'])==totals[i], 'Normalization totals differ')
                require(int(r[name+'_score_denominator'])==max(counts[i][p-1],1),'TR ratio denominator differs')
                if totals[0]>0 and totals[i]>0:near(r[name+'_normalization_factor'],totals[0]/totals[i])
                else:require(r[name+'_normalization_factor']=='NA','Unavailable factor must be NA')
                if status=='evaluated':
                    require(int(r[name+'_count_corr'])==math.ceil(counts[i][p]*totals[0]/totals[i]),'Corrected count differs')
                    near(r[name+'_score'],old[name+'_score'],.00501)
                else:require(r[name+'_count_corr']==r[name+'_score']=='NA','Excluded values must be NA')
        called=old['FC_criteria']==old['Score_criteria']=='meet'
        require((r['Candidate']=='true')==called,'rRNA candidate decision differs')
        if status=='evaluated':
            for new,legacy in [('TotalScore','TotalScore'),('TotalFoldChange','FoldChange')]:near(r[new],old[legacy],.00501)
            require(r['FC_criteria']==old['FC_criteria'] and r['Score_criteria']==old['Score_criteria'],'rRNA criteria differ')
            near(float(r['TotalScore'])/(len(inputs)-1),r['AveScore'])
            near(float(r['TotalFoldChange'])/(len(inputs)-1),r['AveFC'])
            corrected=sum(int(r[s+'_count_corr']) for s in names[1:])
            require(int(r['TotalCorrectedCount'])==corrected,'Corrected sum differs')
            near(r['TotalFoldChange'],corrected/max(counts[0][p],1))
        else:
            require(r['FC_criteria']==r['Score_criteria']=='not_evaluated' and not called,'Excluded row called')
            require(r['TotalScore']==r['AveScore']==r['TotalFoldChange']==r['AveFC']=='NA','Unavailable aggregates must be NA')
        if called:selected.append(p)
    return dict(status='passed',candidate_sites=len(selected),candidates=selected,rows=len(data),
                scope='Independent Perl calls and rounded scores; raw signal/adjacent counts, normalization and unified output semantics')
