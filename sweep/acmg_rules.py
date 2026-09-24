def acmg2015(counts):
    """counts: dict with keys PVS,PS,PM,PP,BA,BS,BP (after strength modifiers)."""
    vs,s,m,p=counts['PVS'],counts['PS'],counts['PM'],counts['PP']
    ba,bs,bp=counts['BA'],counts['BS'],counts['BP']
    P = (vs>=1 and (s>=1 or m>=2 or (m==1 and p>=1) or p>=2)) or s>=2 or \
        (s==1 and (m>=3 or (m==2 and p>=2) or (m==1 and p>=4))) or vs>=2
    LP = (vs>=1 and m==1) or (s==1 and 1<=m<=2) or (s==1 and p>=2) or m>=3 or (m==2 and p>=2) or (m==1 and p>=4)
    B = ba>=1 or bs>=2
    LB = (bs==1 and bp>=1) or bp>=2
    path = 'Pathogenic' if P else 'Likely pathogenic' if LP else None
    ben = 'Benign' if B else 'Likely benign' if LB else None
    if path and ben: return 'Uncertain significance'
    return path or ben or 'Uncertain significance'
