"""Verbatim CRISPOR MIT/Hsu hit score (crispor.py, maximilianh/crisporWebsite master, retrieved 2026-09-24)."""
# aka Matrix "M"
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]

def calcHitScore(string1,string2):
    " see 'Scores of single hits' on http://crispr.mit.edu/about "
    # The Patrick Hsu weighting scheme
    # S. aureus requires 21bp long guides. We fudge by using only last 20bp
    matrixStart = 0
    maxDist = 19

    assert(string1[0].isupper())
    assert(len(string1)==len(string2))
    #for nmCas9 and a few others with longer guides, we limit ourselves to 20bp
    if len(string1)>20:
        string1 = string1[-20:]
        string2 = string2[-20:]
    # for 19bp guides, we fudge a little, but first pos has no weight anyways
    elif len(string1)==19:
        string1 = "A"+string1
        string2 = "A"+string2
    # for shorter guides, I'm not sure if this score makes sense anymore, we force things
    elif len(string1)<19:
        matrixStart = 20-len(string1)
        maxDist = len(string1)-1

    assert(len(string1)==len(string2))

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(matrixStart, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            if lastMmPos!=None:
                dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos
    # 2nd part of the score
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((maxDist-avgDist)/float(maxDist)) * 4 + 1)
    # 3rd part of the score
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

def calcMitGuideScore(hitSum):
    """ Sguide defined on http://crispr.mit.edu/about
    Input is the sum of all off-target hit scores. Returns the specificity of the guide.
    """
    score = 100 / (100+hitSum)
    score = int(round(score*100))
    return score


