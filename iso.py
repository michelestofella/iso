# =============================================================================
#
# iso.py
#
# Evaluates the theoretical isotopic envelope given the intrinsic rates and 
# protection factors of a specific peptide.
#
# The fully protonated envelope must be evaluated from 
#     prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msisotope
# and saved as N.prospector.ucsf.edu, N being the number of the peptide considered.
#
# Arguments to be parsed:
#   --pep: number of peptide to be considered
#   --ass: file containing assignments (must be .ass file)
#   --kint: file containing intrinsic rates (must be .kint file)
#   --pfact: file containinf protection factors (must be .pfact file)
#   --times: file containing times at which isotopic envelope has to be evaluated (must be .times file)
#
# =============================================================================

''' Libraries '''
import pandas as pd
from iso_fun import *
import os
import argparse

''' Parse input arguments '''
parser = argparse.ArgumentParser()

parser.add_argument("--ass")
parser.add_argument("--pep")
parser.add_argument("--kint")
parser.add_argument("--pfact")
parser.add_argument("--times")

config = {}
opts = parser.parse_args()

if opts.ass:
    ass_file = str(opts.ass)+'.ass'
if opts.pep:
    pep = int(opts.pep)
if opts.kint:
    kint_file = str(opts.kint)+'.kint'
if opts.pfact:
    lnP_file = str(opts.pfact)+'.pfact'
if opts.times:
    times_file = str(opts.times)+'.times'
        
pi0_file  = str(pep)+'.prospector.ucsf.edu' 

times = []
with open(times_file,'r') as f:
    for line in f.readlines():
        times.append(float(line))
f.close()
    
''' Select residues involving the selected peptide '''
ass  = pd.read_csv(ass_file,header=None,delim_whitespace=True)
start_res = ass.iloc[int(pep)-1][1]
end_res   = ass.iloc[int(pep)-1][2]

''' Upload kint and lnP values '''
kint = list(pd.read_csv(kint_file,header=None,delim_whitespace=True)[1])[start_res:end_res]
lnP  = list(pd.read_csv(lnP_file,header=None,delim_whitespace=True)[1])[start_res:end_res]

''' Upload fully protonated isotopic envelope '''
pi0  = pd.read_csv(pi0_file,header=None,skiprows=1,sep='\t')
mass = list(pi0[1])
fr0  = list(pi0[2])
while len(mass) != len(kint):
    mass.append(mass[-1]+1.0)
    fr0.append(0)

''' Calculate isotopic envelopes at different times '''
for i in range(len(times)):
    f1 = centered_isotopic_envelope(times[i],kint,lnP,fr0)
    with open(lnP_file+'.'+str(pep)+'.'+str(i)+'.isot','w+') as f:
        for j in range(len(f1)):
            f.write('%5.5f ' % mass[j])
            if j == len(f1)-1:
                f.write('%5.5f' % f1[j])
            else:
                f.write('%5.5f\n' % f1[j])
    f.close()
   
# %%
