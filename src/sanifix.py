#! /bin/env/python
from rdkit.Chem import AllChem
import pandas as pd
from rdkit import Chem

##### Sanifix 4 by James Davidson ######
""" sanifix4.py

  Contribution from James Davidson
"""


def _FragIndicesToMol(oMol, indices):
    em = Chem.rdchem.EditableMol(Chem.rdchem.Mol())

    newIndices = {}
    for i, idx in enumerate(indices):
        em.AddAtom(oMol.GetAtomWithIdx(idx))
        newIndices[idx] = i

    for i, idx in enumerate(indices):
        at = oMol.GetAtomWithIdx(idx)
        for bond in at.GetBonds():
            if bond.GetBeginAtomIdx() == idx:
                oidx = bond.GetEndAtomIdx()
            else:
                oidx = bond.GetBeginAtomIdx()
            # make sure every bond only gets added once:
            if oidx < idx:
                continue
            em.AddBond(newIndices[idx], newIndices[oidx], bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    Chem.rdmolops.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap = newIndices
    return res


def _recursivelyModifyNs(mol, matches, indices=None):
    if indices is None:
        indices = []
    res = None
    while len(matches) and res is None:
        tIndices = indices[:]
        nextIdx = matches.pop(0)
        tIndices.append(nextIdx)
        nm = Chem.rdchem.Mol(mol)
        nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
        nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(0)
        cp = Chem.rdchem.Mol(nm)
        try:
            Chem.rdmolops.SanitizeMol(cp)
            res, indices = _recursivelyModifyNs(nm, matches, indices=tIndices)
        except ValueError:
            indices = tIndices
            res = cp
    return res, indices


def AdjustAromaticNs(m,nitrogenPattern='[n&D2&H0;r5,r6]'):
    """
       default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
       to fix: O=c1ccncc1
    """
    Chem.GetSymmSSSR(m)
    m.UpdatePropertyCache(False)

    # break non-ring bonds linking rings:
    em = Chem.EditableMol(m)
    linkers = m.GetSubstructMatches(Chem.MolFromSmarts('[r]!@[r]'))
    plsFix=set()
    for a,b in linkers:
        em.RemoveBond(a,b)
        plsFix.add(a)
        plsFix.add(b)
    nm = em.GetMol()
    for at in plsFix:
        at=nm.GetAtomWithIdx(at)
        if at.GetIsAromatic() and at.GetAtomicNum()==7:
            at.SetNumExplicitHs(1)
            at.SetNoImplicit(True)

    # build molecules from the fragments:
    fragLists = Chem.GetMolFrags(nm)
    frags = [_FragIndicesToMol(nm,x) for x in fragLists]

    # loop through the fragments in turn and try to aromatize them:
    ok=True
    for i,frag in enumerate(frags):
        cp = Chem.Mol(frag)
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
            lres,indices=_recursivelyModifyNs(frag,matches)
            if not lres:
                #print 'frag %d failed (%s)'%(i,str(fragLists[i]))
                ok=False
                break
            else:
                revMap={}
                for k,v in frag._idxMap.items():
                    revMap[v]=k
                for idx in indices:
                    oatom = m.GetAtomWithIdx(revMap[idx])
                    oatom.SetNoImplicit(True)
                    oatom.SetNumExplicitHs(1)
    if not ok:
        return None
    return m

    
if __name__ == "__main__":
    mol = Chem.MolFromSmiles('CC(C)(O)C1:C:C:C(C2:N:C(=O):C3:C:C(F):C:C:C:3:N:2):C:C:1', sanitize=False)
    m = AdjustAromaticNs(mol)
    print(Chem.MolToSmiles(m))
