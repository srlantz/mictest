from ROOT import TTree

t = TTree
t.ReadFile('/home/matevz/mictest/mkFit/zzzz-host.rtt')

for e in t:
   print e.sumWaste
