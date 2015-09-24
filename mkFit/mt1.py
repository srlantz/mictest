from ROOT import *

# gStyle.SetPalette(53)

t = TTree()

t.ReadFile("/home/matevz/mictest/mkFit/zzzz-host.rtt")
VW = 8

cntSumSize  = [0] * 11
cntSumWaste = [0] * 11
cntTotOps   = [0] * 11

for e in t:
    if e.N_proc != VW : continue

    cntSumSize [e.layer] += e.sumSize
    cntSumWaste[e.layer] += e.sumWaste
    cntTotOps  [e.layer] += VW * e.maxSize


for l in range(3, 10):
    print ("%2d %6d %6d %6d" %
           ( l, cntSumSize[l], cntSumWaste[l], cntTotOps[l]) )
print
    
for l in range(3, 10):

    print ("%2d %6d %6d %6d %.5f" %
           ( l, cntSumSize[l], cntSumWaste[l], cntTotOps[l],
             float(cntSumWaste[l]) / cntTotOps[l]) )

    cntSumSize [0] += cntSumSize[l]
    cntSumWaste[0] += cntSumWaste[l]
    cntTotOps  [0] += cntTotOps[l]

print "\nTOTALS:"
print ("   %6d %6d %6d %.5f" %
       ( cntSumSize[0], cntSumWaste[0], cntTotOps[0],
         float(cntSumWaste[0]) / cntTotOps[0]) )

    
h = TH2I("maxSize","maxSize:layer;layer;maxSize",
         7, 2.5, 9.5, 26, -0.5, 25.5)
h.SetStats(0)
t.Draw("maxSize:layer>>maxSize", "N_proc==%d" % VW, "goff")

c = TCanvas("MaxSize", "MaxSize", 1024, 480)
c.Divide(2,1)
c.cd(1)
h.Draw("colz")
c.cd(2)
h.Draw("lego2")
gPad.SetPhi(125)
gPad.SaveAs("maxSize.png")

c = TCanvas("N_proc", "N_proc", 512, 480)
t.Draw("N_proc:layer", "", "lego2")
gPad.SetPhi(61)
gPad.SetTheta(10.5)
gPad.SaveAs("N_proc.png")
