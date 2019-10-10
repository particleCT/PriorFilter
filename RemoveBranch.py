from ROOT import *
import glob,sys

f = TFile(sys.argv[1],"update")
t = f.Get("phase")
b = t.GetBranch("WET_prob")
t.GetListOfBranches().Remove(b);
l= t.GetLeaf("WET_prob");
t.GetListOfLeaves().Remove(l);
t.Write("",TObject.kOverwrite);
f.Close()


f = TFile(sys.argv[1],"update")
t = f.Get("phase")
c = t.GetBranch("Y_prob")
t.GetListOfBranches().Remove(c);
d= t.GetLeaf("Y_prob");
t.GetListOfLeaves().Remove(d);
t.Write("",TObject.kOverwrite);
f.Close()


f = TFile(sys.argv[1],"update")
t = f.Get("phase")
e = t.GetBranch("Y_prob")
t.GetListOfBranches().Remove(e);
f= t.GetLeaf("Y_prob");
t.GetListOfLeaves().Remove(f);
t.Write("",TObject.kOverwrite);
f.Close()



