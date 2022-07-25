import ROOT as r
from optparse import OptionParser,OptionGroup

def setlabel(h):
    stations = ["4/2","4/1","3/2","3/1","2/2","2/1","1/3","1/2","1/1"]
    for i,station in enumerate(stations):
        label_minus = "ME-"+station
        label_plus = "ME+"+station
        h.GetXaxis().SetBinLabel(i+1,label_minus)
        h.GetXaxis().SetBinLabel(18-i,label_plus)
    return

def main(options):

    #f = r.TFile.Open("plots_351405.root")
    f = r.TFile.Open(options.inputFile)
    shower_types = ["alct",'clct','lct']
    shower_thres = ['nom','tight']
    r.gROOT.SetBatch(1)
 
    c1=r.TCanvas("c1","c1",800,600)
    r.gStyle.SetOptStat(0)
    h_count = f.Get("cscTriggerPrimitivesAnalyzer/Counters")
    h_count.Draw("HIST TEXT0")
    h_count.GetXaxis().SetRangeUser(0,15)
    c1.SaveAs(options.odir+"counters.pdf")
    for thres in shower_thres:
        c1=r.TCanvas("c1","c1",2400,1200)
        r.gStyle.SetOptStat(0)
        c1.Divide(3,2)
        for i,shower_type in enumerate(shower_types):
            hname_denom =  shower_type+"_cscshower_data_"+thres+"_summary_denom"
            hname_num =  shower_type+"_cscshower_data_"+thres+"_summary_num"
            h_denom    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_denom)
            h_num    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_num)
            if h_denom.Integral()!=0:
                print(" %s , %s : %s/%s ~ %.2f" %(shower_type,thres,h_num.Integral(),h_denom.Integral(),h_num.Integral()/h_denom.Integral()))
            else:
                print(" %s , %s : %s/%s " %(shower_type,thres,h_num.Integral(),h_denom.Integral()))
            c1.cd(i+1)
            h_denom.Draw("colz")
            c1.cd(i+4)
            h_num.Draw("colz")
        c1.SaveAs(options.odir+"summary_%s.pdf"%thres)
    c1=r.TCanvas("c1","c1",2400,1200)
    r.gStyle.SetOptStat(1)
    c1.Divide(3,2)
    for i,shower_type in enumerate(shower_types):
        hname_denom =  shower_type+"_cscshower_data_nom_summary_denom"
        h_1d    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_denom).ProjectionY("_py",0,-1,"o")
        h_1d.GetXaxis().SetRangeUser(0,18)
        setlabel(h_1d)
        c1.cd(i+1)
        h_1d.Draw()
        c1.cd(i+4)
        hname_denom =  shower_type+"_cscshower_data_tight_summary_denom"
        h_1d_tight    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_denom).ProjectionY("_py",0,-1,"o")
        h_1d_tight.GetXaxis().SetRangeUser(0,18)
        setlabel(h_1d_tight)
        h_1d_tight.Draw()
    c1.SaveAs(options.odir+"Occupancy_summary.pdf")

    c1=r.TCanvas("c1","c1",2400,1200)
    r.gStyle.SetOptStat(0)
    c1.Divide(3,2)
    for i,shower_type in enumerate(shower_types):
        hname_denom =  shower_type+"_nom_Bx_v_station"
        h_2d    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_denom)
        c1.cd(i+1)
        h_2d.Draw("colz")
        c1.cd(i+4)
        hname_denom =  shower_type+"_tight_Bx_v_station"
        h_2d_tight    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_denom)
        h_2d_tight.Draw("colz")
    c1.SaveAs(options.odir+"Bx_station_summary.pdf")

            
            
    for thres in shower_thres:
        c1=r.TCanvas("c1","c1",2400,1200)
        r.gStyle.SetOptStat(0)
        c1.Divide(3,2)
        for i,shower_type in enumerate(shower_types):
            hname_denom =  shower_type+"_cscshower_emul_"+thres+"_summary_denom"
            hname_num =  shower_type+"_cscshower_emul_"+thres+"_summary_num"
            h_denom    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_denom)
            h_num    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_num)
            if h_denom.Integral()!=0:
                print("Unmatched emulation:  %s , %s : %s/%s ~ %.2f" %(shower_type,thres,h_num.Integral(),h_denom.Integral(),h_num.Integral()/h_denom.Integral()))
            else:
                print("Unmatched emulation:  %s , %s : %s/%s " %(shower_type,thres,h_num.Integral(),h_denom.Integral(),))
            c1.cd(i+1)
            h_denom.Draw("colz")
            c1.cd(i+4)
            h_num.Draw("colz")
        c1.SaveAs(options.odir+"summary_%s_emul.pdf"%thres)

    for thres in shower_thres:
        for i,shower_type in enumerate(shower_types):
            c1=r.TCanvas("c1","c1",800,600)
            r.gStyle.SetOptStat(0)
            hname_denom =  shower_type+"_cscshower_data_"+thres+"_summary_denom"
            hname_num =  shower_type+"_cscshower_data_"+thres+"_summary_num"
            h_denom    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_denom)
            h_num    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_num)
            h_num.Divide(h_denom)
            h_num.Draw("colz")
            c1.SaveAs(options.odir+"Eff_%s_%s.pdf"%(shower_type,thres))
    
    for thres in shower_thres:
        for i,shower_type in enumerate(shower_types):
            c1=r.TCanvas("c1","c1",800,600)
            r.gStyle.SetOptStat(1)
            if thres=="nom": thres = "norm"
            hname_data =  "%sBx_data_"%shower_type+thres
            hname_emul =  "%sBx_emul_"%shower_type+thres
            h_data    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_data)
            h_emul    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_emul)
            
            if h_data.Integral()!=0: h_data.Scale(1./h_data.Integral())
            if h_emul.Integral()!=0: h_emul.Scale(1./h_emul.Integral())
            print(hname_data, h_data.GetMean())
            print(hname_emul, h_emul.GetMean())
            h_emul.SetLineColor(r.kRed)
            h_data.Draw("HIST")
            h_emul.Draw("same HIST")
            c1.SaveAs(options.odir+"Bx_%s_%s_eff.pdf"%(shower_type,thres))

    for i,shower_type in enumerate(shower_types):
        c1=r.TCanvas("c1","c1",800,600)
        r.gStyle.SetOptStat(0)
        hname_Bx =  "%s_Bx_data_v_emul"%shower_type
        hname_thres =  "%s_thres_data_v_emul"%shower_type
        h_Bx    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_Bx)
        h_thres  = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_thres)
        h_Bx.Draw("colz TEXT0")
        c1.SaveAs(options.odir+"%s_Bx_data_v_emul.pdf"%(shower_type))
        c1.Clear()
        h_thres.Draw("colz TEXT0")
        c1.SaveAs(options.odir+"%s_thres_data_v_emul.pdf"%(shower_type))

    c1=r.TCanvas("c1","c1",2400,1200)
    r.gStyle.SetOptStat(0)
    c1.Divide(3,2)
    for i,shower_type in enumerate(shower_types):
        hname_Bx =  "%s_Bx_data_v_emul"%shower_type
        hname_thres =  "%s_thres_data_v_emul"%shower_type
        h_Bx    = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_Bx)
        h_thres  = f.Get("cscTriggerPrimitivesAnalyzer/"+hname_thres)
        c1.cd(i+1)
        h_Bx.Draw("colz TEXT0")
        c1.cd(i+4)
        h_thres.Draw("colz TEXT0")
    c1.SaveAs(options.odir+"summary_Bx_thres.pdf")



if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-i', '--inputFile', dest='inputFile', default='./plots.root', help='plots output', metavar='inputFile')
    parser.add_option('-o', '--odir', dest='odir', default='./', help='output dir', metavar='odir')
    parser.add_option('-r', '--run', dest='run', default=353018, help='run number', metavar='run')

    (options, args) = parser.parse_args()

    main(options)
