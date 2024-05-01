import json
import os
import sys
import ROOT
from types import SimpleNamespace


def get_attr(data, value):
    if hasattr(data, value):
        return getattr(data, value)
    else:
        print("Error: No attribute named '" + value + "' found")
        return None


def process_functions(data, obj, atr='functions'):
    if obj is not None:
        functions = get_attr(data, atr)
        if functions and isinstance(functions, list):
            for df in functions:
                name = get_attr(df, 'name')
                params = get_attr(df, 'params')
                if name and params:
                    try:
                        method = getattr(obj, df.name)
                        method(*tuple(df.params))
                    except Exception as e:
                        print("Error: Failed to call the function " + df.name + " with params" + df.params)
                        print(e)
    else:
        print("Error: object None!")


class Plot:
    def __init__(self, name, cut, sample, region, color):
        self.histogram = None
        self.name = name
        self.cut = cut
        self.sample = sample
        self.color = color
        self.region = region
        self.full_name = name + "_" + cut + "_" + sample + "_" + region
        self.shortened_name = name + "_" + cut + "_" + region

    def get_histogram(self, files):
        if self.histogram is not None:
            return self.histogram

        for f in files:
            h = f.Get(self.full_name)
            if isinstance(h, ROOT.TH1):
                self.histogram = h
                return self.histogram
        return None

    def functions(self, data):
        process_functions(data, self.histogram)

    def process(self, data):
        if self.histogram is not None:
            self.histogram.SetLineColor(self.color)
            self.histogram.SetLineColor(self.color)
            self.histogram.SetFillColor(self.color)
            self.histogram.SetMarkerColor(self.color)
            self.functions(data)
            return self.histogram
        else:
            print("Error: histogram " + self.full_name + " not initialised!")
            return None


def sort_plots(plots):
    plot_tree = {}

    for plot in plots:
        if plot.shortened_name in plot_tree:
            plot_tree[plot.shortened_name].append(plot)
        else:
            h = [plot]
            plot_tree[plot.shortened_name] = h
    return plot_tree


def plotter():
    output_root_file = "stack.root"
    output_directory = "stack"

    default_colors = [40, 50, 60, 70, 80, 90, 100]
    default_index = 0

    files = []
    plots = []
    cuts = []
    names = []
    regions = []
    samples = []
    colors = {}
    labels = {}

    with open('config.json') as jsonfile:
        data = json.load(jsonfile, object_hook=lambda d: SimpleNamespace(**d))

    color_data = get_attr(data, 'colors')
    file_data = get_attr(data, 'files')
    if file_data:
        for d in file_data:
            try:
                f = ROOT.TFile.Open(d, "OPEN")
                try:
                    for a in f.GetListOfKeys():
                        h = a.ReadObj()
                        s = h.GetName().split("_")
                        if len(s) == 4:
                            name = s[0]
                            cut = s[1]
                            sample = s[2]
                            region = s[3]

                            if name not in names:
                                names.append(name)
                            if cut not in cuts:
                                cuts.append(cut)
                            if region not in regions:
                                regions.append(region)
                            if sample not in samples:
                                samples.append(sample)
                                if color_data:
                                    if hasattr(color_data, sample):
                                        colors[sample] = getattr(color_data, sample)[0]
                                        labels[sample] = getattr(color_data, sample)[1]

                            if sample not in colors:
                                colors[sample] = default_colors[default_index]
                                default_index += 1
                                labels[sample] = sample

                            plot = Plot(name, cut, sample, region, colors[sample])
                            plots.append(plot)
                    files.append(f)
                except:
                    print("Error: Loading plots in root file")
                    continue
            except:
                print("Error: Could not open root file: " + d.file)
    else:
        sys.exit(0)

    if len(files) == 0:
        print("Error: No input files found!")
        sys.exit(0)

    if len(plots) == 0:
        print("Error: No input plots found!")
        sys.exit(0)

    if not os.path.isdir(_directory):
        os.mkdir(output_directory)
    os.chdir(output_directory)

    oF = ROOT.TFile.Open(output_root_file, "RECREATE")

    plot_tree = sort_plots(plots)

    data_plots = get_attr(data, 'plots')
    if data_plots:
        for d_p in data_plots:
            sps = samples if not hasattr(d_p, "samples") else d_p.samples
            cts = cuts if not hasattr(d_p, "cuts") else d_p.cuts
            rgs = regions if not hasattr(d_p, "cuts") else d_p.regions
            nms = names if not hasattr(d_p, "names") else d_p.names
            legend_data = get_attr(d_p, "legend")

            for ct in cts:
                for r in rgs:
                    for n_i in range(len(nms)):
                        n = nms[n_i]
                        prf = n + "_" + ct + "_" + r
                        if prf in plot_tree:
                            stack = ROOT.THStack()
                            if legend_data:
                                if hasattr(legend_data, "location"):
                                    leg = ROOT.TLegend(*tuple(legend_data.location))
                                else:
                                    leg = ROOT.TLegend(*tuple(20, 20, 50, 50))
                                if hasattr(legend_data, "functions"):
                                    process_functions(legend_data, leg)
                            hists = []
                            for s in sps:
                                for p in plot_tree[prf]:
                                    if p.sample == s:
                                        p.get_histogram(files)
                                        h = p.process(d_p)
                                        if h is not None:
                                            stack.Add(h)
                                            a = h.Clone()
                                            a.Scale(1. / h.Integral())
                                            a.SetFillColor(0)
                                            hists.append(a)
                                            if legend_data:
                                                leg.AddEntry(h, labels[s], "L" if not hasattr(legend_data, 'entries') else legend_data.entries)
                            plot_name = get_attr(d_p, "plot_name")
                            c = ROOT.TCanvas("" if not plot_name else plot_name + "_" + p.shortened_name)
                            stack.Draw("hist")

                            subplot = get_attr(d_p, 'subplot')
                            if subplot:
                                subplot_sample = get_attr(subplot, 'sample')
                                if subplot_sample and subplot_sample in samples:
                                    for p in plot_tree[prf]:
                                        if p.sample == subplot_sample:
                                            p.get_histogram(files)
                                            h = p.process(d_p)
                                            if h is not None:
                                                h.Draw("same" + ("" if not hasattr(subplot, "draw") else subplot.draw))
                                                if legend_data:
                                                    leg.AddEntry(h, labels[subplot_sample], "L" if not hasattr(legend_data,
                                                                                          'entries') else legend_data.entries)

                            if legend_data:
                                leg.Draw("" if not hasattr(legend_data, 'draw') else legend_data.draw)

                            process_functions(d_p, stack, atr='stack_functions')

                            if hasattr(d_p, "title_x"):
                                stack.GetXaxis().SetTitle(d_p.title_x[n_i])

                            if hasattr(d_p, "title_y"):
                                stack.GetYaxis().SetTitle(d_p.title_y[n_i])

                            label = getattr(d_p, "label")
                            if label:
                                if hasattr(label, "x") and hasattr(label, "y") and hasattr(label, "text"):
                                    latex = ROOT.TLatex(label.x, label.y, label.text)
                                    latex.SetTextSize(0.03)
                                    latex.Draw("same")

                            c.Print(".png")
                            c.Print(".pdf")
                            oF.cd()
                            stack.Write()

                            c = ROOT.TCanvas("" if not plot_name else plot_name + "_" + p.shortened_name + "_normalised")

                            maxValue = 0
                            maxHistoSample = ""
                            for h in hists:
                                maxHistValue = h.GetMaximum()
                                print(maxValue, maxHistValue)
                                if maxValue < maxHistValue:
                                    maxValue = maxHistValue
                                    maxHistoSample = h

                            maxHistoSample.GetYaxis().SetTitle("Normalised")
                            if hasattr(d_p, "title_x"):
                                maxHistoSample.GetXaxis().SetTitle(d_p.title_x[n_i])

                            maxHistoSample.SetTitle("")
                            maxHistoSample.Draw("hist")
                            maxHistoSample.SetStats(0)
                            leg.Draw("")
                            if label:
                                if hasattr(label, "x") and hasattr(label, "y") and hasattr(label, "text"):
                                    latex = ROOT.TLatex(label.x, 0.042, label.text)
                                    latex.SetTextSize(0.03)
                                    latex.Draw("same")
                            for h in hists:
                                if h != maxHistoSample: h.Draw("same hist")
                            c.Print(".png")
                            c.Print(".pdf")
                            oF.cd()


if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True)
    plotter()
