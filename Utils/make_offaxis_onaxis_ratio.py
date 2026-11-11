#!/usr/bin/env python3
"""
make_offaxis_onaxis_ratio.py

POT-normalizes the off-axis and on-axis event-rate histograms to a chosen
target exposure, sums them, plots the combined 2D event rates and 1D
projections (with a ratio band), and writes everything to a ROOT file.

POT normalization
------------------
Each histogram in the input files was already scaled by the production
framework to the POT declared in the sample config that made it (see
OFFAXIS_GENERATION_POT / ONAXIS_GENERATION_POT below). Critically, this
is NOT the same between the two files even for the "same" 0m position:
the off-axis config declares 0m at 5.5e20 POT, while the on-axis-only
config declares 0m (from the same underlying MC) at 1.1e21 POT. This
script rescales each histogram, TH2::Scale(target_POT / generation_POT),
before summing, so that:

  - the on-axis file's 0m histogram is scaled to represent --onaxis-pot
    POT (default 1.1e21) - with the correct generation POT this is a
    no-op (factor 1), since it was already produced at 1.1e21
  - the off-axis file's 0m histogram is scaled to represent
    --offaxis-0m-pot POT (default 5.5e20) - also a no-op by default
  - each of the other off-axis positions is scaled to represent an equal
    share of --offaxis-other-pot POT (default 5.5e20, split N ways)

By default the on-axis target (1.1e21) equals the off-axis targets summed
(5.5e20 + 5.5e20 = 1.1e21), so the two combined histograms represent the
same total exposure and the ratio is a fair equal-POT comparison of the
two run strategies.

Usage:
    python3 make_offaxis_onaxis_ratio.py \
        --offaxis eventrate_offaxis_template.root \
        --onaxis  eventrate_onaxis_template.root \
        --outdir  plots/

Requires PyROOT (i.e. `python3 -c "import ROOT"` must work).
"""

import argparse
import json
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kViridis)
ROOT.gStyle.SetNumberContours(255)  # smooth viridis gradient, not banded
ROOT.gStyle.SetTextFont(132)   # Times + LaTeX-style symbol font
ROOT.gStyle.SetLabelFont(132, "XYZ")
ROOT.gStyle.SetTitleFont(132, "XYZ")

# Generation POT for each sample - the "POT:" field from the sample config
# used to produce that histogram. Each output ROOT file's histograms are
# already scaled by the production framework to their config's declared
# POT, so this table must match the config that made THAT SPECIFIC file -
# it is NOT necessarily the same between the off-axis and on-axis files,
# even for the "same" 0m position, if the two configs declared different
# POT values for it (as is the case here: the off-axis config's 0m entry
# is declared at 5.5e20 POT, while the on-axis-only config's 0m entry -
# using the same underlying MC files - is declared at 1.1e21 POT).
#
# Override/extend with --offaxis-generation-pot-json / 
# --onaxis-generation-pot-json {"histogram_name": pot, ...} if your
# production used different values than these defaults.
OFFAXIS_GENERATION_POT = {
    "hOffAxis0m_numuCC_numode": 5.5e20,
    "hOffAxis1.75m_numuCC_numode": 3.24e19,
    "hOffAxis2m_numuCC_numode": 3.24e19,
    "hOffAxis4m_numuCC_numode": 3.24e19,
    "hOffAxis5.75m_numuCC_numode": 3.24e19,
    "hOffAxis8m_numuCC_numode": 3.24e19,
    "hOffAxis9.75m_numuCC_numode": 3.24e19,
    "hOffAxis12m_numuCC_numode": 3.24e19,
    "hOffAxis13.75m_numuCC_numode": 3.24e19,
    "hOffAxis16m_numuCC_numode": 3.24e19,
    "hOffAxis17.75m_numuCC_numode": 3.24e19,
    "hOffAxis20m_numuCC_numode": 3.24e19,
    "hOffAxis21.75m_numuCC_numode": 3.24e19,
    "hOffAxis24m_numuCC_numode": 3.24e19,
    "hOffAxis25.75m_numuCC_numode": 3.24e19,
    "hOffAxis26.25m_numuCC_numode": 3.24e19,
    "hOffAxis28.25m_numuCC_numode": 3.24e19,
    "hOffAxis28.5m_numuCC_numode": 3.24e19,
}

# From the on-axis-only sample config: 0m is declared at 1.1e21 POT, and
# every other position is declared at 0 POT (they're empty in this file).
ONAXIS_GENERATION_POT = {
    "hOffAxis0m_numuCC_numode": 1.1e21,
    "hOffAxis1.75m_numuCC_numode": 0.0,
    "hOffAxis2m_numuCC_numode": 0.0,
    "hOffAxis4m_numuCC_numode": 0.0,
    "hOffAxis5.75m_numuCC_numode": 0.0,
    "hOffAxis8m_numuCC_numode": 0.0,
    "hOffAxis9.75m_numuCC_numode": 0.0,
    "hOffAxis12m_numuCC_numode": 0.0,
    "hOffAxis13.75m_numuCC_numode": 0.0,
    "hOffAxis16m_numuCC_numode": 0.0,
    "hOffAxis17.75m_numuCC_numode": 0.0,
    "hOffAxis20m_numuCC_numode": 0.0,
    "hOffAxis21.75m_numuCC_numode": 0.0,
    "hOffAxis24m_numuCC_numode": 0.0,
    "hOffAxis25.75m_numuCC_numode": 0.0,
    "hOffAxis26.25m_numuCC_numode": 0.0,
    "hOffAxis28.25m_numuCC_numode": 0.0,
    "hOffAxis28.5m_numuCC_numode": 0.0,
}
ZERO_POSITION_KEY = "hOffAxis0m_numuCC_numode"


def load_generation_pot(defaults, json_path):
    """Start from the given built-in defaults dict, override/extend with
    any entries supplied in json_path."""
    merged = dict(defaults)
    if json_path is not None:
        with open(json_path) as f:
            merged.update(json.load(f))
    return merged


def get_sample_keys(root_file, key_prefix):
    """Sorted, de-duplicated list of histogram names starting with
    key_prefix (de-duplicated in case of multiple ROOT 'cycles')."""
    keys = sorted({k.GetName() for k in root_file.GetListOfKeys()
                   if k.GetName().startswith(key_prefix)})
    if not keys:
        raise RuntimeError(f"No histograms starting with '{key_prefix}' found "
                            f"in {root_file.GetName()}")
    return keys


def build_offaxis_target_pot(keys, pot_0m, pot_other_total):
    """0m -> pot_0m; every other key -> pot_other_total / n_other, split
    equally across whichever off-axis positions are present besides 0m."""
    others = [k for k in keys if k != ZERO_POSITION_KEY]
    if not others:
        raise RuntimeError("Only found the 0m off-axis histogram - can't "
                            "split the remaining POT across zero positions.")
    per_other = pot_other_total / len(others)
    target = {ZERO_POSITION_KEY: pot_0m}
    target.update({k: per_other for k in others})
    return target


def build_onaxis_target_pot(keys, pot_0m, generation_pot):
    """0m -> pot_0m; every other key -> its own generation POT (scale
    factor 1, or 0/0 -> 0 via the guard in scale_and_sum), which is
    harmless since those histograms should be empty in the on-axis-only
    file anyway."""
    target = {ZERO_POSITION_KEY: pot_0m}
    for k in keys:
        if k != ZERO_POSITION_KEY:
            target[k] = generation_pot.get(k, 0.0)
    return target


def scale_and_sum(root_file, keys, generation_pot, target_pot, label):
    """Scale each named TH2D in root_file from its generation POT to its
    target POT, then sum them all. generation_pot / target_pot: dicts
    keyed by histogram name. Returns the summed TH2D (detached from
    file)."""
    missing_gen = [k for k in keys if k not in generation_pot]
    if missing_gen:
        raise RuntimeError(f"No generation POT known for: {missing_gen} in "
                            f"'{label}'. Pass --offaxis-generation-pot-json or "
                            f"--onaxis-generation-pot-json to supply it.")

    combined = None
    print(f"[{label}] scaling {len(keys)} histograms:")
    for key in keys:
        h = root_file.Get(key)
        if not h:
            raise RuntimeError(f"Histogram '{key}' not found in {root_file.GetName()}")
        h = h.Clone(f"{key}_scaled")
        h.SetDirectory(0)
        gen = generation_pot[key]
        tgt = target_pot[key]
        factor = (tgt / gen) if gen != 0 else 0.0
        h.Scale(factor)
        print(f"    {key:36s} gen={gen:9.4g} POT  ->  target={tgt:9.4g} POT  "
              f"(x{factor:.4f})")
        if combined is None:
            combined = h
        else:
            combined.Add(h)

    total_target = sum(target_pot[k] for k in keys)
    print(f"    total target POT for {label} = {total_target:.4g}, "
          f"summed integral = {combined.Integral():.3f}")
    return combined


def normalize_1d_by_bin_width(hist, new_name):
    """Return a clone with every bin content/error divided by that bin's
    width, i.e. counts -> counts/width. Used only for the 1D projection
    overlay plots, since their binning is non-uniform."""
    h = hist.Clone(new_name)
    h.SetDirectory(0)
    for ix in range(1, h.GetNbinsX() + 1):
        w = h.GetXaxis().GetBinWidth(ix)
        if w <= 0:
            continue
        h.SetBinContent(ix, hist.GetBinContent(ix) / w)
        h.SetBinError(ix, hist.GetBinError(ix) / w)
    return h


def format_pot(pot):
    """1.1e21 -> '1.1#times10^{21} POT' (ROOT TLatex syntax)."""
    exponent = int(f"{pot:.1e}".split("e")[1])
    mantissa = pot / (10 ** exponent)
    return f"{mantissa:.1f}#times10^{{{exponent}}} POT"


def make_projection_ratio_plot(h_on_raw, h_total_raw, xlabel, ylabel_top, outpath,
                                col_on, col_total, col_ratio, pot_total=None):
    """Two-pad canvas: top = on-axis and on+off-axis overlaid (bin-width
    divided for display), bottom = ratio band (total/on-axis), computed
    from the RAW (undivided) counts so the bin width cancels exactly."""
    h_ratio = h_total_raw.Clone(h_total_raw.GetName() + "_ratio")
    h_ratio.Divide(h_on_raw)

    h_on = normalize_1d_by_bin_width(h_on_raw, h_on_raw.GetName() + "_density")
    h_total = normalize_1d_by_bin_width(h_total_raw, h_total_raw.GetName() + "_density")

    c = ROOT.TCanvas("c_ratio", "c_ratio", 900, 800)
    pad_top = ROOT.TPad("pad_top", "pad_top", 0, 0.3, 1, 1)
    pad_bot = ROOT.TPad("pad_bot", "pad_bot", 0, 0.0, 1, 0.3)
    pad_top.SetBottomMargin(0.02)
    pad_top.SetLeftMargin(0.14)
    pad_top.SetTopMargin(0.08)
    pad_bot.SetTopMargin(0.03)
    pad_bot.SetBottomMargin(0.35)
    pad_bot.SetLeftMargin(0.14)
    pad_top.Draw()
    pad_bot.Draw()

    pad_top.cd()
    h_total.SetLineColor(col_total)
    h_total.SetLineWidth(2)
    h_on.SetLineColor(col_on)
    h_on.SetLineWidth(2)
    h_total.SetTitle("")
    h_total.GetXaxis().SetLabelSize(0)
    h_total.GetXaxis().SetTitleSize(0)
    h_total.GetYaxis().SetTitle(ylabel_top)
    h_total.GetYaxis().SetTitleSize(0.05)
    h_total.GetYaxis().SetLabelSize(0.04)
    h_total.GetYaxis().SetTitleOffset(1.3)
    h_total.SetMaximum(1.3 * max(h_total.GetMaximum(), h_on.GetMaximum()))
    h_total.SetMinimum(0)
    h_total.Draw("HIST")
    h_on.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.045)
    leg.AddEntry(h_total, "On + Off-Axis (Combined)", "l")
    leg.AddEntry(h_on, "On-Axis Only", "l")
    leg.Draw()

    pot_latex = None
    if pot_total is not None:
        pot_latex = ROOT.TLatex()
        pot_latex.SetNDC()
        pot_latex.SetTextSize(0.04)
        pot_latex.DrawLatex(0.16, 0.83, format_pot(pot_total))

    pad_bot.cd()
    h_ratio.SetLineColor(col_ratio)
    h_ratio.SetLineWidth(2)
    h_ratio.SetMarkerColor(col_ratio)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.8)
    h_ratio.SetTitle("")
    h_ratio.GetXaxis().SetTitle(xlabel)
    h_ratio.GetYaxis().SetTitle("Combined / On-Axis")
    h_ratio.GetXaxis().SetTitleSize(0.11)
    h_ratio.GetXaxis().SetLabelSize(0.09)
    h_ratio.GetYaxis().SetTitleSize(0.11)
    h_ratio.GetYaxis().SetLabelSize(0.09)
    h_ratio.GetYaxis().SetTitleOffset(0.55)
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.Draw("PE")

    xmin = h_ratio.GetXaxis().GetXmin()
    xmax = h_ratio.GetXaxis().GetXmax()
    line = ROOT.TLine(xmin, 1, xmax, 1)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray + 2)
    line.Draw()

    c.SaveAs(outpath)
    print(f"  wrote {outpath}")
    return h_on, h_total, h_ratio


def make_1d_projections(h_off, h_on, outdir, fout, pot_total=None):
    """Project both combined (already POT-scaled) 2D histograms onto X
    (TrueNeutrinoEnergy) and Y (Enubias). h_off is the on+off-axis
    combined total; h_on is the on-axis-only histogram. Produces one
    overlay+ratio canvas per axis and writes the underlying histograms to
    fout."""
    col_total = ROOT.TColor.GetColor(68, 1, 84)      # viridis dark purple
    col_on = ROOT.TColor.GetColor(253, 231, 37)       # viridis yellow
    col_ratio = ROOT.TColor.GetColor(33, 145, 140)    # viridis mid-teal

    for axis, xlabel, ylabel_top, proj in (
            ("X", "E^{True}_{#nu} [GeV]", "Events / GeV", "ProjectionX"),
            ("Y", "E^{Bias}_{#nu}", "Events / bin width", "ProjectionY")):
        h_total_raw = getattr(h_off, proj)(f"h_total_combined_proj{axis}_raw")
        h_on_raw = getattr(h_on, proj)(f"h_onaxis_only_proj{axis}_raw")

        h_on_density, h_total_density, h_ratio = make_projection_ratio_plot(
            h_on_raw, h_total_raw, xlabel, ylabel_top,
            os.path.join(outdir, f"proj{axis}_onaxis_vs_combined_ratio.pdf"),
            col_on, col_total, col_ratio, pot_total=pot_total)

        fout.cd()
        h_total_raw.Write()
        h_on_raw.Write()
        h_total_density.Write()
        h_on_density.Write()
        h_ratio.Write()


def style_and_save_2d(hist, title, zlabel, outpath):
    c = ROOT.TCanvas("c", "c", 900, 700)
    c.SetRightMargin(0.16)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.12)

    hist.SetTitle(title)
    hist.GetXaxis().SetTitle("E^{True}_{#nu} [GeV]")
    hist.GetYaxis().SetTitle("E^{Bias}_{#nu}")
    hist.GetZaxis().SetTitle(zlabel)
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetZaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetZaxis().SetTitleOffset(1.2)

    hist.Draw("COLZ")
    c.SaveAs(outpath)
    print(f"  wrote {outpath}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--offaxis", required=True, help="eventrate_offaxis_template.root")
    p.add_argument("--onaxis", required=True, help="eventrate_onaxis_template.root")
    p.add_argument("--prefix", default="hOffAxis",
                    help="common name prefix of the histograms to sum (default: hOffAxis)")
    p.add_argument("--outdir", default=".", help="output directory")
    p.add_argument("--onaxis-pot", type=float, default=1.1e21,
                    help="target POT for the on-axis file's 0m histogram (default 1.1e21)")
    p.add_argument("--offaxis-0m-pot", type=float, default=5.5e20,
                    help="target POT for the off-axis file's 0m histogram (default 5.5e20)")
    p.add_argument("--offaxis-other-pot", type=float, default=5.5e20,
                    help="target POT to be split equally across the other "
                         "off-axis positions (default 5.5e20)")
    p.add_argument("--offaxis-generation-pot-json", default=None,
                    help="optional JSON file {histogram_name: generation_POT} "
                         "to override the built-in off-axis-file defaults")
    p.add_argument("--onaxis-generation-pot-json", default=None,
                    help="optional JSON file {histogram_name: generation_POT} "
                         "to override the built-in on-axis-file defaults")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    offaxis_generation_pot = load_generation_pot(OFFAXIS_GENERATION_POT,
                                                  args.offaxis_generation_pot_json)
    onaxis_generation_pot = load_generation_pot(ONAXIS_GENERATION_POT,
                                                 args.onaxis_generation_pot_json)

    f_off = ROOT.TFile.Open(args.offaxis)
    f_on = ROOT.TFile.Open(args.onaxis)
    if not f_off or f_off.IsZombie():
        raise RuntimeError(f"Could not open {args.offaxis}")
    if not f_on or f_on.IsZombie():
        raise RuntimeError(f"Could not open {args.onaxis}")

    off_keys = get_sample_keys(f_off, args.prefix)
    on_keys = get_sample_keys(f_on, args.prefix)

    offaxis_target = build_offaxis_target_pot(off_keys, args.offaxis_0m_pot,
                                               args.offaxis_other_pot)
    onaxis_target = build_onaxis_target_pot(on_keys, args.onaxis_pot, onaxis_generation_pot)

    h_off = scale_and_sum(f_off, off_keys, offaxis_generation_pot, offaxis_target,
                           label="off-axis file (on+off combined)")
    h_on = scale_and_sum(f_on, on_keys, onaxis_generation_pot, onaxis_target,
                          label="on-axis file (on-axis only)")

    h_off.SetName("h_offaxis_combined")
    h_on.SetName("h_onaxis_combined")

    if (h_off.GetNbinsX() != h_on.GetNbinsX()
            or h_off.GetNbinsY() != h_on.GetNbinsY()):
        raise RuntimeError("Off-axis and on-axis histograms have different "
                            "binning - cannot take a bin-by-bin ratio.")

    total_offaxis_pot = args.offaxis_0m_pot + args.offaxis_other_pot
    print(f"[POT normalization] on-axis target total  = {args.onaxis_pot:.4g} POT")
    print(f"[POT normalization] off-axis target total = {total_offaxis_pot:.4g} POT "
          f"({args.offaxis_0m_pot:.4g} at 0m + {args.offaxis_other_pot:.4g} split "
          f"across {len(off_keys) - 1} other positions)")
    if abs(total_offaxis_pot - args.onaxis_pot) > 1e-6 * max(total_offaxis_pot, args.onaxis_pot):
        print("  WARNING: on-axis and off-axis target totals do not match - "
              "the ratio plot will partly reflect a POT difference, not "
              "purely spectral shape.")

    h_ratio = h_off.Clone("h_ratio_offaxis_over_onaxis")
    h_ratio.Divide(h_on)

    style_and_save_2d(h_off, "Combined Off-Axis Event Rate", "Events",
                       os.path.join(args.outdir, "offaxis_combined.pdf"))
    style_and_save_2d(h_on, "Combined On-Axis Event Rate", "Events",
                       os.path.join(args.outdir, "onaxis_combined.pdf"))
    style_and_save_2d(h_ratio, "Off-Axis / On-Axis Ratio", "Ratio",
                       os.path.join(args.outdir, "ratio_offaxis_over_onaxis.pdf"))

    out_root = os.path.join(args.outdir, "combined_and_ratio.root")
    fout = ROOT.TFile.Open(out_root, "RECREATE")
    h_off.Write()
    h_on.Write()
    h_ratio.Write()

    make_1d_projections(h_off, h_on, args.outdir, fout, pot_total=args.onaxis_pot)

    fout.Close()
    print(f"  wrote {out_root}")

    f_off.Close()
    f_on.Close()


if __name__ == "__main__":
    main()
