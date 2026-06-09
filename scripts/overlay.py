#!/usr/bin/env python3
"""

python3 scripts/overlay.py \
    -a plots/performance/energy_res_no_bkg.root -b plots/performance/energy_res.root \
    --label-a 'no-bkg' --label-b 'bkg' \
    -o plots/performance/energy_res_overlay.root

"""

import argparse
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Overlay matching canvases from two ROOT files."
    )
    parser.add_argument("-a", "--input-a", required=True, help="First ROOT file")
    parser.add_argument("-b", "--input-b", required=True, help="Second ROOT file")
    parser.add_argument("-o", "--output", required=True, help="Output ROOT file")
    parser.add_argument("--label-a", default="A", help="Legend label for first file")
    parser.add_argument("--label-b", default="B", help="Legend label for second file")
    return parser.parse_args()


def collect_canvases(directory, prefix=""):
    canvases = {}
    for key in directory.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        path = f"{prefix}/{name}" if prefix else name
        if obj.InheritsFrom("TDirectory"):
            canvases.update(collect_canvases(obj, path))
        elif obj.InheritsFrom("TCanvas"):
            canvases[path] = obj
    return canvases


def ensure_dir(root_file, directory):
    if not directory:
        return root_file
    current = root_file
    for part in directory.split("/"):
        next_dir = current.GetDirectory(part)
        if not next_dir:
            next_dir = current.mkdir(part)
        current = next_dir
    return current


def plottables(canvas):
    out = []
    for obj in canvas.GetListOfPrimitives():
        if obj.InheritsFrom("TH1") or obj.InheritsFrom("TGraph"):
            out.append(obj.Clone())
    return out


def set_style(obj, color):
    if hasattr(obj, "SetLineColor"):
        obj.SetLineColor(color)
    if hasattr(obj, "SetMarkerColor"):
        obj.SetMarkerColor(color)
    if hasattr(obj, "SetMarkerStyle"):
        obj.SetMarkerStyle(20)
    if hasattr(obj, "SetLineWidth"):
        obj.SetLineWidth(2)


def draw_overlay(path, canvas_a, canvas_b, out_file, label_a, label_b):
    objs_a = plottables(canvas_a)
    objs_b = plottables(canvas_b)
    if not objs_a and not objs_b:
        return

    out_canvas = ROOT.TCanvas(
        Path(path).name,
        canvas_a.GetTitle() or canvas_b.GetTitle(),
        1000,
        800,
    )
    out_canvas.SetLogx(canvas_a.GetLogx() or canvas_b.GetLogx())
    out_canvas.SetLogy(canvas_a.GetLogy() or canvas_b.GetLogy())
    out_canvas.SetLogz(canvas_a.GetLogz() or canvas_b.GetLogz())

    legend = ROOT.TLegend(0.68, 0.75, 0.88, 0.88)
    legend.SetBorderSize(1)

    first = True
    for obj in objs_a:
        set_style(obj, ROOT.kBlue + 1)
        draw_opt = "hist" if obj.InheritsFrom("TH1") else "AP"
        if not first:
            draw_opt += " same"
        obj.Draw(draw_opt)
        if first:
            first = False
        legend.AddEntry(obj, label_a, "l")
        break

    for obj in objs_b:
        set_style(obj, ROOT.kRed + 1)
        draw_opt = "hist" if obj.InheritsFrom("TH1") else "P"
        if not first:
            draw_opt += " same"
        obj.Draw(draw_opt)
        if first:
            first = False
        legend.AddEntry(obj, label_b, "l")
        break

    legend.Draw()

    directory = str(Path(path).parent)
    out_dir = ensure_dir(out_file, "" if directory == "." else directory)
    out_dir.cd()
    out_canvas.Write()
    out_file.cd()


def main():
    args = parse_args()

    input_a = ROOT.TFile.Open(args.input_a, "READ")
    input_b = ROOT.TFile.Open(args.input_b, "READ")
    if not input_a or input_a.IsZombie():
        raise SystemExit(f"Failed to open {args.input_a}")
    if not input_b or input_b.IsZombie():
        raise SystemExit(f"Failed to open {args.input_b}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output = ROOT.TFile(args.output, "RECREATE")
    if not output or output.IsZombie():
        raise SystemExit(f"Failed to open {args.output}")

    canvases_a = collect_canvases(input_a)
    canvases_b = collect_canvases(input_b)
    common_paths = sorted(set(canvases_a) & set(canvases_b))

    for path in common_paths:
        draw_overlay(
            path,
            canvases_a[path],
            canvases_b[path],
            output,
            args.label_a,
            args.label_b,
        )

    output.Close()
    input_a.Close()
    input_b.Close()


if __name__ == "__main__":
    main()
