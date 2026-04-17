#!/usr/bin/env python3

import os
import sys
import glob
import re
import json
import subprocess
import argparse


def run_multiqc(base_dir):
    print("\n[Dashboard] Running MultiQC locally...")
    out_dir = os.path.abspath(os.path.join(base_dir, "MultiQC_Report"))
    target_dir = os.path.abspath(base_dir)
    
    # Run the native MultiQC installed via Conda
    cmd = f"multiqc \"{target_dir}\" -o \"{out_dir}\""
    
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"✅ MultiQC Report generated at: {os.path.join(out_dir, 'multiqc_report.html')}")
    except Exception as e:
        print(f"⚠️ MultiQC run failed: {e}")

def parse_busco(summary_file):
    stats = {"Complete": 0, "Single": 0, "Duplicated": 0, "Fragmented": 0, "Missing": 0}
    try:
        with open(summary_file, 'r') as f:
            content = f.read()
            match = re.search(r"C:(.*?)%\[S:(.*?)%,D:(.*?)%\],F:(.*?)%,M:(.*?)%", content)
            if match:
                stats["Complete"] = float(match.group(1))
                stats["Single"] = float(match.group(2))
                stats["Duplicated"] = float(match.group(3))
                stats["Fragmented"] = float(match.group(4))
                stats["Missing"] = float(match.group(5))
    except Exception:
        pass
    return stats

def parse_quast(report_file, accession):
    stats = {"N50": 0, "Total_Length": 0, "GC": 0}
    try:
        with open(report_file, 'r') as f:
            lines = f.readlines()
            if not lines: return stats
            
            # Find which column belongs to this specific assembly
            header = lines[0].strip().split('\t')
            col_idx = 1
            for i, col in enumerate(header):
                if accession in col:
                    col_idx = i
                    break
                    
            for line in lines[1:]:
                parts = line.strip().split('\t')
                if len(parts) > col_idx:
                    if parts[0] == "N50":
                        stats["N50"] = int(parts[col_idx])
                    elif parts[0] == "Total length":
                        stats["Total_Length"] = int(parts[col_idx])
                    elif parts[0] == "GC (%)":
                        stats["GC"] = float(parts[col_idx])
    except Exception:
        pass
    return stats

def generate_custom_html(data, base_dir):
    print("[Dashboard] Generating Custom Interactive HTML Dashboard...")
    labels = list(data.keys())
    n50_data = [d["QUAST"]["N50"] for d in data.values()]
    busco_single = [d["BUSCO"]["Single"] for d in data.values()]
    busco_dup = [d["BUSCO"]["Duplicated"] for d in data.values()]
    busco_frag = [d["BUSCO"]["Fragmented"] for d in data.values()]
    busco_miss = [d["BUSCO"]["Missing"] for d in data.values()]

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AquaG Interactive Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #0f172a; color: #f8fafc; margin: 0; padding: 2rem; }}
        h1 {{ text-align: center; color: #38bdf8; margin-bottom: 2rem; }}
        .container {{ display: flex; flex-wrap: wrap; gap: 2rem; justify-content: center; }}
        .chart-wrapper {{ background: #1e293b; padding: 1.5rem; border-radius: 12px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); width: 45%; min-width: 400px; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 2rem; background: #1e293b; border-radius: 12px; overflow: hidden; }}
        th, td {{ padding: 1rem; text-align: left; border-bottom: 1px solid #334155; }}
        th {{ background-color: #38bdf8; color: #0f172a; }}
        tr:hover {{ background-color: #334155; }}
    </style>
</head>
<body>
    <h1>AquaG Genomic Assembly Dashboard</h1>
    <div class="container">
        <div class="chart-wrapper"><canvas id="buscoChart"></canvas></div>
        <div class="chart-wrapper"><canvas id="n50Chart"></canvas></div>
    </div>
    <table>
        <thead>
            <tr>
                <th>Accession</th>
                <th>Total Length (bp)</th>
                <th>N50 (bp)</th>
                <th>GC (%)</th>
                <th>BUSCO Complete (%)</th>
            </tr>
        </thead>
        <tbody>
            {''.join([f"<tr><td>{k}</td><td>{v['QUAST']['Total_Length']:,}</td><td>{v['QUAST']['N50']:,}</td><td>{v['QUAST']['GC']}</td><td>{v['BUSCO']['Complete']}</td></tr>" for k, v in data.items()])}
        </tbody>
    </table>
    <script>
        const labels = {json.dumps(labels)};
        new Chart(document.getElementById('buscoChart'), {{
            type: 'bar',
            data: {{
                labels: labels,
                datasets: [
                    {{ label: 'Single', data: {json.dumps(busco_single)}, backgroundColor: '#22c55e' }},
                    {{ label: 'Duplicated', data: {json.dumps(busco_dup)}, backgroundColor: '#3b82f6' }},
                    {{ label: 'Fragmented', data: {json.dumps(busco_frag)}, backgroundColor: '#eab308' }},
                    {{ label: 'Missing', data: {json.dumps(busco_miss)}, backgroundColor: '#ef4444' }}
                ]
            }},
            options: {{ responsive: true, plugins: {{ title: {{ display: true, text: 'BUSCO Completeness (%)', color: '#f8fafc' }}, legend: {{ labels: {{ color: '#f8fafc' }} }} }}, scales: {{ x: {{ stacked: true, ticks: {{ color: '#94a3b8' }} }}, y: {{ stacked: true, max: 100, ticks: {{ color: '#94a3b8' }} }} }} }}
        }});
        new Chart(document.getElementById('n50Chart'), {{
            type: 'bar',
            data: {{ labels: labels, datasets: [{{ label: 'N50 Length (bp)', data: {json.dumps(n50_data)}, backgroundColor: '#a855f7', borderRadius: 4 }}] }},
            options: {{ responsive: true, plugins: {{ title: {{ display: true, text: 'Assembly Contiguity (N50)', color: '#f8fafc' }}, legend: {{ labels: {{ color: '#f8fafc' }} }} }}, scales: {{ x: {{ ticks: {{ color: '#94a3b8' }} }}, y: {{ ticks: {{ color: '#94a3b8' }} }} }} }}
        }});
    </script>
</body>
</html>"""
    
    out_file = os.path.join(base_dir, "AquaG_Dashboard.html")
    with open(out_file, "w") as f:
        f.write(html_content)
    print(f"✅ Custom Dashboard generated at: {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate AquaG Dashboard")
    parser.add_argument("-d", "--dir", required=True, help="Base output directory")
    args = parser.parse_args()

    base_dir = os.path.abspath(args.dir)
    if not os.path.isdir(base_dir):
        print(f"❌ Error: Directory {base_dir} does not exist.")
        sys.exit(1)

    print(f"\n==========================================")
    print(f" Building AquaG Reports")
    print(f"==========================================")

    run_multiqc(base_dir)

    compiled_data = {}
    
    # 1. Use GLOB to recursively find files regardless of folder structure!
    busco_files = glob.glob(os.path.join(base_dir, "**", "short_summary.*.txt"), recursive=True)
    quast_reports = glob.glob(os.path.join(base_dir, "**", "quast_output", "report.tsv"), recursive=True)

    for b_file in busco_files:
        # Extract exact accession ID from the BUSCO filename itself
        acc = "Unknown"
        m = re.search(r"short_summary\.[^.]+\.[^.]+\.(.+)\.txt", os.path.basename(b_file))
        if m:
            acc = m.group(1)
            
        compiled_data[acc] = {
            "QUAST": {"N50": 0, "Total_Length": 0, "GC": 0}, 
            "BUSCO": parse_busco(b_file)
        }
        
        # 2. Extract matching QUAST stats
        if quast_reports:
            compiled_data[acc]["QUAST"] = parse_quast(quast_reports[0], acc)

    if compiled_data:
        generate_custom_html(compiled_data, base_dir)
    else:
        print("⚠️ No QUAST or BUSCO data found to visualize. Check directory paths.")

    print(f"==========================================\n")

if __name__ == "__main__":
    main()
