import sys
import os
import re
import argparse
import shutil
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Attempt to import the user-provided library
try:
    import AIPSTBOUTLib
except ImportError:
    print("[ERROR] AIPSTBOUTLib.py not found. Please ensure it is in the same directory.")
    sys.exit(1)

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def parse_log_file(log_path):
    """
    Parses the datacheck log file to extract Runlevels, JMFIT stats, and Table references.
    """
    runlevels = []
    current_level = None
    
    # regex for starting a runlevel
    rx_runlevel = re.compile(r"Runlevel\s+(\d+):\s+(.*)")
    
    # regex for table loading (TBIN)
    # Note: AIPS logs often split filenames across lines. We'll capture the start and look ahead.
    rx_tbin_start = re.compile(r"TBIN\s+\d+:\s+ZTXOP2:\s+using translated file name =")
    
    # regex for JMFIT
    rx_jmfit_source = re.compile(r"JMFIT1:\s+Source=\s+(\S+)")
    rx_jmfit_rms = re.compile(r"JMFIT1:\s+True rms.*=\s+(\S+)")
    rx_jmfit_peak = re.compile(r"JMFIT1:\s+Peak intensity\s*=\s+(\S+)\s+\+/-\s+(\S+)")
    rx_jmfit_flux = re.compile(r"JMFIT1:\s+Integral intensity\s*=\s+(\S+)\s+\+/-\s+(\S+)")
    rx_jmfit_res = re.compile(r"JMFIT1:\s+Component appears unresolved")

    with open(log_path, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # 1. Detect Runlevel
        m_run = rx_runlevel.match(line)
        if m_run:
            # Save previous level if exists
            if current_level:
                runlevels.append(current_level)
            
            current_level = {
                'id': m_run.group(1),
                'desc': m_run.group(2),
                'tables': set(),
                'jmfit': []
            }
            i += 1
            continue

        if current_level is None:
            i += 1
            continue

        # 2. Detect Tables (TBIN)
        if rx_tbin_start.search(line):
            # The filename follows in subsequent lines starting with "TBIN ... ZTXOP2:"
            full_path = ""
            j = i + 1
            while j < len(lines):
                next_line = lines[j].strip()
                if "ZTXOP2:" in next_line:
                    # Extract part after ZTXOP2:
                    part = next_line.split("ZTXOP2:")[-1].strip()
                    full_path += part
                    # A heuristic to stop: if it looks like a file extension or next line isn't ZTXOP2
                    if part.endswith(('.sn', '.bp', '.tab')): 
                        break
                else:
                    break
                j += 1
            
            # Clean up path
            if "/tables/" in full_path:
                fname = full_path.split("/tables/")[-1].strip()
                current_level['tables'].add(fname)
            elif full_path.strip():
                current_level['tables'].add(os.path.basename(full_path.strip()))
            
            i = j # Advance main loop
            continue

        # 3. Detect JMFIT
        m_src = rx_jmfit_source.search(line)
        if m_src:
            # Start a new source entry
            current_level['jmfit'].append({
                'Source': m_src.group(1),
                'RMS': 'N/A', 'Peak': 'N/A', 'Flux': 'N/A', 'Note': ''
            })
        
        if current_level['jmfit']:
            curr_fit = current_level['jmfit'][-1]
            
            m_rms = rx_jmfit_rms.search(line)
            if m_rms: curr_fit['RMS'] = m_rms.group(1)
            
            m_peak = rx_jmfit_peak.search(line)
            if m_peak: curr_fit['Peak'] = f"{m_peak.group(1)} +/- {m_peak.group(2)}"
            
            m_flux = rx_jmfit_flux.search(line)
            if m_flux: curr_fit['Flux'] = f"{m_flux.group(1)} +/- {m_flux.group(2)}"
            
            if rx_jmfit_res.search(line):
                curr_fit['Note'] = "Unresolved"

        i += 1

    if current_level:
        runlevels.append(current_level)

    return runlevels

def determine_plot_params(table_name):
    """
    Decides what to plot based on the table filename.
    Returns list of dicts args for AIPSTBOUTTable.plot()
    """
    # Defaults
    plots_to_make = []
    
    if 'bpass.bp' in table_name:
        plots_to_make.append({'plot_type': 'bandpass', 'xaxis': 'frequency'})
    elif 'fring' in table_name:
        plots_to_make.append({'plot_type': 'delay', 'xaxis': 'time'})
        plots_to_make.append({'plot_type': 'phase', 'xaxis': 'time'})
    elif 'accor' in table_name:
        plots_to_make.append({'plot_type': 'amplitude', 'xaxis': 'time'})
    elif 'apcal' in table_name:
        plots_to_make.append({'plot_type': 'amplitude', 'xaxis': 'time'})
        plots_to_make.append({'plot_type': 'phase', 'xaxis': 'time'})
    elif 'calib' in table_name or 'phs' in table_name:
        plots_to_make.append({'plot_type': 'phase', 'xaxis': 'time'})
    else:
        # Fallback
        plots_to_make.append({'plot_type': 'phase', 'xaxis': 'time'})
        plots_to_make.append({'plot_type': 'amplitude', 'xaxis': 'time'})
        
    return plots_to_make

def process_and_plot(runlevels, tables_dir, output_plot_dir):
    """
    Iterates through runlevels, finds tables, generates plots.
    Returns a dictionary mapping table_name -> list of generated image files (relative path)
    """
    table_map = {}
    ensure_dir(output_plot_dir)
    
    # Get unique tables across all levels to avoid re-plotting
    unique_tables = set()
    for rl in runlevels:
        unique_tables.update(rl['tables'])
    
    print(f"[INFO] Found {len(unique_tables)} unique tables to process.")

    for tbl in unique_tables:
        file_path = os.path.join(tables_dir, tbl)
        if not os.path.exists(file_path):
            print(f"[WARN] Table file not found: {file_path}")
            continue
            
        print(f"[INFO] Processing {tbl}...")
        
        try:
            # Use the library
            aips_obj = AIPSTBOUTLib.AIPSTBOUTTable()
            aips_obj.read(file_path)
            
            # Identify antennas present in data to link them later
            # (The library plots separate files per antenna)
            # We will grab the 'ANTENNA' column if available to know count
            ant_col = aips_obj._get_column_fuzzy(['ANTENNA', 'ANTENNA NO', 'STATION'])
            if ant_col is not None:
                present_ants = np.unique(ant_col[~np.isnan(ant_col)]).astype(int)
            else:
                present_ants = []

            plot_configs = determine_plot_params(tbl)
            generated_files = []

            for config in plot_configs:
                # Construct a prefix that puts files in our output dir
                # AIPSTBOUTLib appends _type_vs_xaxis_antX.png
                safe_tbl_name = tbl.replace('.', '_')
                prefix = os.path.join(output_plot_dir, safe_tbl_name)
                
                # Suppress library print output slightly
                aips_obj.plot(
                    plot_type=config['plot_type'], 
                    xaxis=config['xaxis'], 
                    out_prefix=prefix
                )
                
                # Predict filenames generated by the library
                for ant in present_ants:
                    fname = f"{prefix}_{config['plot_type']}_vs_{config['xaxis']}_ant{ant}.png"
                    if os.path.exists(fname):
                        generated_files.append(fname)
            
            table_map[tbl] = generated_files

        except Exception as e:
            print(f"[ERROR] Failed to plot {tbl}: {e}")
            import traceback
            traceback.print_exc()

    return table_map

def generate_html(runlevels, table_images, log_name, html_path):
    """
    Generates the HTML report.
    """
    
    css = """
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; color: #333; }
        h1 { border-bottom: 2px solid #007bff; padding-bottom: 10px; color: #007bff; }
        h2 { background-color: #f8f9fa; padding: 10px; border-left: 5px solid #28a745; margin-top: 30px; }
        .timestamp { float: right; color: #777; font-size: 0.9em; }
        .meta { margin-bottom: 30px; background: #eee; padding: 10px; border-radius: 5px; }
        
        table { width: 100%; border-collapse: collapse; margin: 15px 0; font-size: 0.9em; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #007bff; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        
        .plot-group { display: flex; flex-wrap: wrap; gap: 10px; justify-content: center; }
        .plot-card { border: 1px solid #ddd; padding: 5px; border-radius: 4px; width: 45%; min-width: 400px; text-align: center; }
        .plot-card img { width: 100%; height: auto; }
        .plot-caption { font-size: 0.8em; color: #555; margin-top: 5px; }
        
        .no-data { font-style: italic; color: #999; }
    </style>
    """

    html = f"""<!DOCTYPE html>
    <html>
    <head>
        <title>Pipeline Diagnostic: {log_name}</title>
        {css}
    </head>
    <body>
        <div class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
        <h1>Pipeline Diagnostic Report</h1>
        <div class="meta">
            <strong>Log File:</strong> {log_name}<br>
            <strong>Working Directory:</strong> {os.getcwd()}
        </div>
    """

    for rl in runlevels:
        has_tables = len(rl['tables']) > 0
        has_jmfit = len(rl['jmfit']) > 0
        
        # Skip empty runlevels to keep report clean
        if not (has_tables or has_jmfit):
            continue

        html += f"""
        <div class="runlevel">
            <h2>Runlevel {rl['id']}: {rl['desc']}</h2>
        """

        # JMFIT Table
        if has_jmfit:
            html += "<h3>Image Statistics (JMFIT)</h3>"
            html += """
            <table>
                <thead>
                    <tr><th>Source</th><th>RMS</th><th>Peak Intensity</th><th>Total Flux</th><th>Notes</th></tr>
                </thead>
                <tbody>
            """
            for fit in rl['jmfit']:
                html += f"""
                <tr>
                    <td>{fit['Source']}</td>
                    <td>{fit['RMS']}</td>
                    <td>{fit['Peak']}</td>
                    <td>{fit['Flux']}</td>
                    <td>{fit['Note']}</td>
                </tr>
                """
            html += "</tbody></table>"

        # Plots
        if has_tables:
            html += "<h3>Calibration Plots</h3>"
            for tbl in sorted(rl['tables']):
                images = table_images.get(tbl, [])
                if not images:
                    html += f"<p class='no-data'>Table {tbl} found, but no plots generated (or file missing).</p>"
                    continue
                
                html += f"<h4>Table: {tbl}</h4>"
                html += "<div class='plot-group'>"
                
                # Sort images to keep antennas in order
                images.sort()
                
                for img_path in images:
                    # Make path relative to html file for viewing
                    rel_path = os.path.relpath(img_path, os.path.dirname(os.path.abspath(html_path)))
                    fname = os.path.basename(img_path)
                    html += f"""
                    <div class="plot-card">
                        <img src="{rel_path}" alt="{fname}">
                        <div class="plot-caption">{fname}</div>
                    </div>
                    """
                html += "</div>"

        html += "</div>"

    html += "</body></html>"

    with open(html_path, 'w') as f:
        f.write(html)
    
    return html

def convert_to_pdf(html_path, pdf_path):
    try:
        from weasyprint import HTML
        print(f"[INFO] Converting HTML to PDF: {pdf_path}...")
        HTML(html_path).write_pdf(pdf_path)
        print(f"[SUCCESS] PDF created.")
    except ImportError:
        print("[WARN] 'weasyprint' not installed. Skipping PDF generation.")
        print("       (You can install it via 'pip install weasyprint' to enable PDF output)")
    except Exception as e:
        print(f"[ERROR] PDF generation failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="Generate VLBI Diagnostic Page from Log")
    parser.add_argument("logfile", help="Path to the datacheck.log file")
    args = parser.parse_args()

    if not os.path.exists(args.logfile):
        print(f"Error: File {args.logfile} not found.")
        sys.exit(1)

    # Configuration
    output_dir = "diagnostic_output"
    plots_dir = os.path.join(output_dir, "plots")
    tables_dir = "tables" # Assumed by prompt
    
    html_filename = os.path.join(output_dir, "diagnostic.html")
    pdf_filename = os.path.join(output_dir, "diagnostic.pdf")

    ensure_dir(output_dir)

    # 1. Parse Log
    print("[INFO] Parsing log file...")
    runlevels = parse_log_file(args.logfile)
    
    # 2. Generate Plots using Library
    print("[INFO] Generating plots from tables...")
    table_images_map = process_and_plot(runlevels, tables_dir, plots_dir)
    
    # 3. Write HTML
    print("[INFO] Generating HTML report...")
    generate_html(runlevels, table_images_map, args.logfile, html_filename)
    print(f"[SUCCESS] Report saved to: {html_filename}")

    # 4. Convert to PDF
    convert_to_pdf(html_filename, pdf_filename)

if __name__ == "__main__":
    main()
