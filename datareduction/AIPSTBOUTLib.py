import numpy as np
import os
import re
import sys
import argparse
import glob
import matplotlib 
matplotlib.use('Agg', force=True) #Non interactive
import matplotlib.pyplot as plt

class AIPSTBOUTTable:
    def __init__(self):
        self.header_lines = []
        self.column_meta = []
        self.data = {}
        self.pass_header_blocks = [] 
        self.pass_col_layouts = [] 

    def _parse_fortran_number(self, s):
        """Parses Fortran-formatted numbers to float. Handles 'INDE'."""
        s = s.strip()
        if not s or s == "''":
            return None 
        
        if 'INDE' in s:
            return np.nan

        try:
            return float(s.replace('D', 'E'))
        except ValueError:
            raise ValueError(f"Unable to parse number: '{s}'")

    def _format_fortran(self, val, form_str):
        """Formats value back to Fortran string based on TFORM."""
        match = re.match(r'([A-Z])(\d+)(\.(\d+))?', form_str)
        if match:
            fmt_type, width, _, precision = match.groups()
            width = int(width)
            precision = int(precision) if precision else 0
        else:
            width = 10 
            fmt_type = 'A'

        if val is None:
            return "" 
        
        if isinstance(val, (float, np.floating)) and np.isnan(val):
            inde_str = " 'INDE'"
            return f"{inde_str:<{width}}"

        if fmt_type == 'I':
            return f"{int(val):{width}d}"
        elif fmt_type in ['E', 'D']:
            fmt = f"{{:{width}.{precision}E}}"
            out_str = fmt.format(val)
            if fmt_type == 'D':
                out_str = out_str.replace('E', 'D')
            return out_str
            
        return f"{str(val):>{width}}"

    def _write_str_at(self, buffer_list, start_idx, string):
        length = len(string)
        if start_idx + length > len(buffer_list):
            buffer_list.extend([" "] * (start_idx + length - len(buffer_list)))
        for i, char in enumerate(string):
            buffer_list[start_idx + i] = char

    def read(self, filepath):
        print(f"[INFO] Reading file: {filepath}")
        with open(filepath, 'r') as f:
            lines = f.readlines()

        header_end_idx = 0
        header_dict = {}
        for i, line in enumerate(lines):
            self.header_lines.append(line)
            if line.strip() == 'END':
                header_end_idx = i
                break
            if '=' in line:
                key, rest = line.split('=', 1)
                val_part = rest.split('/', 1)[0]
                header_dict[key.strip()] = val_part.strip().strip("'")

        n_fields = int(header_dict.get('TFIELDS', 0))
        n_rows = int(header_dict.get('NAXIS2', 0))
        
        for k in range(1, n_fields + 1):
            col_info = {
                'index': k,
                'name': header_dict.get(f'TTYPE{k}', '').strip(),
                'form': header_dict.get(f'TFORM{k}', '').strip(),
                'dim': int(header_dict.get(f'TFDIM{k}', 1)),
                'raw_tbcol': int(header_dict.get(f'TBCOL{k}', 0)),
            }
            match = re.match(r'[A-Z](\d+)', col_info['form'])
            col_info['width'] = int(match.group(1)) if match else 10
            self.column_meta.append(col_info)

        for col in self.column_meta:
            if not col['name']:
                col['name'] = f"FIELD_{col['index']}"
            shape = (n_rows, col['dim']) if col['dim'] > 1 else (n_rows,)
            self.data[col['name']] = np.zeros(shape, dtype=np.float64)

        current_idx = header_end_idx + 1
        
        while current_idx < len(lines):
            line = lines[current_idx]
            if line.strip().startswith("COL. NO."):
                parts = line.replace("COL. NO.", "").split()
                col_indices = tuple(int(x) for x in parts)
                self.pass_col_layouts.append(col_indices)
                
                block_lines = []
                while current_idx < len(lines):
                    curr_line = lines[current_idx]
                    if curr_line.strip().startswith("***BEGIN*PASS***"):
                        break
                    block_lines.append(curr_line)
                    current_idx += 1
                self.pass_header_blocks.append(block_lines)
                current_idx += 1 
                
                active_cols = [self.column_meta[i-1] for i in col_indices]
                if active_cols:
                    min_tbcol = min(c['raw_tbcol'] for c in active_cols)
                    pass_offset = (min_tbcol // 1000) * 1000
                else:
                    pass_offset = 0

                pass_dim = max(c['dim'] for c in active_cols) if active_cols else 1
                row_counter = 0
                
                while current_idx < len(lines):
                    if lines[current_idx].strip().startswith("***END*PASS***"):
                        current_idx += 1
                        break
                    if row_counter >= n_rows: break

                    for sub_idx in range(pass_dim):
                        if current_idx >= len(lines): break
                        line_txt = lines[current_idx]
                        if line_txt.strip().startswith("***END*PASS***"): break 

                        for col in active_cols:
                            phy_tbcol = col['raw_tbcol'] - pass_offset
                            start = phy_tbcol - 1
                            end = start + col['width']
                            
                            if start < len(line_txt):
                                val_str = line_txt[start:end]
                                parsed_val = self._parse_fortran_number(val_str)
                                
                                if col['dim'] == 1:
                                    if sub_idx == 0 and parsed_val is not None:
                                        self.data[col['name']][row_counter] = parsed_val
                                else:
                                    if parsed_val is not None:
                                        if sub_idx < col['dim']:
                                            self.data[col['name']][row_counter, sub_idx] = parsed_val
                        current_idx += 1
                    if current_idx < len(lines) and lines[current_idx].strip().startswith("***END*PASS***"):
                        current_idx += 1
                        break
                    row_counter += 1
            else:
                current_idx += 1

    def write(self, filename):
        print(f"[INFO] Writing to {filename}")
        with open(filename, 'w') as f:
            for line in self.header_lines:
                f.write(line)
            n_rows = list(self.data.values())[0].shape[0]
            for pass_idx, col_indices in enumerate(self.pass_col_layouts):
                for h_line in self.pass_header_blocks[pass_idx]:
                    f.write(h_line)
                f.write("***BEGIN*PASS***\n")
                cols = [self.column_meta[i-1] for i in col_indices]
                if cols:
                    min_tbcol = min(c['raw_tbcol'] for c in cols)
                    pass_offset = (min_tbcol // 1000) * 1000
                else:
                    pass_offset = 0
                pass_dim = max(c['dim'] for c in cols) if cols else 1

                for r in range(n_rows):
                    row_num = r + 1
                    for sub_idx in range(pass_dim):
                        line_chars = [" "] * 132
                        self._write_str_at(line_chars, 0, f"{row_num:8d}")
                        for col in cols:
                            phy_tbcol = col['raw_tbcol'] - pass_offset
                            start = phy_tbcol - 1
                            width = col['width']
                            has_value = False
                            if col['dim'] == 1:
                                if sub_idx == 0: has_value = True
                            else:
                                if sub_idx < col['dim']: has_value = True
                            if has_value:
                                if col['dim'] == 1:
                                    val = self.data[col['name']][r]
                                else:
                                    val = self.data[col['name']][r, sub_idx]
                                val_str = self._format_fortran(val, col['form'])
                            else:
                                val_str = "''".rjust(width)
                            self._write_str_at(line_chars, start, val_str)
                        f.write("".join(line_chars).rstrip() + "\n")
                f.write("***END*PASS***\n")

    def _get_column_fuzzy(self, candidates):
        for c in candidates:
            if c in self.data: return self.data[c]
        for key in self.data.keys():
            for c in candidates:
                if c in key: return self.data[key]
        return None

    def plot(self, plot_type, xaxis='time', antennas='all', freqs='all', pols='all', out_prefix='plot', plotSeparateIFs=True):
        valid_types = ['amplitude', 'phase', 'delay', 'weight', 'bandpass']
        if plot_type not in valid_types:
            print(f"[ERROR] Invalid plot_type '{plot_type}'.")
            return

        time_data = self._get_column_fuzzy(['TIME', 'TIME MJD'])
        if time_data is None: 
            print("[ERROR] 'TIME' column missing from table.")
            return
        
        # --- Time Handling ---
        # No subtraction of t0, time is referenced to midnight of the current day in AIPS conventions
        time_hrs = time_data * 24.0

        ant_col = self._get_column_fuzzy(['ANTENNA', 'ANTENNA NO', 'ANTENNA ID', 'STATION'])
        if ant_col is None:
            print(f"[ERROR] 'ANTENNA' column missing. Keys: {list(self.data.keys())}")
            return

        freq_col = self._get_column_fuzzy(['FREQ ID', 'FQ ID'])
        if freq_col is None: 
            freq_col = np.ones(len(time_data), dtype=np.float64)

        ant_col_safe = ant_col.astype(float)
        freq_col_safe = freq_col.astype(float)

        target_ants = np.unique(ant_col_safe[~np.isnan(ant_col_safe)]).astype(int) if antennas == 'all' else np.array(antennas)
        target_freqs = np.unique(freq_col_safe[~np.isnan(freq_col_safe)]).astype(int) if freqs == 'all' else np.array(freqs)
        target_pols = [1, 2] if pols == 'all' else pols

        print(f"[INFO] Generating {plot_type} vs {xaxis} plots for {len(target_ants)} antennas...")

        # --- Harmonised Font Sizes ---
        plt.rcParams.update({
            'font.size': 12, 
            'axes.labelsize': 14, 
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12
        })
        
        # --- Polarization Color Map ---
        pol_colors = {1: 'tab:blue', 2: 'tab:orange'}
        
        # --- Pre-calculation for Layout ---
        relevant_cols_map = {
            'amplitude': ['REAL', 'CORR'], 'phase': ['REAL', 'CORR'], 'bandpass': ['REAL', 'CORR'],
            'delay': ['DELAY'], 'weight': ['WEIGHT']
        }
        search_keys = relevant_cols_map.get(plot_type, [])
        max_ifs = 1
        for meta in self.column_meta:
            for key in search_keys:
                if key in meta['name']: 
                    max_ifs = max(max_ifs, meta['dim'])
        
        do_separate_ifs = (plotSeparateIFs and max_ifs > 1 and xaxis == 'time' and plot_type != 'bandpass')

        for ant in target_ants:
            
            # --- Figure Setup ---
            if do_separate_ifs:
                fig, axes = plt.subplots(max_ifs, 1, sharex=True, figsize=(10, 2.5 * max_ifs + 1))
                if max_ifs == 1: axes = [axes] 
            elif plot_type == 'bandpass':
                fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 10))
                ax_main = ax1 
            else:
                fig, ax_main = plt.subplots(figsize=(10, 6))

            plotted_something = False

            for freq in target_freqs:
                mask = (ant_col_safe == ant) & (freq_col_safe == freq)
                if not np.any(mask): continue

                for pol in target_pols:
                    curr_color = pol_colors.get(pol, 'gray')

                    def get_data_col(base, p_idx, try_alternates=None):
                        k1 = f"{base} {p_idx}"
                        if k1 in self.data: return self.data[k1]
                        k2 = f"{base}{p_idx}"
                        if k2 in self.data: return self.data[k2]
                        if p_idx == 1 and base in self.data: return self.data[base]
                        if try_alternates:
                            for alt in try_alternates:
                                k_alt1 = f"{alt} {p_idx}"
                                if k_alt1 in self.data: return self.data[k_alt1]
                                k_alt2 = f"{alt}{p_idx}"
                                if k_alt2 in self.data: return self.data[k_alt2]
                                if p_idx == 1 and alt in self.data: return self.data[alt]
                        return None

                    x_vals = None
                    y_vals_1 = None 
                    y_vals_2 = None 

                    if plot_type == 'delay':
                        y_raw = get_data_col("DELAY", pol)
                        if y_raw is not None: 
                            y_vals_1 = y_raw[mask] * 1e9 # s -> ns

                    elif plot_type == 'weight':
                        y_raw = get_data_col("WEIGHT", pol)
                        if y_raw is not None: y_vals_1 = y_raw[mask]

                    elif plot_type in ['amplitude', 'phase', 'bandpass']:
                        r_d = get_data_col("REAL", pol, try_alternates=["CORR"])
                        i_d = get_data_col("IMAG", pol)
                        
                        if r_d is not None:
                            r_d = r_d[mask]
                            if i_d is None: i_d = np.zeros_like(r_d)
                            else: i_d = i_d[mask]

                            amp = np.sqrt(r_d**2 + i_d**2)
                            phs = np.degrees(np.arctan2(i_d, r_d))

                            if plot_type == 'amplitude': y_vals_1 = amp
                            elif plot_type == 'phase': y_vals_1 = phs
                            elif plot_type == 'bandpass':
                                y_vals_1 = amp
                                y_vals_2 = phs

                    if y_vals_1 is None: continue

                    # --- Plotting Logic ---
                    if xaxis == 'time':
                        x_vals = time_hrs[mask]
                        
                        if y_vals_1.ndim > 1:
                            # Multi-IF Data
                            n_ifs_actual = y_vals_1.shape[1]
                            for if_idx in range(n_ifs_actual):
                                if if_idx >= max_ifs and do_separate_ifs: break 
                                y_sub = y_vals_1[:, if_idx]
                                
                                if do_separate_ifs:
                                    lbl = f"Freq {freq} Pol {pol}"
                                else:
                                    lbl = f"Freq {freq} Pol {pol} IF {if_idx+1}"

                                target_ax = axes[if_idx] if do_separate_ifs else ax_main

                                if len(x_vals) == len(y_sub):
                                    if plot_type == 'bandpass':
                                        y_sub2 = y_vals_2[:, if_idx]
                                        ax1.scatter(x_vals, y_sub, label=lbl, s=35, alpha=0.7, color=curr_color)
                                        ax2.scatter(x_vals, y_sub2, label=lbl, s=35, alpha=0.7, color=curr_color)
                                    else:
                                        target_ax.scatter(x_vals, y_sub, label=lbl, s=35, alpha=0.7, color=curr_color)
                                    
                                    plotted_something = True

                        else:
                            # Scalar Data
                            lbl = f"Freq {freq} Pol {pol}"
                            target_ax = axes[0] if do_separate_ifs else ax_main
                            
                            if len(x_vals) == len(y_vals_1):
                                if plot_type == 'bandpass':
                                    ax1.scatter(x_vals, y_vals_1, label=lbl, s=35, alpha=0.7, color=curr_color)
                                    ax2.scatter(x_vals, y_vals_2, label=lbl, s=35, alpha=0.7, color=curr_color)
                                else:
                                    target_ax.scatter(x_vals, y_vals_1, label=lbl, s=35, alpha=0.7, color=curr_color)
                                plotted_something = True

                    elif xaxis == 'frequency':
                        if y_vals_1.ndim > 1:
                            n_rows = y_vals_1.shape[0]
                            channels = np.arange(y_vals_1.shape[1]) + 1
                            lbl_base = f"Freq {freq} Pol {pol}"
                            for r in range(n_rows):
                                lbl = lbl_base if r==0 else "" 
                                if plot_type == 'bandpass':
                                    ax_main.plot(channels, y_vals_1[r, :], label=lbl, alpha=0.5, linewidth=0.8, color=curr_color)
                                    ax2.plot(channels, y_vals_2[r, :], label=lbl, alpha=0.5, linewidth=0.8, color=curr_color)
                                else:
                                    ax_main.plot(channels, y_vals_1[r, :], label=lbl, alpha=0.5, linewidth=0.8, color=curr_color)
                            plotted_something = True

            if plotted_something:
                # --- Post-Plot Formatting ---
                ylabel_map = {
                    'amplitude': "Amp", 
                    'phase': "Phase (deg)", 
                    'delay': "Delay (ns)", 
                    'weight': "Weight"
                }
                base_ylab = ylabel_map.get(plot_type, "Val")
                
                # --- Helper: Colored Text Annotation ---
                def add_colored_pol_labels(target_ax, target_fig=None, is_multi=False):
                    # Add "Pol 1" in Blue and "Pol 2" in Orange as text objects
                    # Positioned top-left with opaque background to prevent border bleed-through
                    
                    # Common bbox properties for opacity
                    box_props = dict(facecolor='white', alpha=1.0, edgecolor='none', pad=2.0)
                    
                    if is_multi and target_fig:
                        # For multi-panel
                        target_ax.text(0.02, 0.96, "Pol 1", color='tab:blue', transform=target_ax.transAxes, 
                                     fontweight='bold', fontsize=15, bbox=box_props)
                        target_ax.text(0.02, 0.86, "Pol 2", color='tab:orange', transform=target_ax.transAxes, 
                                     fontweight='bold', fontsize=15, bbox=box_props)
                    else:
                        target_ax.text(0.02, 0.96, "Pol 1", color='tab:blue', transform=target_ax.transAxes, 
                                     fontweight='bold', fontsize=15, bbox=box_props)
                        target_ax.text(0.02, 0.88, "Pol 2", color='tab:orange', transform=target_ax.transAxes, 
                                     fontweight='bold', fontsize=15, bbox=box_props)

                if do_separate_ifs:
                    axes[0].set_title(f"Antenna {ant}: {plot_type.capitalize()} vs Time")
                    
                    # Add Colored Annotation
                    add_colored_pol_labels(axes[0], fig, is_multi=True)

                    for i, ax in enumerate(axes):
                        ax.grid(True, alpha=0.3)
                        ax.set_ylabel(f"{base_ylab}")
                        # IF Label
                        ax.text(0.99, 0.95, f"IF {i+1}", transform=ax.transAxes, 
                                ha='right', va='top', fontweight='bold', 
                                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
                        
                        if i < len(axes) - 1:
                            plt.setp(ax.get_xticklabels(), visible=False)
                    
                    axes[-1].set_xlabel("Time (Hours)")
                    
                    # Common Legend (Deduplicated)
                    lines_labels = [ax.get_legend_handles_labels() for ax in axes]
                    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
                    unique_dict = dict(zip(labels, lines))
                    fig.legend(unique_dict.values(), unique_dict.keys(), 
                               bbox_to_anchor=(1.0, 1), loc='upper left', markerscale=1.5)

                else:
                    ax_main.set_title(f"Antenna {ant}")
                    
                    # Add Colored Annotation
                    add_colored_pol_labels(ax_main, is_multi=False)

                    if xaxis == 'time':
                        xlabel_text = "Time (Hours)"
                    else:
                        xlabel_text = "Channel Index"

                    if plot_type == 'bandpass':
                        ax_main.set_ylabel("Amplitude")
                        ax2.set_ylabel("Phase (deg)")
                        ax2.set_xlabel(xlabel_text)
                        ax_main.grid(True, alpha=0.3)
                        ax2.grid(True, alpha=0.3)
                        # LEGEND REMOVED
                    else:
                        ax_main.set_xlabel(xlabel_text)
                        ax_main.set_ylabel(base_ylab)
                        ax_main.grid(True, alpha=0.3)
                        h, l = ax_main.get_legend_handles_labels()
                        by_label = dict(zip(l, h))
                        ax_main.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=1.5)

                plt.tight_layout()
                out_name = f"{out_prefix}_{plot_type}_vs_{xaxis}_ant{ant}.png"
                plt.savefig(out_name, dpi=150)
                print(f"  Saved {out_name}")
                plt.close(fig)
            else:
                plt.close(fig)

def test_reproduction(input_file):
    output_file = "test.tbout.txt"
    if not os.path.exists(input_file):
        print(f"[ERROR] Input file {input_file} not found.")
        sys.exit(1)
    sn = AIPSTBOUTTable()
    try:
        sn.read(input_file)
        sn.write(output_file)
    except Exception as e:
        print(f"[FATAL] {e}")
        return
    with open(input_file, 'r') as f1, open(output_file, 'r') as f2:
        c1 = f1.readlines()
        c2 = f2.readlines()
    if len(c1) != len(c2):
        print(f"[FAIL] Line count: {len(c1)} vs {len(c2)}")
        return
    mismatches = 0
    for i, (l1, l2) in enumerate(zip(c1, c2)):
        if l1.rstrip() != l2.rstrip():
            mismatches += 1
    if mismatches == 0:
        print("SUCCESS: Files are character identical.")
    else:
        print(f"FAILURE: {mismatches} mismatches.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read, Verify and Plot AIPS TBOUT.")
    parser.add_argument('filename', nargs='?', help='Path to file')
    args = parser.parse_args()
    target_file = args.filename
    if not target_file:
        found = glob.glob("*.sn") + glob.glob("*.bp") + glob.glob("*.tbout")
        if not found: sys.exit(1)
        target_file = found[0]
        print(f"[INFO] Auto-detected: {target_file}")
    else:
        print(f"[INFO] Using: {target_file}")
    test_reproduction(target_file)
