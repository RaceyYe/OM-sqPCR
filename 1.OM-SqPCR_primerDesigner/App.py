import tkinter as tk
from tkinter import ttk, filedialog, messagebox, Toplevel
from integrated_execution import execution
from modules.pairing import generate_combinations
import argparse
import re
import os
import pandas as pd
import csv
import traceback
import threading
import signal
import multiprocessing as mp
from contextlib import contextmanager
from outputs_integration import integrate_outputs

# Ensure working directory is the script's directory
type_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(type_dir)

# Handle broken pipe errors
def signal_handler(signal, frame):
    pass

signal.signal(signal.SIGPIPE, signal_handler)

@contextmanager
def process_pool(max_workers=None):
    pool = mp.Pool(processes=max_workers)
    try:
        yield pool
    finally:
        pool.close()
        pool.join()


def parse_fasta(text):
    entries = []
    header = None
    seq_lines = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header and seq_lines:
                entries.append((header, ''.join(seq_lines)))
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line)
    if header and seq_lines:
        entries.append((header, ''.join(seq_lines)))
    return entries


class PrimerAnalysisBatchApp:
    def __init__(self, root, parameter_option=False):
        self.root = root
        self.parameter_option = parameter_option  # control visibility of parameter button/icon
        self.root.title("Batch Primer Designer")
        self.root.geometry("820x600")

        # Styling
        self.style = ttk.Style()
        self.style.configure("TFrame", background="#f9f9f9")
        self.style.configure("TLabel", background="#f9f9f9", font=("Arial", 10))
        self.style.configure("Header.TLabel", font=("Arial", 12, "bold"))
        self.style.configure("TButton", font=("Arial", 10), padding=5)

        # Data holders
        self.srna_data = []
        self.stemloop_data = []
        self.output_dir = None
        self.cancel_processing = False
        self.error_log = []

        # Default parameters
        self.parameters = {
            "tm_base": 59,
            "tm_var": 2,
            "tm_base_probe_1": 59,
            "tm_var_probe_1": 2,
            "tm_base_probe_2": 71,
            "tm_var_probe_2": 2,
            "dg_min": -4.5, # (gs), # -5 #-5, #
            "dg_max": 10,
            "dg_min_FR": -3.9,
            "dg_max_FR": 10,
            "dg_min_PP":  -3.9, # (gs), # -4, 3.5X #-4.5, #
            "dg_max_PP": 10,
            "temp": 37,
            "min_5ext": 6,
            "max_5ext": 6,
            "GAP": 30,
            "PR_gap_max": 15,
            "min_RP_len": 18,
            "max_RP_len": 30,
            "min_probe_len": 18,
            "max_probe_len": 30
        }

        self.create_widgets()

    def create_widgets(self):
        main_frame = ttk.Frame(self.root, style="TFrame", padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Top controls
        top_frame = ttk.Frame(main_frame, style="TFrame")
        top_frame.pack(fill=tk.X, pady=(0, 10))
        # Only show Parameters button/icon if parameter_option is True
        if self.parameter_option:
            param_btn = ttk.Button(top_frame, text="Parameters", command=self.open_parameter_window)
            param_btn.pack(side=tk.LEFT, padx=5)
        # always show help
        help_btn = ttk.Button(top_frame, text="Help", command=self.show_help)
        help_btn.pack(side=tk.LEFT, padx=5)
        out_btn = ttk.Button(top_frame, text="Select Output Dir", command=self.select_output_directory)
        out_btn.pack(side=tk.LEFT, padx=5)
        self.out_label = ttk.Label(top_frame, text="No output directory selected", style="TLabel")
        self.out_label.pack(side=tk.LEFT, padx=10)

        # Input sections
        ttk.Label(main_frame, text="sRNA FASTA:", style="Header.TLabel").pack(anchor=tk.W, pady=(10, 0))
        btn_srna_file = ttk.Button(main_frame, text="Load sRNA File", command=self.load_srna_file)
        btn_srna_file.pack(anchor=tk.W, pady=5)
        self.srna_text = tk.Text(main_frame, height=6)
        self.srna_text.pack(fill=tk.X)

        ttk.Label(main_frame, text="Stemloop FASTA:", style="Header.TLabel").pack(anchor=tk.W, pady=(10, 0))
        btn_stem_file = ttk.Button(main_frame, text="Load Stemloop File", command=self.load_stemloop_file)
        btn_stem_file.pack(anchor=tk.W, pady=5)
        self.stemloop_text = tk.Text(main_frame, height=6)
        # insert default sequences
        # default = ">stemloop_1\nGTTGGCTCTGGTGCAGGGTCCGAGGTATTCGCACCAGAGCCAAC\n>stemloop_2\nGTCGTATCCAGTGCGTGTCGTGGAGTCGGCAATTGCACTGGATACGAC\n>stemloop_3\nCTCAACTGGTGTCGTGGAGTCAGTGCAATTCCAGTTGAG\n>stemloop_4\nGTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC"
        default = ">stemloop_1\nGTTGGCTCTGGTGCAGGGTCCGAGGTATTCGCACCAGAGCCAAC\n>stemloop_2\nGAAAGAAGGCGAGGAGCAGATCGAGGAAGAAGACGGAAGAATGTGCGTCTCGCCTTCTTTC\n>stemloop_3\nCTCAACTGGTGTCGTGGAGTCAGTGCAATTCCAGTTGAG\n>stemloop_4\nGTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC"
        self.stemloop_text.insert('1.0', default)
        self.stemloop_text.pack(fill=tk.X)

        # Criteria selection frame
        ttk.Label(main_frame, text="Selection Criteria:", style="Header.TLabel").pack(anchor=tk.W, pady=(10, 0))

        criteria_frame = ttk.Frame(main_frame, padding=10, style="TFrame")
        criteria_frame.pack(fill=tk.X)

        # Create Boolean variables for each checkbox
        self.criteria_var_tm = tk.BooleanVar(value=True)
        self.criteria_var_lp = tk.BooleanVar(value=False)
        self.criteria_var_dg = tk.BooleanVar(value=True)

        # Create checkboxes without mutual exclusivity
        ttk.Checkbutton(criteria_frame, text="TM", variable=self.criteria_var_tm).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(criteria_frame, text="LP", variable=self.criteria_var_lp).pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(criteria_frame, text="DG", variable=self.criteria_var_dg).pack(side=tk.LEFT, padx=5)

        # Forbidden criteria frame
        forbidden_frame = ttk.Frame(main_frame, padding=10, style="TFrame")
        forbidden_frame.pack(fill=tk.X, pady=10)

        self.forbidden_var = tk.BooleanVar(value=True)
        forbidden_check = ttk.Checkbutton(forbidden_frame, text="Apply forbidden criteria", variable=self.forbidden_var)
        forbidden_check.pack(side=tk.LEFT, anchor=tk.W)

        self.forbidden_entry = ttk.Entry(forbidden_frame, width=75)
        self.forbidden_entry.pack(side=tk.LEFT, padx=10)
        self.forbidden_entry.insert(0, "GTATGTGGGAGCCCACACTCTGGTCGACAGATACGAATATCTGGACCCGACCGCTCCCACATAC")

        def toggle_forbidden_entry(*args):
            if self.forbidden_var.get():
                self.forbidden_entry.configure(state='normal')
            else:
                self.forbidden_entry.configure(state='disabled')

        self.forbidden_var.trace_add('write', toggle_forbidden_entry)
        toggle_forbidden_entry()

        # Run button
        run_btn = ttk.Button(main_frame, text="Run Batch Analysis", command=self.run_batch)
        run_btn.pack(pady=20)

    def open_parameter_window(self):
        win = Toplevel(self.root)
        win.title("Parameter Settings")
        win.geometry("450x800")
        win.grab_set()
        instr = ttk.Label(win, text="Adjust numerical parameters:", style="Header.TLabel")
        instr.pack(pady=10)
        entries = {}
        for key, val in self.parameters.items():
            frame = ttk.Frame(win, padding=5)
            frame.pack(fill=tk.X)
            lbl = ttk.Label(frame, text=key.replace('_', ' ').capitalize() + ":", width=20)
            lbl.pack(side=tk.LEFT)
            var = tk.StringVar(value=str(val))
            ent = ttk.Entry(frame, textvariable=var, width=10)
            ent.pack(side=tk.LEFT)
            entries[key] = var
        btn_frame = ttk.Frame(win)
        btn_frame.pack(pady=15)


        def save_params():
            try:
                for k, v in entries.items():
                    self.parameters[k] = float(v.get())
                win.destroy()
            except ValueError:
                messagebox.showerror("Error", "Invalid parameter value.")

        ttk.Button(btn_frame, text="Save", command=save_params).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=win.destroy).pack(side=tk.LEFT, padx=5)

    def get_selected_criteria(self):
        """Collects all selected criteria checkboxes."""
        criteria = []

        # Add each selected criterion to the list
        if self.criteria_var_tm.get():
            criteria.append("TM")
        if self.criteria_var_lp.get():
            criteria.append("LP")
        if self.criteria_var_dg.get():
            criteria.append("DG")

        # If no criteria selected, default to TM
        if not criteria:
            criteria = ["TM"]

        # Get forbidden sequence if enabled
        gf_seq = None
        if self.forbidden_var.get():
            gf_seq = self.forbidden_entry.get().strip()

        return criteria, gf_seq

    def show_help(self):
        win = Toplevel(self.root)
        win.title("Help")
        win.geometry("600x400")
        text = (
            "# Batch Primer Designer Help\n"
            "\n"
            "Load your sRNA and stemloop sequences in FASTA format (either paste or load files),\n"
            "then click **Run Batch Analysis**.\n"
            "\n"
            "**Parameters** allows you to tweak TM, ΔG, lengths, and temperature ranges.\n"
            "**Output Dir** sets where each Task_* subfolder will be created.\n"
            "The program will pair the shortest sRNA with the longest stemloop sequences.\n"
            "After completion, an Excel file with integrated results will be generated.\n"
        )
        txt = tk.Text(win, wrap=tk.WORD)
        txt.insert('1.0', text)
        txt.config(state=tk.DISABLED)
        txt.pack(fill=tk.BOTH, expand=True)

    def load_srna_file(self):
        path = filedialog.askopenfilename(filetypes=[('FASTA', '*.fa *.fasta'), ('All', '*.*')])
        if path:
            with open(path) as f:
                self.srna_text.delete('1.0', tk.END)
                self.srna_text.insert('1.0', f.read())

    def load_stemloop_file(self):
        path = filedialog.askopenfilename(filetypes=[('FASTA', '*.fa *.fasta'), ('All', '*.*')])
        if path:
            with open(path) as f:
                self.stemloop_text.delete('1.0', tk.END)
                self.stemloop_text.insert('1.0', f.read())

    def select_output_directory(self):
        d = filedialog.askdirectory()
        if d:
            self.output_dir = d
            self.out_label.config(text=f"Output: {os.path.basename(d)}")

    def prioritize_sequences(self, srna_entries, stem_entries):
        """
        Sort sRNAs by length (ascending) and stemloops by length (descending).
        Return both the prioritized sequence lists.
        """
        # Create tuples of (header, sequence, length) for sorting
        srna_with_len = [(header, seq, len(seq)) for header, seq in srna_entries]
        stem_with_len = [(header, seq, len(seq)) for header, seq in stem_entries]

        # Sort sRNAs by length (ascending) and stemloops by length (descending)
        sorted_srna = sorted(srna_with_len, key=lambda x: x[2])
        sorted_stem = sorted(stem_with_len, key=lambda x: x[2], reverse=True)

        return sorted_srna, sorted_stem

    def generate_summary_excel(self, tasks_run):
        """
        Generate a summary Excel file from all the Integrated_df.csv files.
        Memory-efficient implementation with proper RT_primer column handling.
        """
        try:
            # Create Excel writer with xlsxwriter engine
            excel_path = os.path.join(self.output_dir, "Summary_Results.xlsx")
            with pd.ExcelWriter(excel_path, engine='xlsxwriter') as writer:
                all_top5_data = []

                # Parse FASTA input at the beginning - we'll use this to build RT_primers
                srna_entries = parse_fasta(self.srna_text.get('1.0', 'end'))
                stem_entries = parse_fasta(self.stemloop_text.get('1.0', 'end'))

                # Create dictionaries for quick lookup
                srna_dict = {header: seq for header, seq in srna_entries}
                stem_dict = {header: seq for header, seq in stem_entries}

                # Process each task's data
                for task, path in tasks_run:
                    csv_path = os.path.join(path, "Integrated_df.csv")
                    if not os.path.exists(csv_path):
                        continue

                    try:
                        # Extract sRNA and stemloop names from task
                        # Task format: "Task_*_sRNA_header_x_stemloop_header"
                        rt_primer = None
                        srna_seq = None
                        stem_seq = None

                        # First find which sRNA and stemloop match this task
                        for srna_header, seq in srna_dict.items():
                            if srna_header in task:
                                srna_seq = seq
                                break

                        for stem_header, seq in stem_dict.items():
                            if stem_header in task:
                                stem_seq = seq
                                break

                        # Create RT_primer if both sequences were found
                        if srna_seq and stem_seq:
                            reversed_srna_tail = srna_seq[-6:][::-1]  # Get last 6 nucleotides and reverse
                            rt_primer = stem_seq + reversed_srna_tail

                        # Read column headers first to identify sort column
                        headers = None
                        sort_col = None
                        with open(csv_path, 'r') as f:
                            reader = csv.reader(f)
                            headers = next(reader, None)
                            if headers:
                                # Determine sorting column
                                for potential_col in ['TM_PP', 'overall_assess', 'OverallScore', 'Score']:
                                    if potential_col in headers:
                                        sort_col = potential_col
                                        break

                        if headers and sort_col:
                            # Read data in chunks for memory efficiency
                            chunk_size = 500
                            top5_rows = []

                            for chunk in pd.read_csv(csv_path, chunksize=chunk_size):
                                if sort_col in chunk.columns:
                                    # Sort chunk and get top 5
                                    sorted_chunk = chunk.sort_values(by=sort_col, ascending=False)
                                    chunk_top5 = sorted_chunk.head(5)

                                    # Update overall top 5
                                    if len(top5_rows) < 5:
                                        top5_rows.extend(chunk_top5.to_dict('records'))
                                        top5_rows = sorted(top5_rows, key=lambda x: x.get(sort_col, 0), reverse=True)[
                                                    :5]
                                    else:
                                        min_value = min(top5_rows, key=lambda x: x.get(sort_col, 0)).get(sort_col, 0)
                                        for _, row in chunk_top5.iterrows():
                                            if row[sort_col] > min_value:
                                                top5_rows.append(row.to_dict())
                                                top5_rows = sorted(top5_rows, key=lambda x: x.get(sort_col, 0),
                                                                   reverse=True)[:5]
                                                min_value = min(top5_rows, key=lambda x: x.get(sort_col, 0)).get(
                                                    sort_col, 0)

                            if top5_rows:
                                # Convert to DataFrame and add task info
                                top5_df = pd.DataFrame(top5_rows)
                                top5_df['Task'] = task

                                # Add RT_primer column with the calculated value
                                top5_df['RT_primer'] = rt_primer if rt_primer else "Not found"

                                # Add to overall results
                                all_top5_data.append(top5_df)

                                # Create valid sheet name
                                sheet_name = f"Top5_{task[-15:]}" if len(task) > 15 else f"Top5_{task}"
                                sheet_name = re.sub(r'[\\/*?:\[\]]', '_', sheet_name)[:31]

                                # Write to Excel
                                top5_df.to_excel(writer, sheet_name=sheet_name, index=False)

                    except Exception as e:
                        self.error_log.append(f"Error processing {csv_path}: {str(e)}")
                        traceback.print_exc()

                # Create summary sheet with all top 5 results
                if all_top5_data:
                    try:
                        # Combine all data frames
                        combined_df = pd.concat(all_top5_data, ignore_index=True)

                        # Make sure the 'RT_primer' column exists
                        if 'RT_primer' not in combined_df.columns:
                            combined_df['RT_primer'] = "Not available"

                        # Ensure RT_primer is the last column for better visibility
                        cols = [col for col in combined_df.columns if col != 'RT_primer'] + ['RT_primer']
                        combined_df = combined_df[cols]

                        # Write the All_Top5 sheet
                        combined_df.to_excel(writer, sheet_name='All_Top5', index=False)

                    except Exception as e:
                        self.error_log.append(f"Error creating summary sheet: {str(e)}")
                        traceback.print_exc()

            return excel_path
        except Exception as e:
            self.error_log.append(f"Overall summary generation error: {str(e)}")
            traceback.print_exc()
            return None

    def execute_task(self, task_info, task_path, task_name, args, update_callback=None):
        """Execute a single task with proper error handling."""
        try:
            # Create task directory
            os.makedirs(task_path, exist_ok=True)

            # Execute task with timeout control
            execution_complete = threading.Event()
            execution_result = {"success": False, "error": None}

            def run_execution():
                try:
                    executor = execution(args, task_name,args.criteria, args.GF_Seq)
                    executor.run_design()
                    execution_result["success"] = True
                except Exception as e:
                    execution_result["error"] = str(e)
                    traceback.print_exc()
                finally:
                    execution_complete.set()

            # Run in separate thread
            thread = threading.Thread(target=run_execution)
            thread.daemon = True
            thread.start()

            # Wait for completion with timeout
            timeout_seconds = 30000  # 5 minutes
            if not execution_complete.wait(timeout_seconds):
                return False, f"Task {task_name} timed out after {timeout_seconds} seconds"

            if not execution_result["success"]:
                return False, f"Task {task_name} failed: {execution_result['error']}"

            return True, task_path

        except Exception as e:
            return False, f"Task {task_name} exception: {str(e)}"

    def run_batch(self):
        # Reset error log
        self.error_log = []
        self.cancel_processing = False

        # Parse input sequences
        srna_entries = parse_fasta(self.srna_text.get('1.0', 'end'))
        stem_entries = parse_fasta(self.stemloop_text.get('1.0', 'end'))

        # Validate inputs
        if not srna_entries or not stem_entries:
            messagebox.showerror("Error", "Please provide both sRNA and stemloop FASTA inputs.")
            return
        if not self.output_dir:
            messagebox.showerror("Error", "Please select an output directory.")
            return

        # Prioritize sequences based on length
        sorted_srna, sorted_stem = self.prioritize_sequences(srna_entries, stem_entries)

        tasks_run = []
        pairs_to_process = []

        # First, pair the shortest sRNA with the longest stemloop
        if sorted_srna and sorted_stem:
            shortest_srna_header, shortest_srna_seq, _ = sorted_srna[0]
            longest_stem_header, longest_stem_seq, _ = sorted_stem[0]

            # Add the special pair
            pairs_to_process.append({
                'sRNA_header': shortest_srna_header,
                'sRNA_sequence': shortest_srna_seq,
                'Stemloop_header': longest_stem_header,
                'Stemloop_sequence': longest_stem_seq,
                'Task': f"Priority_Short{len(shortest_srna_seq)}_Long{len(longest_stem_seq)}"
            })

            # Remove these from further combinations
            sorted_srna = sorted_srna[1:]
            sorted_stem = sorted_stem[1:]

        # Generate remaining combinations
        if sorted_srna and sorted_stem:
            srna_seqs = {seq for _, seq, _ in sorted_srna}
            stem_seqs = {seq for _, seq, _ in sorted_stem}
            combo_df = generate_combinations(srna_seqs, stem_seqs)

            # Add the regular combinations
            for _, row in combo_df.iterrows():
                pairs_to_process.append({
                    'sRNA_sequence': row['sRNA_sequence'],
                    'Stemloop_sequence': row['Stemloop_sequence'],
                    'Task': row['Task']
                })

        # Create progress window
        progress_window = Toplevel(self.root)
        progress_window.title("Processing Batch")
        progress_window.geometry("400x150")

        # Progress indicators
        progress_label = ttk.Label(progress_window, text="Starting batch analysis...")
        progress_label.pack(pady=10)
        progress_bar = ttk.Progressbar(progress_window, length=300, mode="determinate")
        progress_bar.pack(pady=10)
        status_label = ttk.Label(progress_window, text="")
        status_label.pack(pady=5)

        # Cancel button
        cancel_btn = ttk.Button(progress_window, text="Cancel",
                                command=lambda: setattr(self, 'cancel_processing', True))
        cancel_btn.pack(pady=5)

        total_tasks = len(pairs_to_process)

        # Thread function for batch processing
        def process_batch():
            successful_tasks = []

            for i, pair_info in enumerate(pairs_to_process):
                if self.cancel_processing:
                    self.root.after(0, lambda: status_label.config(text="Processing cancelled"))
                    break

                # Update progress in the GUI thread
                self.root.after(0, lambda i=i, total=total_tasks: [
                    progress_bar.config(value=((i / total) * 100)),
                    progress_label.config(text=f"Processing task {i + 1} of {total}...")
                ])

                # Determine task name
                task_name = str(pair_info['Task'])
                if 'Priority' in task_name:
                    task_prefix = "Priority"
                else:
                    task_prefix = f"Task_{i + 1}"

                # Get header information for better naming
                srna_header = pair_info.get('sRNA_header', '')
                stemloop_header = pair_info.get('Stemloop_header', '')

                # If headers aren't available (for regular combinations), try to find them
                if not srna_header or not stemloop_header:
                    # Find headers by matching sequences
                    srna_seq = pair_info['sRNA_sequence']
                    stemloop_seq = pair_info['Stemloop_sequence']

                    # Find matching headers from original entries
                    for header, seq in srna_entries:
                        if seq == srna_seq:
                            srna_header = header
                            break

                    for header, seq in stem_entries:
                        if seq == stemloop_seq:
                            stemloop_header = header
                            break

                # Clean headers for directory naming (remove problematic characters)
                srna_header_clean = re.sub(r'[\\/*?:"<>|]', '_', srna_header)[:15] if srna_header else 'sRNA'
                stemloop_header_clean = re.sub(r'[\\/*?:"<>|]', '_', stemloop_header)[:15] if stemloop_header else 'SL'

                # Create directory name with headers
                task = f"{task_prefix}_{srna_header_clean}_x_{stemloop_header_clean}"

                out_path = os.path.join(self.output_dir, task)

                # Get selected criteria and forbidden sequence
                criteria, gf_seq = self.get_selected_criteria()

                # Prepare arguments
                args = argparse.Namespace(
                    upper_strands=pair_info['sRNA_sequence'],
                    lower_strands=pair_info['Stemloop_sequence'],
                    output=out_path,
                    criteria=criteria,  # Pass selected criteria
                    GF_Seq=gf_seq,  # Pass forbidden sequence if specified
                    **self.parameters
                )

                # Update status
                self.root.after(0, lambda task=task: status_label.config(text=f"Running: {task}"))

                # Execute task
                success, result = self.execute_task(pair_info, out_path, task, args)

                if success:
                    successful_tasks.append((task, result))
                else:
                    self.error_log.append(result)

            # Generate summary Excel if tasks were completed
            summary_path = None
            if successful_tasks:
                self.root.after(0, lambda: [
                    progress_label.config(text="Generating summary report..."),
                    status_label.config(text="Please wait...")
                ])
                summary_path = self.generate_summary_excel(successful_tasks)

            # Show final report
            self.root.after(0, lambda: [
                progress_window.destroy(),
                self.show_completion_report(successful_tasks, summary_path)
            ])

        # Start processing in a separate thread
        processing_thread = threading.Thread(target=process_batch)
        processing_thread.daemon = True
        processing_thread.start()

    def show_completion_report(self, tasks_run, excel_path=None):
        """Display a report of the completed tasks."""
        rpt = Toplevel(self.root)
        rpt.title("Batch Analysis Report")
        rpt.geometry("600x400")

        # Create scrollable text area
        txt_frame = ttk.Frame(rpt)
        txt_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        scrollbar = ttk.Scrollbar(txt_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        txt = tk.Text(txt_frame, wrap=tk.WORD, yscrollcommand=scrollbar.set)
        txt.pack(fill=tk.BOTH, expand=True)
        scrollbar.config(command=txt.yview)

        # Summary information
        txt.insert('end', f"Batch completed: {len(tasks_run)} tasks run\n")
        txt.insert('end', f"Output directory: {self.output_dir}\n\n")

        if excel_path:
            txt.insert('end', f"Summary Excel file: {os.path.basename(excel_path)}\n\n")

        # List of completed tasks
        txt.insert('end', "Completed Tasks:\n")
        for t, p in tasks_run:
            txt.insert('end', f"• {t}\n")

        # Display errors if any
        if self.error_log:
            txt.insert('end', "\n\nErrors encountered:\n")
            for error in self.error_log:
                txt.insert('end', f"• {error}\n")

        txt.config(state=tk.DISABLED)
        messagebox.showinfo("Batch Complete", f"Batch analysis finished with {len(tasks_run)} successful tasks.")

        # after you’ve collected `tasks_run` and know your output folder:
        csv_path, excel_path = integrate_outputs(self.output_dir, self.output_dir)
        messagebox.showinfo("Integration Complete",
                            f"Summary written to:\n  {csv_path}\nExcel written to:\n  {excel_path}")


if __name__ == '__main__':
    # Set up multiprocessing to work with tkinter
    mp.set_start_method('spawn', force=True)

    # Create and start application
    root = tk.Tk()
    app = PrimerAnalysisBatchApp(root)
    root.mainloop()
