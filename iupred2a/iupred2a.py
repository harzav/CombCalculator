#!/usr/bin/python3

import sys
import os
import iupred2a_lib

# Define your input sequence and mode here
input_sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQVTDSQMTKTAYIAKQRQISFVK"
prediction_mode = "glob"  # Can be 'long', 'short', or 'glob'
enable_anchor = True  # Set to True if you want to enable ANCHOR2 prediction

def run_iupred(sequence, mode="long", anchor=False):
    PATH = os.path.dirname(os.path.realpath(__file__))

    if not os.path.isdir(PATH):
        sys.exit(f'Data directory not found at {PATH}!')

    if mode not in ['short', 'long', 'glob']:
        sys.exit(f'Invalid IUPred2 mode: {mode}. Must be "short", "long", or "glob".')

    # Perform IUPred2 analysis
    iupred2_result = iupred2a_lib.iupred(sequence, mode)
    anchor2_res = None

    if anchor and mode == "long":
        anchor2_res = iupred2a_lib.anchor2(sequence)

    # Output results
    print("""# IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding
# Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi
# Nucleic Acids Research 2018;46(W1):W329-W337.
#
# Prediction type: {}
# Prediction output""".format(mode))

    if mode == 'glob':
        print(iupred2_result[1])

    if anchor:
        print("# POS\tRES\tIUPRED2\tANCHOR2")
    else:
        print("# POS\tRES\tIUPRED2")

    for pos, residue in enumerate(sequence):
        print(f"{pos + 1}\t{residue}\t{iupred2_result[0][pos]:.4f}", end="")
        if anchor:
            print(f"\t{anchor2_res[pos]:.4f}", end="")
        print()

# Run the prediction with the specified sequence, mode, and anchor option
run_iupred(input_sequence, mode=prediction_mode, anchor=enable_anchor)
