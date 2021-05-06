import numpy as np
import pandas as pd
import struct
import argparse
import os
import glob
import pickle
import csv

# Accepts in input stream to read from, which should point to a binary config file
def unpack_config_binary(infile):
    out = []
    size = struct.unpack('I', infile.read(4))[0]
    for idx in range(size):
        out.append(struct.unpack('I', infile.read(4))[0])
    return out

# Accepts an input stream to read from, which should point to a binary times file
def unpack_times_binary(infile):
    out = []
    size = struct.unpack('I', infile.read(4))[0]
    for idx in range(size):
        out.append(struct.unpack('d', infile.read(8))[0])
    return out

# Write each set of values per iteration into a row of a csv
def save_to_csv(out, opath, fname):
    with open(opath + fname + '.csv', 'w', newline='') as outfile:
        wr = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        wr.writerow(out)

def collect_file_names(stamp, path):
    config_fnames = []
    times_fnames = []

    os.chdir(path)
    for fname in glob.glob('*' + stamp + '.bin'):
        if "config" in fname:
            config_fnames.append(fname)
        elif "times" in fname:
            times_fnames.append(fname)

    return config_fnames, times_fnames

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GPU Decoder SSA - Config Convert Script")
    parser.add_argument('--ipath', type=str, help="Path to the folder where the binaries are stored.")
    parser.add_argument('--opath', type=str, help="Path to the folder where the converted files will be stored.")
    parser.add_argument('--stamp', type=str, help="Only configs with this timestamp will be converted. Otherwise all configs are converted.")
    args = parser.parse_args()

    stamp = args.stamp
    ipath = args.ipath
    opath = args.opath

    if opath is None:
        opath = ipath

    if stamp is None:
        stamp = ""

    config_fnames, times_fnames = collect_file_names(stamp, ipath)

    print("Converting config binaries.")

    for fname in config_fnames:
        with open(ipath + fname, "rb") as infile:
            save = unpack_config_binary(infile)
            save_to_csv(save, opath, os.path.splitext(fname)[0])

    print("Converting times binaries.")
    for fname in times_fnames:
        with open(ipath + fname, "rb") as infile:
            save = unpack_times_binary(infile)
            save_to_csv(save, opath, os.path.splitext(fname)[0])

    print("All done.")