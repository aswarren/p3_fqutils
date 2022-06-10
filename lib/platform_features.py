#!/usr/bin/env python
"""
Extract features from the phred quality string.

:Authors:
    Jacob Porter <jsporter@virginia.edu>
"""

import argparse
import datetime
import math
import sys
from collections import Counter

from SeqIterator import SeqReader

# String of all of the quality characters.
QUAL_STR = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

# List of features.
HEADER = [
    "mean", "max", "min", "variance", "skewness", "kurtosis", "mean_diff",
    "interval_diff"
]


def transform_phred_to_prob(c, offset=33):
    """Transform a quality character into a probability."""
    return 10**((ord(c) - offset) / (-10.0))


def quality_features(vector, reduced=False):
    """Extract features from a probability vector."""
    # sys.stderr.write(" ".join(map(str, vector)) + "\n")
    mean = float(sum(vector)) / len(vector)
    my_max = max(vector)
    my_min = min(vector)
    variance = float(sum([math.pow(item - mean, 2)
                          for item in vector])) / len(vector)
    # Division by zero error
    if reduced:
        return [mean, my_max, my_min, variance]
    try:
        skewness = float(
            sum([
                math.pow((item - mean) / (variance + 0.0), 3)
                for item in vector
            ])) / len(vector)
    except ZeroDivisionError:
        skewness = 0
    try:
        kurt = float(
            sum([
                math.pow((item - mean) / (variance + 0.0), 4)
                for item in vector
            ])) / len(vector)
    except ZeroDivisionError:
        kurt = 9 / 5.0
    diff = []
    for i in range(1, len(vector)):
        diff.append(vector[i] - vector[i - 1])
    mean_diff = float(sum(diff)) / len(diff)
    return [
        mean, my_max, my_min, variance, skewness, kurt, mean_diff,
        vector[0] - vector[len(vector) - 1]
    ]


def get_offset(qual):
    """Get the quality offset from the quality character string."""
    qual_counter = Counter(qual)
    count_list = [qual_counter[char] for char in QUAL_STR]
    offset_33 = sum(count_list[0:25])
    offset_64 = sum(count_list[42:72])
    # offset_inb = sum(count_list[25:42])
    if offset_64 == 0 and offset_33 == 0:
        return 64
    elif offset_33 == 0:
        return 64
    return 33


def get_features(fastq_input,
                 label=None,
                 subportions=3,
                 positions=(0, 0),
                 header=False,
                 reduced=True,
                 output=sys.stdout,
                 debug=False):
    """
    Get the features from sequences in the fastq file.

    Parameters
    ----------
    fastq_input: str
        The location of a fastq file.
    label: str or None
        The platform of the fastq file.
        Use None if the platform is unknown.  For example, for prediction.
    subportions: int
        The number of subportions to extract features.
    positions: (int, int)
        The beginning and ending of positions of fastq sequences.
        Use (0, 0) to get the whole file.
    header: bool
        If True, a header labeling the features will be given.
        If False, no header will be given.
    reduced: bool
        If True, use a reduced feature set. (recommended)
        If False, use all features.  This could result in bad numbers.
    output: writeable or None
        If a writeable is given, write the output to this 'file'.
        If None is given, the features will be returned in a list.
    debug: bool
        Controls additional debug information.

    Returns
    -------
    count: int
        An integer of the number of sequences processed.
    all_features: lst
        A list of features.  (output must be None)
    all_labels: lst
        A list of labels.  (output must be None)
    my_header: lst
        A list of strings labeling each column.  (output must be None.)

    """
    def get_qual_features(qual_ascii, reduced=False):
        seq_qual_prob = list(
            map(lambda x: transform_phred_to_prob(x, offset=offset),
                qual_ascii))
        return quality_features(seq_qual_prob, reduced=reduced)

    reader = SeqReader(fastq_input, file_type='fastq')
    position = 0
    x, y = positions[0], positions[1]
    if x and y and y < x:
        x, y = y, x
    count = 0
    my_header = []
    if header:
        if reduced:
            my_header = HEADER[0:4]
        else:
            my_header = HEADER.copy()
        for i in range(subportions + subportions - 1):
            for item in my_header:
                my_header.append(item + "_" + str(i + 1))
        if label:
            my_header.append("label")
        if output:
            print("\t".join(my_header), file=sys.stdout)
    all_features = []
    all_labels = []
    for record in reader:
        features = []
        _, read, qual, _ = record
        position += 1
        if len(qual) == 0 or (x and (position < x or position > y)):
            continue
        count += 1
        offset = get_offset(qual)
        if debug:
            print("{} Qual len:{} Offset: {}".format(count, len(qual), offset),
                  file=sys.stderr)
        features += get_qual_features(qual, reduced=reduced)
        totallength = int(len(read) / subportions)
        halflength = int(totallength / 2)
        for i in range(subportions + subportions - 1):
            if i == subportions + subportions - 2:
                finallength = len(read)
            else:
                finallength = i * halflength + totallength
            if debug:
                sys.stderr.write(
                    "{} Read len:{} Offset: {} Begin: {}, End: {} {} \n".
                    format(count, len(read), offset, i * halflength,
                           finallength, qual[i * halflength:finallength]))
                sys.stderr.flush()
            features += get_qual_features(qual[i * halflength:finallength],
                                          reduced)
        if label and output:
            features.append(label)
        elif label and not output:
            all_labels.append(label)
        if output:
            print("\t".join(map(str, features)), file=output)
        else:
            all_features.append(features)
    if not output:
        return count, all_features, all_labels, my_header
    return count


def nonnegative(value):
    """Check that value is a non-negative integer."""
    my_error = argparse.ArgumentTypeError(
        "%s is not a non-negative integer value" % value)
    try:
        my_value = int(value)
    except ValueError:
        raise my_error
    if my_value < 0:
        raise my_error
    return my_value


def main():
    """Parse the arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__)
    parser.add_argument("fastq_input", type=str, help=('An input fastq file.'))
    parser.add_argument("label",
                        type=str,
                        help=("The string to label the rows."))
    parser.add_argument("--subportions",
                        "-s",
                        type=int,
                        help=('The number of subportions to divide into.'),
                        default=3)
    parser.add_argument("--range",
                        "-r",
                        nargs=2,
                        type=nonnegative,
                        help=('The range of reads to sample.  '
                              'To process the whole file, use "0 0".'),
                        default=(0, 0))
    parser.add_argument(
        "--header",
        "-d",
        action="store_true",
        help=("Print a header at the top of the feature files."),
        default=False)
    parser.add_argument("--reduced",
                        "-r",
                        action="store_true",
                        help=("Reduce the number of features."),
                        default=False)
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help=("The location of the file to write the output to."),
        default="stdout")
    parser.add_argument("--debug",
                        action="store_true",
                        help=("Print more info to stderr."),
                        default=False)
    args = parser.parse_args()
    print("Extracting features...", file=sys.stderr)
    print("Started at: {}".format(tick), file=sys.stderr)
    print(args, file=sys.stderr)
    if args.output == "stdout":
        output = sys.stdout
    else:
        output = open(args.output, "w")
    count = get_features(args.fastq_input,
                         args.label,
                         args.subportions,
                         positions=args.range,
                         header=args.header,
                         reduced=args.reduced,
                         output=output,
                         debug=args.debug)
    tock = datetime.datetime.now()
    print("There were {} records processed.".format(count), file=sys.stderr)
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == "__main__":
    main()
