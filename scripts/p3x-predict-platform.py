#!/usr/bin/env python
"""
Predicts the sequencing platform for a fastq file.

:Authors:
    Jacob S. Porter <jsporter@virginia.edu>
"""
import argparse
import datetime
import os
import statistics
import sys

from numpy import asarray
from platform_features import get_features
from simple_estimator import load_model, predict

MODEL_PATH = os.getenv("PLATFORM_MODEL")
if not MODEL_PATH:
    MODEL_PATH = "./models/reduced/RandomForestClassifier/"


def predict_platform(fastq_input, positions=(1, 1000), model_path=MODEL_PATH):
    """
    Predict the sequencing platform from a fastq file.

    Parameters
    ----------
    fastq_input: str
        The location of the fastq file.
    positions: (int, int)
        The beginning and ending record positions of reads
        to be used for prediction.
    model_path: str
        The location of the Simple_Estimator trained model directory.

    Returns
    -------
    platform: str
        A string representing the platform.
    prob_1: float
        The probability that the prediction is correct.
        This is the frequency of reads predicted that
        resulted in the most common platform.
    prob_2: float
        An alternate probability of prediction correctness.
        This is the average probability of the reads as determined
        by the classifier.
    count: int
        The number of feature records created.

    """
    count, features, _, _ = get_features(
        fastq_input, label=None, positions=positions, output=None
    )
    features = asarray(features)
    model, encoder, name = load_model(model_path)
    responses, proba, order = predict(model, name, features, encoder)
    responses = responses.tolist()
    platform = max(set(encoder.classes_), key=responses.count)
    order_p = 0
    for i, item in enumerate(order):
        if item == platform:
            order_p = i
            break
    preds = list(filter(lambda x: x[0] == platform, zip(responses, proba)))
    prob_1 = len(preds) / len(responses)
    prob_2 = 0.0
    stdev = []
    for item in preds:
        prob_2 += item[1][order_p]
        stdev.append(item[1][order_p])
    prob_2 /= len(preds)
    return fastq_input, platform, prob_1, prob_2, statistics.stdev(stdev), count


def nonnegative(value):
    """Check that value is a non-negative integer."""
    my_error = argparse.ArgumentTypeError(
        "%s is not a non-negative integer value" % value
    )
    try:
        my_value = int(value)
    except ValueError as non_neg:
        raise my_error from non_neg
    if my_value < 0:
        raise my_error
    return my_value


def main():
    """Parse the arguments."""
    tick = datetime.datetime.now()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__
    )
    parser.add_argument(
        "fastq_input", type=str, nargs="+", help=("An input fastq file.")
    )
    parser.add_argument(
        "--range",
        "-r",
        type=nonnegative,
        nargs=2,
        help=(
            "The beginning and ending positions of reads to evaluate.  "
            "Use '0 0' to specify the whole file."
        ),
        default=[1, 1000],
    )
    parser.add_argument(
        "--model",
        "-m",
        help=("The path to the Simple_Estimator model to use for prediction."),
        default=MODEL_PATH,
    )
    args = parser.parse_args()
    print("Predicting platform...", file=sys.stderr)
    print("Started at: {}".format(tick), file=sys.stderr)
    print(args, file=sys.stderr)
    print(
        "\t".join(
            ("File", "Platform", "Maj. Vote", "Avg. Prob.", "Stdev Prob.", "Reads")
        )
    )
    for fastq_file in args.fastq_input:
        print("\t".join(map(str, predict_platform(fastq_file, args.range, args.model))))
    tock = datetime.datetime.now()
    print("The process took time: {}".format(tock - tick), file=sys.stderr)


if __name__ == "__main__":
    main()
